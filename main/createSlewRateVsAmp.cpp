#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include "interface/Convolution.h"
#include "interface/Functions.h"

#include <TROOT.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>
#include <TMath.h>
#include <TLegend.h>
#include <TFile.h>
#include <TLine.h>
#include <TBox.h>
#include <TLine.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TVirtualFFT.h>
#include <TText.h>
#include <TMinuit.h>
#include "TProfile.h"
#include "TLatex.h"
#include "TVirtualFitter.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <sstream>
#include <ctime>
#include <map>
#include <algorithm>
#include <math.h>
#include <vector>
#include <fstream>


// Define global variables
std::vector<SiPMParams> SiPMParamsVec;

// Define time axis and binning
const int npt     = int(pow(2,13));
const double tmax = 409.6;
const double freq = npt/tmax;

// TGraphs and TProfile for data
TH1D* hSiPM;
TH1D* hLYSO;
TH1D* hInTot;
TH1D* hOuTot;;
TH1D* hLowP;
TH1D* highP;
TH1D* hDisc;
TH1D* hDeLP;

TF1* f_pole;



// ******************************
// SiPM model on an external load
// ******************************
void myInSignal(const float& Npe, const float& gain)
{
  // define histograms
  hSiPM  = new TH1D(Form("hSiPM_Npe%.0f",Npe), "",npt,0,tmax); // SiPM model
  hLYSO  = new TH1D(Form("hLYSO_Npe%.0f",Npe), "",npt,0,tmax); // LYSO scintillation
  hInTot = new TH1D(Form("hInTot_Npe%.0f",Npe),"",npt,0,tmax); // TOFHIR input
  hOuTot = new TH1D(Form("hOuTot_Npe%.0f",Npe),"",npt,0,tmax); // TOFHIR response
  hLowP  = new TH1D(Form("hLowP_Npe%.0f",Npe), "",npt,0,tmax); // preamp input low-pass filter
  highP  = new TH1D(Form("highP_Npe%.0f",Npe), "",npt,0,tmax); // high-pass at the pre input
  hDisc  = new TH1D(Form("hDisc_Npe%.0f",Npe), "",npt,0,tmax); // 2^ order low-pass (Disc input)
  hDeLP  = new TH1D(Form("hDeLP_Npe%.0f",Npe), "",npt,0,tmax); // delay line input low-pass filter
  
  // LYSO scintillation
  double tau_r = 0.1; // ns
  double tau_d = 38.5; // ns
  float ampli = Npe * 3.674E-06 / (tmax/npt);
  
  // TOFHIR bandwidth
  float x = Npe/9500. * gain/3.72;   // normalize to 9500 p.e. and a gain of 3.72 (E05), which is ~ the gain at 3.5 Vov
  float tau1 = 2.5;                                   // average - precise for Np = 4700 would be 3.65909;   // TOFHIR LP
  float tau2 = 4/30.;                                 // effect of inductance
  float tau3 = 1. / f_pole -> Eval(x) / 6.28 * 1000.; // high pass
  float tauE = 0.5736;                                // delay line - Elmore pole  
  
  // SiPM parameters
  SiPMParams* sipmPars = &(SiPMParamsVec.at(0));

  for(int ipt = 0; ipt <= npt; ++ipt)
    { 
      double tt = (double)ipt * (tmax/npt);
      
      //hLYSO -> SetBinContent(ipt+1, ampli * funcScint(tt,tau_r,tau_d));
      hLYSO -> SetBinContent(ipt+1, 1./20. *ampli * funcLP(tt,tau_d));   // 1/20 comes from backward FFT normalization in Tommaso's code
      
      hSiPM -> SetBinContent(ipt+1, 1/(20.*20.*20.*20.) * 1./1.602E-10*SiPMPulseShape(tt,*sipmPars,1.,0.));
      
      hLowP -> SetBinContent( ipt+1,  (1.-0.5/(freq*tau1)) * funcLP(tt,tau1));
      hDisc -> SetBinContent( ipt+1,  (1.-0.5/(freq*tau2)) * funcLP(tt,tau2));
      if( ipt == 0 ) highP -> SetBinContent( ipt+1, freq -  (1.-0.5/(freq*tau3)) * funcLP(tt,tau3) );
      if( ipt >  0 ) highP -> SetBinContent( ipt+1, 0    -  (1.-0.5/(freq*tau3)) * funcLP(tt,tau3) );
      if( ipt == 0 ) hDeLP -> SetBinContent( ipt+1, freq -  (1.-0.5/(freq*tauE)) * funcLP(tt,tauE) );
      if( ipt >  0 ) hDeLP -> SetBinContent( ipt+1, 0    -  (1.-0.5/(freq*tauE)) * funcLP(tt,tauE) );
    }
  
  // convolutions
  hConvol(hSiPM,hLYSO,hInTot);
  
  hConvol(hInTot,highP,hOuTot);
  hConvol(hOuTot,hLowP,hOuTot);
  hConvol(hOuTot,hDisc,hOuTot);
  hConvol(hOuTot,hDeLP,hOuTot);
  
  return;
} 






int main(int argc, char** argv)
{
  //----------------------
  // parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  
  std::string label = opts.GetOpt<std::string>("Input.label");
  std::vector<float> Npes = opts.GetOpt<std::vector<float> >("Input.Npes");
  
  std::string plotFolder = opts.GetOpt<std::string>("Output.plotFolder");
  
  TFile outFile(Form("%s/slewRate_%s.root",plotFolder.c_str(),label.c_str()), "RECREATE");
  TGraph* g_SR = new TGraph();
  
  
  //--------------------------------
  // moving pole in TOFHIR bandwidth
  f_pole = new TF1("f_pole","[0]*x+[1]*x*x+[2]*x*x*x",0.,2.);
  f_pole -> SetParameters(-1.67421,34.7462,-15.7621);
  
  
  //------------------
  // get the SiPM list
  GetSiPMParsFromCfg(argv[1],SiPMParamsVec);
  
  
  //----------------
  // run the analsys
  for(auto Npe : Npes)
    {
      float gainScale = (SiPMParamsVec.at(0).Cq+SiPMParamsVec.at(0).Cd) / (12.4e-15+3.2e-15);
      myInSignal(Npe,gainScale*3.72);
      
      //--- get the slew rate
      float x1 = -1.;
      float x2 = -1.;
      float y1 = -1.;
      float y2 = -1.;
      for(int bin = 1; bin <= hOuTot->GetNbinsX(); ++bin)
	{
	  float x = hOuTot->GetBinCenter(bin);
	  float y = hOuTot->GetBinContent(bin);
	  
	  if( y >= 3. && x1 == -1 )
	    {
	      x1 = x;
	      y1 = y;
	      x2 = hOuTot->GetBinCenter(bin+1);
	      y2 = hOuTot->GetBinContent(bin+1);
	      break;
	    }
	}
      
      float x = Npe/9500.*gainScale;
      float SR = (y2-y1)/(x2-x1);
      
      std::cout << "Npe: " << Npe << "   gain scale: " << gainScale << "   x: " << x << "   slew rate: " << SR << std::endl;
      g_SR -> SetPoint(g_SR->GetN(),x,SR);
      
      outFile.cd();
      hSiPM -> Write();
      hLYSO -> Write();
      hInTot -> Write();
      hOuTot -> Write();
   }
  
  g_SR -> Write("g_SR");
  outFile.Close();
  
  return 0;
}
