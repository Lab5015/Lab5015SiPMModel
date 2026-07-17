#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include "interface/Convolution.h"
#include "interface/Functions.h"

//#****************************************************************
// Fit SiPM and LYSO Data with a model to filter and clip the pulse
//*****************************************************************

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

float BW = -1.;
float RF = -1.;
float Npe_LYSO = -1.;
float Npe_laser = -1.;

// Define time axis and binning
const int npt     = int(pow(2,16));
const double tmax = 1000.;

// TGraphs and TProfile for data
std::vector<TH1D*> hSiPM;
std::vector<TH1D*> hSiPM_1;
std::vector<TH1D*> hSiPM_2;
std::vector<TH1D*> hSiPM_3;
std::vector<TH1D*> hSiPM_4;
std::vector<TH1D*> hSiPM_laser;
std::vector<TH1D*> hSiPM_LYSO;
std::vector<TH1D*> hAdvanSid_laser;
std::vector<TH1D*> hAdvanSid_LYSO;
std::vector<TH1D*> hTOFHIR_laser;
std::vector<TH1D*> hTOFHIR_LYSO;
TH1D* hLYSO;
TH1D* hLaser;
TH1D* hBandAdvanSid;
TH1D* hBandTOFHIR;






// ******************************
// SiPM model on an external load
// ******************************
void myInSignal()
{
  // define histograms
  hLYSO = new TH1D(Form("hLYSO"),"",npt,0,tmax); // LYSO scintillation
  hLaser = new TH1D(Form("hLaser"),"",npt,0,tmax); // Laser scintillation
  hBandAdvanSid = new TH1D(Form("hBandAdvanSid"),"",npt,0,tmax); // AdvanSid bandwidth
  
  for(unsigned int SiPMIt = 0; SiPMIt < SiPMParamsVec.size(); ++SiPMIt)
    {
      SiPMParams* sipmPars = &(SiPMParamsVec.at(SiPMIt));

      hSiPM.push_back( new TH1D(Form("hSiPM_%s",(*sipmPars).label.c_str()), "",npt,0,tmax) );
      hSiPM_1.push_back( new TH1D(Form("hSiPM_1_%s",(*sipmPars).label.c_str()), "",npt,0,tmax) );
      hSiPM_2.push_back( new TH1D(Form("hSiPM_2_%s",(*sipmPars).label.c_str()), "",npt,0,tmax) );
      hSiPM_3.push_back( new TH1D(Form("hSiPM_3_%s",(*sipmPars).label.c_str()), "",npt,0,tmax) );
      hSiPM_4.push_back( new TH1D(Form("hSiPM_4_%s",(*sipmPars).label.c_str()), "",npt,0,tmax) );
      hSiPM_laser.push_back( new TH1D(Form("hSiPM_laser_%s",(*sipmPars).label.c_str()), "",npt,0,tmax)  );
      hSiPM_LYSO.push_back( new TH1D(Form("hSiPM_LYSO_%s",(*sipmPars).label.c_str()), "",npt,0,tmax)  );
      hAdvanSid_laser.push_back( new TH1D(Form("hAdvanSid_laser%s",(*sipmPars).label.c_str()), "",npt,0,tmax)  );
      hAdvanSid_LYSO.push_back( new TH1D(Form("hAdvanSid_LYSO_%s",(*sipmPars).label.c_str()), "",npt,0,tmax)  );
      hTOFHIR_laser.push_back( new TH1D(Form("hTOFHIR_laser_%s",(*sipmPars).label.c_str()), "",npt,0,tmax)  );
      hTOFHIR_LYSO.push_back( new TH1D(Form("hTOFHIR_LYSO_%s",(*sipmPars).label.c_str()), "",npt,0,tmax)  );      
    }
  
  // fill bandwidth histograms
  // laser response
  double sigma = 0.021; // 21 ps sigma for our laser
  double mu = 20.; // arbitrary
  // LYSO scintillation
  double tau_r = 0.1; // ns
  double tau_d = 38.5; // ns
  // AdvanSid bandwidth low-pass filter
  double tau_band = 0.35/(BW*2.2)*1e3; // nanosenconds (BW is in GHz)
  // TOFHIR
  //...
  
  for(int ipt = 0; ipt <= npt; ++ipt)
    { 
      double tt = (double)ipt * (tmax/npt);
      
      hLYSO -> SetBinContent(ipt+1, funcScint(tt,tau_r,tau_d));
      
      if( TMath::Gaus(tt,mu,sigma,true) > 1e-5) hLaser -> SetBinContent(ipt+1, TMath::Gaus(tt,mu,sigma,true));
      else                                      hLaser -> SetBinContent(ipt+1, 1e-10);
      
      hBandAdvanSid -> SetBinContent(ipt+1, funcLP(tt,tau_band));
      
      // TOFHIR...
    }
  
  hLYSO -> Scale(Npe_LYSO/hLYSO->Integral());
  hLaser -> Scale(Npe_laser/hLaser->Integral());
  
  
  // Hitograms for FFT analysis 
  for(unsigned int SiPMIt = 0; SiPMIt < SiPMParamsVec.size(); ++SiPMIt)
    {
      SiPMParams* sipmPars = &(SiPMParamsVec.at(SiPMIt));
      
      for(int ipt = 0; ipt <= npt; ++ipt)
	{ 
	  double tt = (double)ipt * (tmax/npt);
	  
	  hSiPM[SiPMIt] -> SetBinContent(ipt+1, 1e6*SiPMPulseShape(tt,*sipmPars,1.,10.));
	  hSiPM_1[SiPMIt] -> SetBinContent(ipt+1, 1e6*SiPMPulseShape(tt,*sipmPars,1.,10.,1));
	  hSiPM_2[SiPMIt] -> SetBinContent(ipt+1, 1e6*SiPMPulseShape(tt,*sipmPars,1.,10.,2));
	  hSiPM_3[SiPMIt] -> SetBinContent(ipt+1, 1e6*SiPMPulseShape(tt,*sipmPars,1.,10.,3));
	  hSiPM_4[SiPMIt] -> SetBinContent(ipt+1, 1e6*SiPMPulseShape(tt,*sipmPars,1.,10.,4));
	  hSiPM_laser[SiPMIt] -> SetBinContent(ipt+1, 1e6*SiPMPulseShape(tt,*sipmPars,1.,10.));
	  hSiPM_LYSO[SiPMIt] -> SetBinContent(ipt+1, 1e6*SiPMPulseShape(tt,*sipmPars,1.,10.));
	  
	  hAdvanSid_laser[SiPMIt] -> SetBinContent(ipt+1, 0.5*RF*SiPMPulseShape(tt,*sipmPars,1.,10.));
	  hAdvanSid_LYSO[SiPMIt] -> SetBinContent(ipt+1, 0.5*RF*SiPMPulseShape(tt,*sipmPars,1.,10.));
	  
	  // TOFHIR...
	}
      
      // convolutions
      hConvol(hSiPM[SiPMIt],hLaser,hSiPM_laser[SiPMIt]);
      hConvol(hSiPM[SiPMIt],hLYSO,hSiPM_LYSO[SiPMIt]);
      
      hConvol(hAdvanSid_laser[SiPMIt],hLaser,hAdvanSid_laser[SiPMIt]);
      hConvol(hAdvanSid_laser[SiPMIt],hBandAdvanSid,hAdvanSid_laser[SiPMIt]);
      hConvol(hAdvanSid_LYSO[SiPMIt],hLYSO,hAdvanSid_LYSO[SiPMIt]);
      hConvol(hAdvanSid_LYSO[SiPMIt],hBandAdvanSid,hAdvanSid_LYSO[SiPMIt]);
    }
  
  return;
} 



// ******************************
// Draw result
// ******************************
void myDrawFit(std::vector<TH1D*> histos, const std::string& label, const std::string& title,
	       const std::string& plotFolder, const std::string& commonLabel)
{
  system(Form("mkdir -p %s",plotFolder.c_str())); 
  
  
  TCanvas* c1 = new TCanvas(Form("c_%s",label.c_str()), Form("c_%s",label.c_str()), 1800., 750.); 
  c1 -> Divide(2,1);
  c1 -> cd(1);
  gPad -> SetLeftMargin(0.2); gPad -> SetRightMargin(0.1);
  gPad -> SetTicks();
  
  float Amax = -1.;
  for(unsigned int SiPMIt = 0; SiPMIt < SiPMParamsVec.size(); ++SiPMIt)
    {
      for(int iBin = 1; iBin <= histos[SiPMIt]->GetNbinsX(); ++iBin)
	{
	  if( Amax <  histos[SiPMIt]->GetBinContent(iBin) ) 
	    Amax = histos[SiPMIt]->GetBinContent(iBin);
	} 
    }
  
  TH1F* hFrame = (TH1F*)gPad->DrawFrame(0.,0.,200.,1.2*Amax);
  hFrame -> SetTitle(title.c_str());
  hFrame -> GetXaxis()->SetLabelSize(0.045); 
  hFrame -> GetYaxis()->SetLabelSize(0.045); 
  hFrame -> GetXaxis()->SetTitleSize(0.050); 
  hFrame -> GetYaxis()->SetTitleSize(0.050);
  
  TLegend* legend = new TLegend(0.22,0.80,0.70,0.80+0.04*SiPMParamsVec.size());
  legend -> SetBorderSize(0); 
  legend -> SetTextSize(0.04);
  
  for(unsigned int SiPMIt = 0; SiPMIt < SiPMParamsVec.size(); ++SiPMIt)
    {
      SiPMParams* sipmPars = &(SiPMParamsVec.at(SiPMIt));
      histos[SiPMIt] -> SetLineColor(1+SiPMIt);
      histos[SiPMIt] -> SetLineWidth(3);
      histos[SiPMIt] -> Draw("same");
      legend -> AddEntry(histos[SiPMIt],Form("%s",(*sipmPars).title.c_str()),"LP"); 
    }
  legend -> Draw("same");
  
  c1 -> cd(2);
  gPad -> SetLeftMargin(0.2); gPad -> SetRightMargin(0.1); 
  gPad -> SetTicks();
  gPad -> SetLogy(); 
  hFrame -> Draw(); 
  hFrame -> SetMinimum(Amax/300.);

  for(unsigned int SiPMIt = 0; SiPMIt < SiPMParamsVec.size(); ++SiPMIt)
    {
      histos[SiPMIt] -> SetLineColor(1+SiPMIt);
      histos[SiPMIt] -> SetLineWidth(3);
      histos[SiPMIt] -> Draw("same");
    }
  
  c1 -> Update();
  c1 -> Print(Form("%s/c_%s_%s.png",plotFolder.c_str(),label.c_str(),commonLabel.c_str()));
  
  return;
}






int main(int argc, char** argv)
{
  //----------------------
  // parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  
  std::string commonLabel = opts.GetOpt<std::string>("Input.commonLabel");
  std::string plotFolder = opts.GetOpt<std::string>("Output.plotFolder");
  
  BW = opts.GetOpt<float>("Input.BW");
  RF = opts.GetOpt<float>("Input.RF");
  Npe_laser = opts.GetOpt<float>("Input.Npe_laser");
  Npe_LYSO = opts.GetOpt<float>("Input.Npe_LYSO");  
  
  //------------------
  // get the SiPM list
  GetSiPMParsFromCfg(argv[1],SiPMParamsVec);
  
  
  //----------------
  // run the analsys
  myInSignal();
  myDrawFit(hSiPM,"SiPM",";time [ns];current [#muA]",plotFolder,commonLabel); // draw the fit result
  myDrawFit(hSiPM_laser,"SiPM_laser",";time [ns];current [#muA]",plotFolder,commonLabel); // draw the fit result
  myDrawFit(hSiPM_LYSO,"SiPM_LYSO",";time [ns];current [#muA]",plotFolder,commonLabel); // draw the fit result
  myDrawFit(hAdvanSid_laser,"AdvanSid_laser",";time [ns];voltage [V]",plotFolder,commonLabel); // draw the fit result
  myDrawFit(hAdvanSid_LYSO,"AdvanSid_LYSO",";time [ns];voltage [V]",plotFolder,commonLabel); // draw the fit result
  
  
  TFile outFile(Form("%s/plots_%s.root",plotFolder.c_str(),commonLabel.c_str()), "RECREATE");
  outFile.cd();
  for(unsigned int it = 0; it < hSiPM.size(); ++it)
    {
      hSiPM[it] -> Write();
      hSiPM_1[it] -> Write();
      hSiPM_2[it] -> Write();
      hSiPM_3[it] -> Write();
      hSiPM_4[it] -> Write();
      hSiPM_laser[it] -> Write();
      hSiPM_LYSO[it] -> Write();
      hAdvanSid_laser[it] -> Write();
      hAdvanSid_LYSO[it] -> Write();
    }
  
  return 0;
}
