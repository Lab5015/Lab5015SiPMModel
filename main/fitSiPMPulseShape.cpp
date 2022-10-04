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
TVirtualFitter* fitter;
std::map<int,int> parIndex;
std::vector<double> fitPars;

std::vector<SiPMParams> SiPMParamsVec;
std::vector<int> runsVec;

// Define time axis and binning
const int npt     = 1024;
const double tmax = 200.;

// TGraphs and TProfile for data
unsigned int nRuns = -1;
std::vector<TGraphErrors*> g_data;
std::vector<TH1D*> hSiPM;
std::vector<TH1D*> hLaser;
std::vector<TH1D*> hInTot;
std::vector<TH1D*> hBand;
std::vector<TH1D*> hOuTot;



void myLoadData(std::vector<int> runs, const std::string& inFolder)
{
  for(unsigned int iRun = 0; iRun < nRuns; ++iRun)
  {
    TFile* inFile = TFile::Open(Form("%s/pulses_run%04d.root",inFolder.c_str(),runsVec[iRun]),"READ");
    TProfile* p_data = (TProfile*)inFile->Get("p_ps_noiseFilter_avg");
    
    double x,y, ey;
    for( int iBin = 1; iBin <= p_data->GetNbinsX(); ++iBin)
    {
      x = p_data -> GetBinCenter(iBin); 
      y = p_data -> GetBinContent(iBin);
      ey = p_data -> GetBinError(iBin);
      g_data[iRun] -> SetPoint(g_data[iRun]->GetN(),x,y);
      g_data[iRun] -> SetPointError(g_data[iRun]->GetN()-1,0.,ey);
    }
    
    inFile->Close();
  }
  
  return;
}



// ******************************
// SiPM model on an external load
// ******************************
void myInSignal(const int& iRun, double* par)
{
  SiPMParams* sipmPars = &(SiPMParamsVec.at(iRun));

  double amp = 5.;
  double t0  = 0.;
  double Rq  = (*sipmPars).Rq;
  double BW  = 60.;
    
  if( parIndex[0+4*iRun] >= 0 ) amp = par[parIndex[0+4*iRun]];
  if( parIndex[1+4*iRun] >= 0 )  t0 = par[parIndex[1+4*iRun]];
  if( parIndex[2+4*iRun] >= 0 )  Rq = par[parIndex[2+4*iRun]];
  if( parIndex[3+4*iRun] >= 0 )  BW = par[parIndex[3+4*iRun]];

  (*sipmPars).Rq = Rq;
  
  fitPars[0+4*iRun] = amp;
  fitPars[1+4*iRun] = t0;
  fitPars[2+4*iRun] = Rq;
  fitPars[3+4*iRun] = BW;
  
  
  // bandwidth low-pass filter
  double tau_band = 0.35/(BW*2.2)*1e3; // nanosenconds (BW is in GHz)
  
  // Hitograms for FFT analysis 
  for(int ipt = 0; ipt <= npt; ++ipt)
  { 
    double tt = (double)ipt * (tmax/npt);
    
    hSiPM[iRun]  -> SetBinContent(ipt+1, (*sipmPars).RL*SiPMPulseShape(tt,*sipmPars,amp,t0));
    hBand[iRun]  -> SetBinContent(ipt+1, funcLP(tt,tau_band));
  }
  
  hConvol(hSiPM[iRun],hLaser[iRun],hInTot[iRun]);
  hConvol(hInTot[iRun],hBand[iRun],hOuTot[iRun]);
  
  return;
} 



// ******************************
// Draw result
// ******************************
void myDrawFit(const std::string& plotFolder, const std::string& commonLabel)
{
  for(int iRun = 0; iRun < nRuns; ++iRun)
  {
    SiPMParams sipmPars = SiPMParamsVec.at(iRun);
    
    // find maximum for plot range
    double Amax = 0.;
    double xx, yy;
    for(int iPoint = 0; iPoint < g_data[iRun]->GetN(); ++iPoint)
    {
      g_data[iRun] -> GetPoint(iPoint,xx,yy);
      if( Amax <  yy ) 
        Amax = yy;
    }
    for(int iBin = 1; iBin <= hOuTot[iRun]->GetNbinsX(); ++iBin)
    {
      if( Amax <  hOuTot[iRun]->GetBinContent(iBin) ) 
        Amax = hOuTot[iRun]->GetBinContent(iBin);
    } 
    
    // Draw
    TCanvas* c1 = new TCanvas(Form("c%d",iRun), Form("c%d",iRun), 1800., 750.); 
    c1->Divide(2,1);
    c1 -> cd(1);
    gPad->SetLeftMargin(0.2); gPad->SetRightMargin(0.1);
    gPad->SetTicks();
    
    TH1F* hFrame = (TH1F*)gPad->DrawFrame(0.,-0.01,200.,1.2*Amax);
    hFrame -> SetTitle(Form("%s",sipmPars.title.c_str())); 
    hFrame -> GetXaxis()->SetTitle("time [ns]");
    hFrame -> GetYaxis()->SetTitle("Amplitude [V]"); 
    hFrame -> GetXaxis()->SetLabelSize(0.045); 
    hFrame -> GetYaxis()->SetLabelSize(0.045); 
    hFrame -> GetXaxis()->SetTitleSize(0.050); 
    hFrame -> GetYaxis()->SetTitleSize(0.050);
    
    g_data[iRun] -> SetMarkerSize(0.5);
    g_data[iRun] -> SetMarkerStyle(8);
    g_data[iRun] -> Draw("P, same");
    
    hOuTot[iRun] ->SetLineColor(kGreen+1); hOuTot[iRun] ->SetLineWidth(2);
    hOuTot[iRun] ->Draw("L,same");
    
    TLegend* legend = new TLegend(0.40,0.75,0.70,0.85);
    legend->SetBorderSize(0); 
    legend->SetTextSize(0.04);
    legend->AddEntry(hOuTot[iRun],"Complete model","L");
    legend->AddEntry(g_data[iRun],"Data","LP"); 
    legend->Draw("same");
    
    TLatex* latex = new TLatex(0.40,0.60,Form("V_{OV} = %.1f V", sipmPars.OV));
    latex -> SetNDC();
    latex -> SetTextFont(42);
    latex -> SetTextSize(0.04);
    latex -> SetTextColor(kBlack);
    latex -> Draw("same");
    
    TLatex* latex1 = new TLatex(0.40,0.55,Form("N_{c} = %.0f   C_{g} = %.0f pF", sipmPars.Nc, sipmPars.Cg*1e12));
    latex1 -> SetNDC();
    latex1 -> SetTextFont(42);
    latex1 -> SetTextSize(0.04);
    latex1 -> SetTextColor(kBlack);
    latex1 -> Draw("same");
    
    TLatex* latex2 = new TLatex(0.40,0.50,Form("C_{d} = %.1f fF, C_{q} = %.1f fF, R_{q} = %.0f k#Omega", sipmPars.Cd*1e15,sipmPars.Cq*1e15, fitPars[2+4*iRun]*1e-3));
    latex2 -> SetNDC();
    latex2 -> SetTextFont(42);
    latex2 -> SetTextSize(0.04);
    latex2 -> SetTextColor(kBlack);
    latex2 -> Draw("same");
    
    TLatex* latex3 = new TLatex(0.40,0.40,Form("BW = %.0f MHz", fitPars[3+4*iRun]));
    latex3 -> SetNDC();
    latex3 -> SetTextFont(42);
    latex3 -> SetTextSize(0.04);
    latex3 -> SetTextColor(kBlack);
    latex3 -> Draw("same");
    
    c1 -> cd(2);
    gPad->SetLeftMargin(0.2); gPad->SetRightMargin(0.1); 
    gPad->SetTicks();
    gPad->SetLogy(); 
    hFrame ->Draw(); 
    hFrame -> SetMinimum(0.002);
    hFrame -> GetXaxis() -> SetRangeUser(0.,200.);
    g_data[iRun] -> Draw("P,same");
    hOuTot[iRun] -> Draw("L,same");
    gPad->Update();
    c1->Update();
    
    c1 -> Print(Form("%s/c_%s_%s.png",plotFolder.c_str(),sipmPars.label.c_str(),commonLabel.c_str()));
    
    TFile outFile(Form("%s/plots_%s_%s.root",plotFolder.c_str(),sipmPars.label.c_str(),commonLabel.c_str()), "RECREATE");
    outFile.cd();
    g_data[iRun] -> Write("g_data");
    hSiPM[iRun] -> Write();
    hLaser[iRun] -> Write();
    hInTot[iRun] -> Write();
    hBand[iRun] -> Write();
    hOuTot[iRun] -> Write();
  }
  
  return;
}



// ************************
// FCN 
// ************************
void fcn(int& npar, double* gin, double& f, double* par, int iflag)
{
  double chisq = 0.;
  double delta = 0.;
  
  for(unsigned int iRun = 0; iRun < nRuns; ++iRun)
  {
    myInSignal(iRun, par);

    // restrict fit range
    double x1 = hOuTot[iRun]->GetXaxis()->GetXmin();
    double x2 = hOuTot[iRun]->GetXaxis()->GetXmax();
    x1 = 10.; 
    x2 = 180.;

    double xx, yy, ey;
    for(int iPoint = 0; iPoint < g_data[iRun]->GetN(); ++iPoint)
    {
      g_data[iRun] -> GetPoint(iPoint,xx,yy);
      if( yy > 0.002 )
      {
        x1 = xx;
        break;
      }
    }
    for(int iPoint = g_data[iRun]->GetN()-1; iPoint >= 0; --iPoint)
    {
      g_data[iRun] -> GetPoint(iPoint,xx,yy);
      if( yy > 0.002 )
      {
        x2 = xx;
        break;
      }
    }

    // compute chisquare
    for(int iPoint = 0; iPoint < g_data[iRun]->GetN(); iPoint++)
    {
      g_data[iRun] -> GetPoint(iPoint,xx,yy);
      ey = g_data[iRun]->GetErrorY(iPoint);
      if( xx > x1 && xx < x2)
      { 
        delta = ( yy - hOuTot[iRun]->Interpolate(xx) ) / ey;
        chisq += delta*delta;
        //std::cout << "iPoint: " << iPoint << "   xx: " << xx << "  yy: " << yy << "   func: " << hOuTot[iRun]->Interpolate(xx) << "   delta: " << delta << "   chisq: " << chisq << std::endl;
      }
    }
  }
  
  f = chisq;
}



// ************************
// fit the pulse shape
// ************************
void myFitPulse(const bool& minimize,
                const int& nPars_amp, const int& nPars_t0, const int& nPars_Rq, const int& nPars_BW)
{
  int nPars = nPars_amp+nPars_t0+nPars_Rq+nPars_BW;
  
  // Define fit model
  TVirtualFitter::SetDefaultFitter("Minuit");
  fitter = TVirtualFitter::Fitter(NULL,nPars);
  fitter->SetFCN(fcn);

  double arglist[100];
  
  // SET PRINT LEVEL
  arglist[0] = 1;
  fitter -> ExecuteCommand("SET PRINT",arglist,1);
  
  // ERROR CHISQUARE
  arglist[0] = 1;
  fitter -> ExecuteCommand("SET ERR", arglist,1);
  
  // INITIALIZE PARAMETERS
  for(int iPar = 0; iPar < nPars_amp; ++iPar) fitter -> SetParameter(iPar, Form("Amp"), 5., 0.2, 0.1,100.);
  for(int iPar = 0; iPar < nPars_t0; ++iPar)  fitter -> SetParameter(nPars_amp+iPar, Form("t0"), 0., 0.1, -10., 10.);
  for(int iPar = 0; iPar < nPars_Rq; ++iPar)  fitter -> SetParameter(nPars_amp+nPars_t0+iPar, Form("Rq"), 450e3, 10e3, 0. ,1000e3);
  for(int iPar = 0; iPar < nPars_BW; ++iPar)  fitter -> SetParameter(nPars_amp+nPars_t0+nPars_Rq+iPar, "BW", 60., 2., 0., 200.);
  
  // INITIALIZE FCN
  arglist[0] = 1; 
  fitter -> ExecuteCommand("CALL FCN", arglist, 1);
  
  // MINIMIZATION
  if( minimize )
  {
    arglist[0] = 5000; // number of function calls
    arglist[1] = 0.05; // tolerance
    fitter -> ExecuteCommand("MINIMIZE",arglist,2);
  }
  
  return;
}






int main(int argc, char** argv)
{
  //----------------------
  // parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  
  int doMinimization = opts.GetOpt<int>("Input.doMinimization");
  std::string dataFolder = opts.GetOpt<std::string>("Input.dataFolder");
  std::string commonLabel = opts.GetOpt<std::string>("Input.commonLabel");
  std::string plotFolder = opts.GetOpt<std::string>("Output.plotFolder");
  
  
  //--------------------------  
  // get the SiPM and run list
  GetSiPMParsFromCfg(argv[1],SiPMParamsVec,runsVec);
  nRuns = runsVec.size();
  fitPars.reserve(4*nRuns);
  

  //-----------------------  
  // get the fit parameters
  int nPars_amp = -1;
  int nPars_t0 = -1;
  int nPars_Rq = -1;
  int nPars_BW = -1;
  GetFitParsFromCfg(argv[1], nRuns, parIndex, nPars_amp, nPars_t0, nPars_Rq, nPars_BW);

  
  //-------------------
  // define histograms
  for(unsigned int iRun = 0; iRun < nRuns; ++iRun)
  {
    g_data.push_back( new TGraphErrors() );
    hSiPM.push_back(  new TH1D(Form("hSiPM_%d",iRun), "",npt,0,tmax) ); // SiPM pulse
    hLaser.push_back( new TH1D(Form("hLaser_%d",iRun),"",npt,0,tmax) ); // laser response
    hInTot.push_back( new TH1D(Form("hInTot_%d",iRun),"",npt,0,tmax) ); // output
    hBand.push_back(  new TH1D(Form("hBand_%d",iRun), "",npt,0,tmax) ); // readout board bandwidth
    hOuTot.push_back( new TH1D(Form("hOuTot_%d",iRun),"",npt,0,tmax) ); // output
  }
  
  
  //-----------------------
  // fill static histograms
  // laser response
  double sigma = 0.021; // 21 ps sigma for our laser
  double mu = 20.; // arbitrary
  for(unsigned int iRun = 0; iRun < nRuns; ++iRun)
  {
    SiPMParams* sipmPars = &(SiPMParamsVec.at(iRun));
    
    for(int ipt = 0; ipt <= npt; ++ipt)
    { 
      double tt = (double)ipt * (tmax/npt);
      if( TMath::Gaus(tt,mu,sigma,true) > 1e-5) hLaser[iRun] -> SetBinContent(ipt+1, TMath::Gaus(tt,mu,sigma,true));
      else                                      hLaser[iRun] -> SetBinContent(ipt+1, 1e-15);
    }
    hLaser[iRun] -> Scale((*sipmPars).Npe/hLaser[iRun]->Integral());
  }
  
  
  //----------------
  // run the analsys
  myLoadData(runsVec,dataFolder); // fill graphs from files
  
  myFitPulse(bool(doMinimization),
             nPars_amp, nPars_t0, nPars_Rq, nPars_BW);
  
  myDrawFit(plotFolder,commonLabel); // draw the fit result
  
  return 0;
}
