#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include "interface/Convolution.h"
#include "interface/Functions.h"

/* *****************************************************************
// Fit SiPM and LYSO Data with a model  to filter and clip the pulse
******************************************************************** */

#include <TROOT.h>
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

int nRuns = 20;

TMinuit* gMinuit;

const float OV = 1.5;   // volt

std::vector<SiPMParams> SiPMParamsVec;
std::vector<int> runsVec;

// Define convolutions
const int npt     = 1024;
const double tmax = 200.;

// TGraphs and TProfile for data
TGraphErrors** g_data = new TGraphErrors*[nRuns];   
TH1D** hBand = new TH1D*[nRuns];
TH1D** hSiPM = new TH1D*[nRuns];
TH1D** hOuTot = new TH1D*[nRuns];

double* pars_val;
double* pars_err;
double* pars_low;
double* pars_hig;



void myLoadData(std::vector<int> runs, const std::string& inFolder)
{
  for(unsigned int iRun = 0; iRun < nRuns; ++iRun)
  {
    TFile* inFile = TFile::Open(Form("%s/pulses_run%04d.root",inFolder.c_str(),runsVec[iRun]),"READ");
    TProfile* p_data = (TProfile*)inFile->Get("p_ps_noiseFilter_avg");

    g_data[iRun] = new TGraphErrors();
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
  double amp = par[0+3*iRun];
  double x0  = par[1+3*iRun];
  double Rq  = par[2+3*iRun];
  double BW  = par[3*nRuns];
  (*sipmPars).Rq = Rq;
  
  // bandwidth low-pass filter
  double tau_band = 0.35/(BW*2.2)*1e3; // nanosenconds (BW is in GHz)
  
  // Hitograms for FFT analysis 
  for(int ipt = 0; ipt <= npt; ++ipt)
  { 
    double tt = (double)ipt * (tmax/npt);
    
    hSiPM[iRun] -> SetBinContent(ipt+1, (*sipmPars).RL*SiPMPulseShape(tt,*sipmPars,OV,amp,x0));
    hBand[iRun] -> SetBinContent(ipt+1, funcLP(tt,tau_band));
  }
  
  hConvol(hSiPM[iRun],hBand[iRun],hOuTot[iRun]);
  
  return;
} 



// Draw
void myDrawFit()
{
  double BW, BWerr;
  gMinuit -> GetParameter(3*nRuns,BW,BWerr);
  
  for(int d = 0; d < nRuns; d++)
  {
    SiPMParams sipmPars = SiPMParamsVec.at(d);
    
    // find maximum for plot range
    double Amax = 0.;
    double xx, yy;
    for(int iPoint = 0; iPoint < g_data[d]->GetN(); ++iPoint)
    {
      g_data[d] -> GetPoint(iPoint,xx,yy);
      if( Amax <  yy ) 
        Amax = yy;
    }
    for(int iBin = 1; iBin <= hOuTot[d]->GetNbinsX(); ++iBin)
    {
      if( Amax <  hOuTot[d]->GetBinContent(iBin) ) 
        Amax = hOuTot[d]->GetBinContent(iBin);
    } 
    
    
    // Draw
    TCanvas* c1 = new TCanvas(Form("c%d", d),Form("c%d", d),1800.,750.); 
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
    
    g_data[d] -> SetMarkerSize(0.5);
    g_data[d] -> SetMarkerStyle(8);
    g_data[d] -> Draw("P, same");
    
    hOuTot[d] ->SetLineColor(kGreen+1); hOuTot[d] ->SetLineWidth(2);
    hOuTot[d] ->Draw("L,same");
    
    TLegend* legend = new TLegend(0.40,0.75,0.70,0.85);
    legend->SetBorderSize(0); 
    legend->SetTextSize(0.04);
    legend->AddEntry(hOuTot[d],"Complete model","L");
    legend->AddEntry(g_data[d],"Data","LP"); 
    legend->Draw("same");
    
    TLatex* latex = new TLatex(0.40,0.60,Form("V_{OV} = %.1f V", OV));
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
    
    TLatex* latex2 = new TLatex(0.40,0.50,Form("C_{d} = %.1f fF, C_{q} = %.1f fF, R_{q} = %.0f k#Omega", sipmPars.Cd*1e15,sipmPars.Cq*1e15, sipmPars.Rq*1e-3));
    latex2 -> SetNDC();
    latex2 -> SetTextFont(42);
    latex2 -> SetTextSize(0.04);
    latex2 -> SetTextColor(kBlack);
    latex2 -> Draw("same");
    
    TLatex* latex3 = new TLatex(0.40,0.40,Form("BW = %.0f MHz", BW));
    latex3 -> SetNDC();
    latex3 -> SetTextFont(42);
    latex3 -> SetTextSize(0.04);
    latex3 -> SetTextColor(kBlack);
    latex3 -> Draw("same");
    
    c1->cd(2);
    gPad->SetLeftMargin(0.2); gPad->SetRightMargin(0.1); 
    gPad->SetTicks();
    gPad->SetLogy(); 
    hFrame ->Draw(); 
    hFrame -> SetMinimum(0.003);
    g_data[d] -> Draw("P,same");
    hOuTot[d] ->Draw("L,same");
    gPad->Update();
    c1->Update();
    
    c1->Print(Form("c_%s.png",sipmPars.label.c_str()));
    
    //outfile[d] = new TFile(Form("out.root"), "RECREATE");
    //outfile[d] -> cd();
    
    //outfile[d] -> WriteObject(g_data[d], Form("data_RL_%.0f_RF_%.0f_SiPM_%s_approximated_avg", RL[d], RF[d]*1e-3,tySiPM[d].c_str()));
    //outfile[d] -> WriteObject(fTot[d], Form("fft_function_RL_%.0f_RF_%.0f_SiPM_%s_approximated_avg", RL[d], RF[d]*1e-3,tySiPM[d].c_str()));
    
  }
  
  return;
}



// ************************
// FCN 
// ************************
void fcn(int& npar, double* gin, double& f, double* par, int iflag)
{
  // calculate chisquare
  double chisq = 0.;
  double delta = 0.;
  for(unsigned int iRun = 0; iRun < nRuns; ++iRun)
  {
    myInSignal(iRun, par);
    
    double x1 = hOuTot[iRun]->GetXaxis()->GetXmin();
    double x2 = hOuTot[iRun]->GetXaxis()->GetXmax();
    x1 = 20.; 
    x2 = 180.;
    
    double xx, yy, ey;
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
    
    // std::cout << "chisq: " << chisq;
    // for(int ipar = 0; ipar < npar; ++ipar)
    //   std::cout << "   ipar: " << ipar << " val: " << par[ipar];
    // std::cout << std::endl;
  }
  
  f = chisq;
}



void myFitPulse(const bool& minimize,
                const bool& fixAmp, const bool& fixT0, const bool& fixRq,
                const bool& fixBW)
{
  // Define fit model
  int nPars = 3*nRuns+1;
  gMinuit = new TMinuit(nPars);
  gMinuit->SetFCN(fcn);

  double arglist[10]; int ierflg = 0;

  
  // ERROR CHISQUARE
  arglist[0] = 1;
  gMinuit -> mnexcm("SET ERR", arglist, 1, ierflg);
  
  
  // INITIALIZE PARAMETERS
  double vstart[nPars]; 
  double vstep[nPars]; 
  for(int iRun = 0; iRun < nRuns; ++iRun)
  {
    SiPMParams sipmPars = SiPMParamsVec.at(iRun);
    
    // amp
    vstart[0+3*iRun] = 2500.;
    vstep[0+3*iRun] = 100.;
    gMinuit -> mnparm(0+3*iRun, Form("Amp"), vstart[0+3*iRun], vstep[0+3*iRun], 100.,100000., ierflg);
    if( fixAmp ) gMinuit->FixParameter(0+3*iRun);
    
    // t0
    vstart[1+3*iRun] = 20.;
    vstep[1+3*iRun] = 1.;
    gMinuit->mnparm(1+3*iRun, Form("t0"), vstart[1+3*iRun], vstep[1+3*iRun], 0. ,100.,ierflg);
    if( fixT0 ) gMinuit->FixParameter(1+3*iRun);
    
    // Rq
    vstart[2+3*iRun] = sipmPars.Rq;
    vstep[2+3*iRun] = 10e3;
    gMinuit->mnparm(2+3*iRun, Form("Rq"), vstart[2+3*iRun], vstep[2+3*iRun], 0. ,1000e3,ierflg);
    if( fixRq ) gMinuit->FixParameter(2+3*iRun);
  }
  
  vstart[3*nRuns] = 50.;
  vstep[3*nRuns] = 10.0;
  gMinuit->mnparm(3*nRuns, "BW", vstart[3*nRuns], vstep[3*nRuns], 0., 500., ierflg);
  if( fixBW ) gMinuit->FixParameter(3*nRuns);
  
  
  // initialize 
  arglist[0] = 1; 
  std::cout << "CALL FCN" << std::endl;
  gMinuit->mnexcm("CALL FCN",arglist ,1,ierflg);

  
  // Now ready for minimization step
  if( minimize )
  {
    arglist[0] = 500; arglist[1] = 1.; 
    gMinuit->mnexcm("MIN", arglist ,2,ierflg);   
  }
  
  
  // // Read parameters
  // pars_val = new double(nPars);
  // pars_err = new double(nPars);
  // pars_low = new double(nPars);
  // pars_hig = new double(nPars);
  // TString pname;
  // int iuint; 
  // for(int iPar = 0; iPar < 3; ++iPar)
  //   gMinuit->mnpout(iPar, pname, pars_val[iPar], pars_err[iPar], pars_low[iPar], pars_hig[iPar], iuint);
  
  
  return;
}






int main(int argc, char** argv)
{
  //----------------------
  // parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  
  int doMinimization = opts.GetOpt<int>("Input.doMinimization");
  int fixAmp = opts.GetOpt<int>("Input.fixAmp");
  int fixT0 = opts.GetOpt<int>("Input.fixT0");
  int fixRq = opts.GetOpt<int>("Input.fixRq");
  int fixBW = opts.GetOpt<int>("Input.fixBW");
  
  std::string dataFolder = opts.GetOpt<std::string>("Input.dataFolder");
  
  // get the SiPM and run list
  GetSiPMParsFromCfg(argv[1],SiPMParamsVec,runsVec);

  
  //------------------
  // define histograms
  nRuns = runsVec.size();
  for(unsigned int iRun = 0; iRun < nRuns; ++iRun)
  {
    hSiPM[iRun]  = new TH1D(Form("hSiPM_%d",iRun), "",npt,0,tmax); // SiPM pulse
    hBand[iRun]  = new TH1D(Form("hBand_%d",iRun), "",npt,0,tmax); // readout board bandwidth
    hOuTot[iRun] = new TH1D(Form("hOuTot_%d",iRun),"",npt,0,tmax); // output
  }
  

  //--------------
  // run the analsys
  myLoadData(runsVec,dataFolder); // fill graphs from files
  
  myFitPulse(bool(doMinimization),
             bool(fixAmp), bool(fixT0), bool(fixRq),
             bool(fixBW));
  
  myDrawFit(); // draw the fit result
  
  return 0;
}
