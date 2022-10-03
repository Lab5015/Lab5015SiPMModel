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

int nRuns = 20;

TVirtualFitter* fitter;
std::map<int,int> parIndex;
std::vector<double> fitPars;

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

  double amp = 25.;
  double t0 = 10.;
  double Rq = (*sipmPars).Rq;
  double BW = 50.;
    
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
    
    hSiPM[iRun] -> SetBinContent(ipt+1, (*sipmPars).RL*SiPMPulseShape(tt,*sipmPars,OV,amp,t0));
    hBand[iRun] -> SetBinContent(ipt+1, funcLP(tt,tau_band));
  }
  
  hConvol(hSiPM[iRun],hBand[iRun],hOuTot[iRun]);
  
  return;
} 



// Draw
void myDrawFit()
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
    TCanvas* c1 = new TCanvas(Form("c%d", iRun),Form("c%d", iRun),1800.,750.); 
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
    
    c1->cd(2);
    gPad->SetLeftMargin(0.2); gPad->SetRightMargin(0.1); 
    gPad->SetTicks();
    gPad->SetLogy(); 
    hFrame ->Draw(); 
    hFrame -> SetMinimum(0.002);
    hFrame -> GetXaxis() -> SetRangeUser(0.,180.);
    g_data[iRun] -> Draw("P,same");
    hOuTot[iRun] ->Draw("L,same");
    gPad->Update();
    c1->Update();
    
    c1->Print(Form("c_%s.png",sipmPars.label.c_str()));
    
    //outfile[iRun] = new TFile(Form("out.root"), "RECREATE");
    //outfile[iRun] -> cd();
    
    //outfile[iRun] -> WriteObject(g_data[iRun], Form("data_RL_%.0f_RF_%.0f_SiPM_%s_approximated_avg", RL[d], RF[d]*1e-3,tySiPM[d].c_str()));
    //outfile[iRun] -> WriteObject(fTot[d], Form("fft_function_RL_%.0f_RF_%.0f_SiPM_%s_approximated_avg", RL[d], RF[d]*1e-3,tySiPM[d].c_str()));
    
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
  arglist[0] = 2;
  fitter -> ExecuteCommand("SET PRINT",arglist,1);
  
  // ERROR CHISQUARE
  arglist[0] = 1;
  fitter -> ExecuteCommand("SET ERR", arglist,1);
  
  // INITIALIZE PARAMETERS
  for(int iPar = 0; iPar < nPars_amp; ++iPar)
    fitter -> SetParameter(iPar, Form("Amp"), 25., 1., 1.,100.);
  for(int iPar = 0; iPar < nPars_t0; ++iPar)
    fitter -> SetParameter(nPars_amp+iPar, Form("t0"), 20., 0.5, 10. ,30.);
  for(int iPar = 0; iPar < nPars_Rq; ++iPar)
    fitter -> SetParameter(nPars_amp+nPars_t0+iPar, Form("Rq"), 450e3, 10e3, 0. ,1000e3);
  for(int iPar = 0; iPar < nPars_BW; ++iPar)
    fitter -> SetParameter(nPars_amp+nPars_t0+nPars_Rq+iPar, "BW", 49., 2., 0., 200.);
  
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

  
  // get the SiPM and run list
  GetSiPMParsFromCfg(argv[1],SiPMParamsVec,runsVec);
  nRuns = runsVec.size();
  

  // get the fit parameters
  fitPars.reserve(4*nRuns);
  std::vector<int> params_amp;
  std::vector<int> params_t0;
  std::vector<int> params_Rq;
  std::vector<int> params_BW;
  std::vector<std::string> SiPMList = opts.GetOpt<std::vector<std::string> >("Input.SiPMList");
  for(auto iSiPM : SiPMList)
  {
    int param_amp = opts.GetOpt<int>(Form("%s.paramAmp",iSiPM.c_str())); params_amp.push_back( param_amp );
    int param_t0  = opts.GetOpt<int>(Form("%s.paramT0",iSiPM.c_str()));  params_t0.push_back( param_t0 );
    int param_Rq  = opts.GetOpt<int>(Form("%s.paramRq",iSiPM.c_str()));  params_Rq.push_back( param_Rq );
    int param_BW  = opts.GetOpt<int>(Form("%s.paramBW",iSiPM.c_str()));  params_BW.push_back( param_BW );
  }
  int nPars_amp = CountUnique(params_amp);
  int nPars_t0 = CountUnique(params_t0);
  int nPars_Rq = CountUnique(params_Rq);
  int nPars_BW = CountUnique(params_BW);
  std::cout << "nPars_amp: " << nPars_amp << std::endl;
  std::cout << "nPars_t0: " << nPars_t0 << std::endl;
  std::cout << "nPars_Rq: " << nPars_Rq << std::endl;
  std::cout << "nPars_BW: " << nPars_BW << std::endl;
  for(unsigned int iRun = 0; iRun < nRuns; ++iRun)
  {
    if( params_amp.at(iRun) >= 0 ) parIndex[0+4*iRun] = params_amp.at(iRun);                                  else parIndex[0+4*iRun] = -1;
    if( params_t0.at(iRun)  >= 0 ) parIndex[1+4*iRun] = nPars_amp + params_t0.at(iRun);                       else parIndex[1+4*iRun] = -1;
    if( params_Rq.at(iRun)  >= 0 ) parIndex[2+4*iRun] = nPars_amp + nPars_t0 + params_Rq.at(iRun);            else parIndex[2+4*iRun] = -1;
    if( params_BW.at(iRun)  >= 0 ) parIndex[3+4*iRun] = nPars_amp + nPars_t0 + nPars_Rq + params_BW.at(iRun); else parIndex[3+4*iRun] = -1;
    std::cout << "iRun: " << iRun << "  parIndex[0] = " << parIndex[0+4*iRun] << "   parIndex[1] = " << parIndex[1+4*iRun] << "   parIndex[2] = " << parIndex[2+4*iRun] << "   parIndex[3] = " << parIndex[3+4*iRun] << std::endl;
  }
  
  
  //------------------
  // define histograms
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
             nPars_amp, nPars_t0, nPars_Rq, nPars_BW);
  
  myDrawFit(); // draw the fit result
  
  return 0;
}
