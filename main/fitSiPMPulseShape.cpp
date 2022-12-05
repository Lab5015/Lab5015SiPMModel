#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include "interface/Convolution.h"
#include "interface/Functions.h"
#include "interface/FitUtils.h"

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
#include "TRatioPlot.h"
#include "TRandom3.h"

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
std::vector<std::string> runsVec;

bool isTOFHIR = false;

float amp0 = -1.;
float t00 = -1.;
float BW0 = -1.;
float BW20 = -1.;


// Define time axis and binning
const int npt     = 1000;
const double tmax = 200.;
const double freq = npt / tmax;

float tMin = 0.;
float tMax = 200.;

// TGraphs and TProfile for data
unsigned int nRuns = -1;
std::vector<TGraphErrors*> g_data;
std::vector<TH1D*> hSiPM;
std::vector<TH1D*> hLaser;
std::vector<TH1D*> hInTot;
std::vector<TH1D*> hBand;
std::vector<TH1D*> hBand2;
std::vector<TH1D*> hLowP;
std::vector<TH1D*> hDisc;
std::vector<TH1D*> hHigP;
std::vector<TH1D*> hDelP;
std::vector<TH1D*> hOuTot;



void myLoadData(std::vector<std::string> runs, const std::string& inFolder, const std::string& baseName, const std::string& psName)
{
  for(unsigned int iRun = 0; iRun < nRuns; ++iRun)
  {
    TFile* inFile = TFile::Open(Form("%s/%s_run%s.root",inFolder.c_str(),baseName.c_str(),runsVec[iRun].c_str()),"READ");
    
    if( !isTOFHIR )
      {
	TProfile* p_data = (TProfile*)inFile->Get(Form("%s",psName.c_str()));
	
	double x,y, ey;
	for( int iBin = 1; iBin <= p_data->GetNbinsX(); ++iBin)
	  {
	    x = p_data -> GetBinCenter(iBin); 
	    y = p_data -> GetBinContent(iBin);
	    ey = p_data -> GetBinError(iBin);
	    g_data[iRun] -> SetPoint(g_data[iRun]->GetN(),x,y);
	    g_data[iRun] -> SetPointError(g_data[iRun]->GetN()-1,0.,ey);
	  }
      }
    else
      {
	g_data[iRun] = (TGraphErrors*)( inFile->Get(Form("%s",psName.c_str())) );
	double x,y, ey;
	for(int iPoint = 0; iPoint < g_data[iRun]->GetN(); ++iPoint)
	  {
	    g_data[iRun] -> GetPoint(iPoint,x,y);
	    g_data[iRun] -> SetPoint(iPoint,x+20.,y*1e-6);
	    g_data[iRun] -> SetPointError(iPoint,0.,0.2*1e-6);
	  }
      }
    inFile->Close();
  }
  
  return;
}



// ******************************
// SiPM model on an external load
// ******************************
void myInSignal(const int& iRun, double* par = NULL)
{
  SiPMParams* sipmPars = &(SiPMParamsVec.at(iRun));

  double amp = amp0;
  double t0  = t00;
  double Rq  = (*sipmPars).Rq;
  double Cd  = (*sipmPars).Cd;
  double BW  = BW0;
  double BW2 = BW20;
  
  if( par )
    {
      if( parIndex[0+6*iRun] >= 0 ) amp = par[parIndex[0+6*iRun]];
      if( parIndex[1+6*iRun] >= 0 )  t0 = par[parIndex[1+6*iRun]];
      if( parIndex[2+6*iRun] >= 0 )  Rq = par[parIndex[2+6*iRun]];
      if( parIndex[3+6*iRun] >= 0 )  Cd = par[parIndex[3+6*iRun]];
      if( parIndex[4+6*iRun] >= 0 )  BW = par[parIndex[4+6*iRun]];
      if( parIndex[5+6*iRun] >= 0 ) BW2 = par[parIndex[5+6*iRun]];
      
      amp0 = amp;
      t00 = t0;
      BW0 = BW;
      BW20 = BW2;
    }
  
  (*sipmPars).Rq = Rq;
  (*sipmPars).Cd = Cd;
  
  fitPars[0+6*iRun] = amp;
  fitPars[1+6*iRun] = t0;
  fitPars[2+6*iRun] = Rq;
  fitPars[3+6*iRun] = Cd;
  fitPars[4+6*iRun] = BW;
  fitPars[5+6*iRun] = BW2;
  
  // bandwidth low-pass filter
  double tau_band = 1./(BW*3.14159*2.)*1e3; // nanosenconds (BW is in MHz)
  
  // L/R low-pass filter
  double tau_band2 = 1./((BW+BW2)*3.14159*2.)*1e3; // nanosenconds (BW is in MHz)
  
  double tau1 = 1./(BW*3.14159*2.)*1e3; // nanosenconds (BW is in MHz)
  double tau2 = 0.133;
  double tau3 = 1./(BW2*3.14159*2.)*1e3; // nanosenconds (BW is in MHz)
  //double tauE = 0.5736; // from Elmore
  double tauE = 0.8736; // from Elmore
  
  // Hitograms for FFT analysis 
  for(int ipt = 0; ipt <= npt; ++ipt)
  { 
    double tt = (double)ipt * (tmax/npt);
    
    hSiPM[iRun] -> SetBinContent(ipt+1, 0.5*((*sipmPars).RF)*SiPMPulseShape(tt,*sipmPars,amp,t0)); // 0.5*RF is the total transimpedance
    
    // AdvanSid bandwidth
    hBand[iRun]     -> SetBinContent(ipt+1, (1.-0.5/(freq*tau_band))*funcLP(tt,tau_band));
    hBand2[iRun]    -> SetBinContent(ipt+1, (1.-0.5/(freq*tau_band2))*funcLP(tt,tau_band2));
    
    // TOFHIR bandwidth
    hLowP[iRun] -> SetBinContent(ipt+1, (1.-0.5/(freq*tau1))*funcLP(tt,tau1) );
    hDisc[iRun] -> SetBinContent(ipt+1, (1.-0.5/(freq*tau2))*funcLP(tt,tau2) );
    if( ipt == 0 )
      {
	hHigP[iRun] -> SetBinContent(ipt+1, freq-(1.-0.5/(freq*tau3))*funcLP(tt,tau3) );
	hDelP[iRun] -> SetBinContent(ipt+1, freq-(1.-0.5/(freq*tauE))*funcLP(tt,tauE) );
      }
    else
      {
	hHigP[iRun] -> SetBinContent(ipt+1, 0.-(1.-0.5/(freq*tau3))*funcLP(tt,tau3) );
	hDelP[iRun] -> SetBinContent(ipt+1, 0.-(1.-0.5/(freq*tauE))*funcLP(tt,tauE) );
      }
  }
  
  hBand[iRun] -> Scale(1./hBand[iRun]->Integral());
  hBand2[iRun] -> Scale(1./hBand2[iRun]->Integral());
  
  // hLowP[iRun] -> Scale(1./hLowP[iRun]->Integral());
  // hDisc[iRun] -> Scale(1./hDisc[iRun]->Integral());
  // hHigP[iRun] -> Scale(1./hHigP[iRun]->Integral());
  // hDelP[iRun] -> Scale(1./hDelP[iRun]->Integral());
  
  
  // convolutions
  if( !isTOFHIR )
    {
      hConvol(hSiPM[iRun],hLaser[iRun],hInTot[iRun]);
      hConvol(hInTot[iRun],hBand[iRun],hOuTot[iRun]);
      hConvol(hOuTot[iRun],hBand2[iRun],hOuTot[iRun]);
    }
  
  else
    {
      hConvol(hSiPM[iRun],hLaser[iRun],hInTot[iRun]);
      hConvol(hInTot[iRun],hHigP[iRun],hOuTot[iRun]);
      hConvol(hOuTot[iRun],hLowP[iRun],hOuTot[iRun]);
      hConvol(hOuTot[iRun],hDisc[iRun],hOuTot[iRun]);
      hConvol(hOuTot[iRun],hDelP[iRun],hOuTot[iRun]);
    }
  
  // TEMP!!!
  //hConvol(hSiPM[iRun],hLaser[iRun],hOuTot[iRun]);
  
  // TF1* f_band_re = new TF1("f_band_re","[0]*sin([1]*x+[2])*exp(-1.*[3]*x*x)",0.,2.5);
  // f_band_re -> SetParameters(2.12208e00,5.12866e01,1.49266e00,8.76131e00);
  // TF1* f_band_im = new TF1("f_band_im","[0]*sin([1]*x+[2])*exp(-1.*[3]*x*x)",0.,2.5);
  // f_band_im -> SetParameters(2.23291e00,5.03354e01,3.25655e00,9.07746e00);
  // TF1* f_band_mag = new TF1("f_band_mag","exp([0]+[1]*x+[2]*x*x+[3]*x*x*x)",0.,2.5);
  // f_band_mag -> SetParameters(8.79094e-01,-2.85437e00,6.23893e00,-1.32199e01);
  //hConvol(hInTot[iRun],f_band_re,f_band_im,hOuTot[iRun]);
  
  // //TFile* inFile_band = TFile::Open("/home/cmsdaq/DRS4/drs-5.0.6-lab5015/reco_pulses/pulses_run0149.root");
  // TFile* inFile_band = TFile::Open("data/effectiveBW_147.root");
  // // TH1F* hBand_time = (TH1F*) inFile_band -> Get("h_ps_entry0");
  // // hBand_time -> Scale(1./hBand_time->Integral());
  // //TH1F* hBand_ph = (TH1F*) inFile_band -> Get("h_DFT_PH_avg");
  // //hConvol(hInTot[iRun],f_band_mag,hBand_ph,hOuTot[iRun]);
  // hConvol(hInTot[iRun],hBand_re[iRun],hBand_im[iRun],hOuTot[iRun]);
  // //hConvol(hInTot[iRun],(TH1D*)hBand_time,hOuTot[iRun]);
  
  // inFile_band -> Close();
  
  return;
} 



// ******************************
// Compute error band
// ******************************

TGraphAsymmErrors* GetConfidenceInterval(const int& iRun, const int& nToys)
{
  TGraphAsymmErrors* ret = new TGraphAsymmErrors();
  
  int N = hOuTot[iRun] -> GetNbinsX();
  for(int iBin = 1; iBin <= N; ++iBin)
    {
      double center = hOuTot[iRun] -> GetBinCenter(iBin);
      double content = hOuTot[iRun] -> GetBinContent(iBin);
      ret -> SetPoint(iBin-1,center,content);
    }
  
  SiPMParams* sipmPars = &(SiPMParamsVec.at(iRun));
  TRandom3 r;
  
  std::map<int,std::vector<double> > vals;
  for(int iToy = 0; iToy < nToys; ++iToy)
    {
      float RqStart = (*sipmPars).Rq;
      float RqErr = (*sipmPars).RqErr;
      (*sipmPars).Rq = RqStart + 1.*r.Gaus(0.,RqErr);
      
      float CqStart = (*sipmPars).Cq;
      float CqErr = (*sipmPars).CqErr;
      (*sipmPars).Cq = CqStart + 1.*r.Gaus(0.,CqErr);
      
      float RdStart = (*sipmPars).Rd;
      float RdErr = (*sipmPars).RdErr;
      (*sipmPars).Rd = RdStart + 1.*r.Gaus(0.,RdErr);
      
      float CdStart = (*sipmPars).Cd;
      float CdErr = (*sipmPars).CdErr;
      (*sipmPars).Cd = CdStart + 1.*r.Gaus(0.,CdErr);
      
      float CgStart = (*sipmPars).Cg;
      float CgErr = (*sipmPars).CgErr;
      (*sipmPars).Cg = CgStart + 1.*r.Gaus(0.,CgErr);
      
      float RLStart = (*sipmPars).RL;
      float RLErr = (*sipmPars).RLErr;
      (*sipmPars).RL = RLStart + 1.*r.Gaus(0.,RLErr);
      
      myInSignal(iRun,NULL);
      
      for(int iBin = 1; iBin <= N; ++iBin)
	{
	  vals[iBin].push_back( hOuTot[iRun]->GetBinContent(iBin) );
	}
      
      (*sipmPars).Rq = RqStart;
      (*sipmPars).Cq = CqStart;
      (*sipmPars).Rd = RdStart;
      (*sipmPars).Cd = CdStart;
      (*sipmPars).Cg = CgStart;
      (*sipmPars).RL = RLStart;
    }

  for(int iBin = 1; iBin <= N; ++iBin)
    {
      std::vector<double> vec = vals[iBin];
      std::pair<double,double> cl = FindSmallestInterval(&vec);
      
      double x,y;
      ret -> GetPoint(iBin-1,x,y);
      ret -> SetPointEYhigh(iBin-1,cl.second-0.5*(cl.first+cl.second));
      ret -> SetPointEYlow(iBin-1,0.5*(cl.first+cl.second)-cl.first);
    }
  
  return ret;
}




// ******************************
// Draw result
// ******************************
void myDrawFit(const std::string& plotFolder, const std::string& commonLabel)
{
  system(Form("mkdir -p %s",plotFolder.c_str())); 
  
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
	if( Amax < hOuTot[iRun]->GetBinContent(iBin) ) 
	  Amax = hOuTot[iRun]->GetBinContent(iBin);
      } 
    
    // Draw
    TCanvas* c1 = new TCanvas(Form("c%d",iRun), Form("c%d",iRun), 1800., 750.); 
    c1->Divide(2,1);
    c1 -> cd(1);
    gPad->SetLeftMargin(0.2); gPad->SetRightMargin(0.1);
    gPad->SetTicks();
    
    TH1F* hFrame = (TH1F*)gPad->DrawFrame(tMin,-0.2*Amax,tMax,1.2*Amax);
    hFrame -> SetTitle(Form("%s",sipmPars.title.c_str())); 
    hFrame -> GetXaxis()->SetTitle("time [ns]");
    if( !isTOFHIR ) hFrame -> GetYaxis()->SetTitle("Amplitude [V]"); 
    else            hFrame -> GetYaxis()->SetTitle("Current [A]"); 
    hFrame -> GetXaxis()->SetLabelSize(0.045); 
    hFrame -> GetYaxis()->SetLabelSize(0.045); 
    hFrame -> GetXaxis()->SetTitleSize(0.050); 
    hFrame -> GetYaxis()->SetTitleSize(0.050);
    hFrame -> Draw();
    
    g_data[iRun] -> SetMarkerSize(0.5);
    g_data[iRun] -> SetMarkerStyle(8);
    g_data[iRun] -> Draw("P, same");
    
    hOuTot[iRun] ->SetLineColor(kGreen+1);
    hOuTot[iRun] ->SetLineWidth(2);
    hOuTot[iRun] ->DrawCopy("L3,same");
    
    TLegend* legend = new TLegend(0.55,0.79,0.85,0.89);
    legend->SetBorderSize(0); 
    legend->SetTextSize(0.035);
    legend->AddEntry(hOuTot[iRun],"complete model","L");
    legend->AddEntry(g_data[iRun],"data","LP"); 
    legend->Draw("same");
    
    TLatex* latex = new TLatex(0.40,0.70,Form("V_{OV} = %.1f V", sipmPars.OV));
    latex -> SetNDC();
    latex -> SetTextFont(42);
    latex -> SetTextSize(0.04);
    latex -> SetTextColor(kBlack);
    latex -> Draw("same");
    
    TLatex* latex1 = new TLatex(0.40,0.65,Form("N_{c} = %.0f   C_{g} = %.0f pF", sipmPars.Nc, sipmPars.Cg*1e12));
    latex1 -> SetNDC();
    latex1 -> SetTextFont(42);
    latex1 -> SetTextSize(0.04);
    latex1 -> SetTextColor(kBlack);
    latex1 -> Draw("same");
    
    TLatex* latex2 = new TLatex(0.40,0.60,Form("C_{d} = %.1f fF, C_{q} = %.1f fF, R_{q} = %.0f k#Omega", fitPars[3+6*iRun]*1e15,sipmPars.Cq*1e15, fitPars[2+6*iRun]*1e-3));
    latex2 -> SetNDC();
    latex2 -> SetTextFont(42);
    latex2 -> SetTextSize(0.04);
    latex2 -> SetTextColor(kBlack);
    latex2 -> Draw("same");
    
    TLatex* latex3;
    if( !isTOFHIR ) latex3 = new TLatex(0.40,0.50,Form("BW_{1} = %.0f MHz   BW_{2} = %.0f MHz", fitPars[4+6*iRun],fitPars[4+6*iRun]+fitPars[5+6*iRun]));
    else            latex3 = new TLatex(0.40,0.50,Form("BW_{1} = %.0f MHz   BW_{2} = %.0f MHz", fitPars[4+6*iRun],fitPars[5+6*iRun]));
    latex3 -> SetNDC();
    latex3 -> SetTextFont(42);
    latex3 -> SetTextSize(0.04);
    latex3 -> SetTextColor(kBlack);
    latex3 -> Draw("same");

    TLatex* latex4 = new TLatex(0.22,0.15,Form("#splitline{amp. = %.2f}{t_{0} = %.2f ns}",fitPars[0+6*iRun],5.+fitPars[1+6*iRun]));
    latex4 -> SetNDC();
    latex4 -> SetTextFont(42);
    latex4 -> SetTextSize(0.03);
    latex4 -> SetTextColor(kBlack);
    latex4 -> Draw("same");
    
    // TH1D* h_data = new TH1D(Form("h_data_%d",iRun), "",npt,0,tmax);
    // for(int iBin = 1; iBin <= hOuTot[iRun]->GetNbinsX(); ++iBin)
    //   {
    // 	hOuTot[iRun] -> SetBinError(iBin,0.01);
    // 	h_data -> SetBinContent(iBin,g_data[iRun]->Eval(hOuTot[iRun]->GetBinCenter(iBin)));
    // 	h_data -> SetBinError(iBin,0.001);
    //   }
    // TRatioPlot* rp = new TRatioPlot(h_data, hOuTot[iRun]);
    // rp -> GetLowerRefYaxis() -> SetRangeUser(-10.,10.);
    // rp -> Draw("same");
    
    c1 -> cd(2);
    TH1F* hFrame2 = (TH1F*)gPad->DrawFrame(tMin,Amax/100.,tMax,2.*Amax);
    hFrame2 -> SetTitle(Form("%s",sipmPars.title.c_str())); 
    hFrame2 -> GetXaxis()->SetTitle("time [ns]");
    if( !isTOFHIR ) hFrame2 -> GetYaxis()->SetTitle("Amplitude [V]"); 
    else            hFrame2 -> GetYaxis()->SetTitle("Current [A]"); 
    hFrame2 -> GetXaxis()->SetLabelSize(0.045); 
    hFrame2 -> GetYaxis()->SetLabelSize(0.045); 
    hFrame2 -> GetXaxis()->SetTitleSize(0.050); 
    hFrame2 -> GetYaxis()->SetTitleSize(0.050);
    hFrame2 -> Draw();
    
    gPad->SetLeftMargin(0.2); gPad->SetRightMargin(0.1); 
    gPad->SetTicks();
    gPad->SetLogy(); 
    
    g_data[iRun] -> Draw("P,same");
    hOuTot[iRun] -> DrawCopy("L,same");
    gPad->Update();
    
    TGraphAsymmErrors* errorBand = GetConfidenceInterval(iRun,100);
    errorBand -> SetLineWidth(0);
    errorBand -> SetFillColor(kGreen+1);
    errorBand -> SetFillStyle(3001);
    c1 -> cd(1);
    errorBand -> Draw("L3,same");
    c1 -> cd(2);
    errorBand -> Draw("L3,same");
    
    c1->Update();
    c1 -> Print(Form("%s/c_%s_Vov%.2f_%s.png",plotFolder.c_str(),sipmPars.label.c_str(),sipmPars.OV,commonLabel.c_str()));
    
    TFile outFile(Form("%s/plots_%s_Vov%.2f_%s.root",plotFolder.c_str(),sipmPars.label.c_str(),sipmPars.OV,commonLabel.c_str()), "RECREATE");
    outFile.cd();
    g_data[iRun] -> Write("g_data");
    hSiPM[iRun] -> Write();
    hLaser[iRun] -> Write();
    hInTot[iRun] -> Write();
    hBand[iRun] -> Write();
    hBand2[iRun] -> Write();
    hOuTot[iRun] -> Write();
    errorBand -> Write("g_errorBand");
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
    double ymax = 0.;
    for(int iPoint = 0; iPoint < g_data[iRun]->GetN(); ++iPoint)
    {
      g_data[iRun] -> GetPoint(iPoint,xx,yy);
      if( yy > ymax )
      {
        ymax = yy;
      }
    }
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
      if( yy > 0.02*ymax )
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
    // std::cout << ">>> chisq: " << chisq;
    // for(int iPar = 0; iPar < npar; ++iPar)
    //   std::cout << "   iPar: " << iPar << " val: " << par[iPar];
    // std::cout << std::endl;
  }
  
  f = chisq;
}



// ************************
// fit the pulse shape
// ************************
void myFitPulse(const bool& minimize,
                const int& nPars_amp, const int& nPars_t0, const int& nPars_Rq, const int& nPars_Cd, const int& nPars_BW, const int& nPars_BW2)
{
  int nPars = nPars_amp+nPars_t0+nPars_Rq+nPars_Cd+nPars_BW+nPars_BW2;
  
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
  for(int iPar = 0; iPar < nPars_amp; ++iPar) fitter -> SetParameter(iPar, Form("Amp"), amp0, 0.01, 1e-6,10.);
  for(int iPar = 0; iPar < nPars_t0; ++iPar)  fitter -> SetParameter(nPars_amp+iPar, Form("t0"), t00, 0.5, -50., 50.);
  for(int iPar = 0; iPar < nPars_Rq; ++iPar)  fitter -> SetParameter(nPars_amp+nPars_t0+iPar, Form("Rq"), 450e3, 10e3, 0. ,1000e3);
  for(int iPar = 0; iPar < nPars_Cd; ++iPar)  fitter -> SetParameter(nPars_amp+nPars_t0+nPars_Rq+iPar, Form("Cd"), 12.4e-15, 0.2e-15, 1e-15 ,30e-15);
  for(int iPar = 0; iPar < nPars_BW; ++iPar)  fitter -> SetParameter(nPars_amp+nPars_t0+nPars_Rq+nPars_Cd+iPar, "BW", BW0, 2., 0., 200.);
  for(int iPar = 0; iPar < nPars_BW2;++iPar)  fitter -> SetParameter(nPars_amp+nPars_t0+nPars_Rq+nPars_Cd+nPars_BW+iPar, "BW2", BW20, 5., 0., 500.);
  
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
  std::string baseName = opts.GetOpt<std::string>("Input.baseName");
  std::string psName = opts.GetOpt<std::string>("Input.psName");
  std::string commonLabel = opts.GetOpt<std::string>("Input.commonLabel");
  std::string plotFolder = opts.GetOpt<std::string>("Output.plotFolder");
  
  isTOFHIR = (bool)(opts.GetOpt<int>("Input.isTOFHIR"));
  amp0 = opts.GetOpt<float>("Input.amp0");
  t00 = opts.GetOpt<float>("Input.t00");
  BW0 = opts.GetOpt<float>("Input.BW0");
  BW20 = opts.GetOpt<float>("Input.BW20");
  
  tMin = opts.GetOpt<float>("Input.tMin");
  tMax = opts.GetOpt<float>("Input.tMax");
  
  
  //--------------------------  
  // get the SiPM and run list
  GetSiPMParsFromCfg(argv[1],SiPMParamsVec,&runsVec);
  nRuns = runsVec.size();
  fitPars.reserve(6*nRuns);
  

  //-----------------------  
  // get the fit parameters
  int nPars_amp = -1;
  int nPars_t0 = -1;
  int nPars_Rq = -1;
  int nPars_Cd = -1;
  int nPars_BW = -1;
  int nPars_BW2 = -1;
  GetFitParsFromCfg(argv[1], nRuns, parIndex, nPars_amp, nPars_t0, nPars_Rq, nPars_Cd, nPars_BW, nPars_BW2);
  
  
  //-------------------
  // define histograms
  for(unsigned int iRun = 0; iRun < nRuns; ++iRun)
    {
    g_data.push_back( new TGraphErrors() );
    hSiPM.push_back(  new TH1D(Form("hSiPM_%d",iRun), "",npt,0,tmax) ); // SiPM pulse
    hLaser.push_back( new TH1D(Form("hLaser_%d",iRun),"",npt,0,tmax) ); // laser response
    hInTot.push_back( new TH1D(Form("hInTot_%d",iRun),"",npt,0,tmax) ); // input pulse
    hBand.push_back(  new TH1D(Form("hBand_%d",iRun), "",npt,0,tmax) ); // readout board bandwidth
    hBand2.push_back( new TH1D(Form("hBand2_%d",iRun),"",npt,0,tmax) ); // readout board bandwidth
    hLowP.push_back(  new TH1D(Form("hLowP_%d",iRun), "",npt,0,tmax) ); // readout board bandwidth
    hDisc.push_back(  new TH1D(Form("hDisc_%d",iRun), "",npt,0,tmax) ); // readout board bandwidth
    hHigP.push_back(  new TH1D(Form("hHigP_%d",iRun), "",npt,0,tmax) ); // readout board bandwidth
    hDelP.push_back(  new TH1D(Form("hDelP_%d",iRun), "",npt,0,tmax) ); // readout board bandwidth
    hOuTot.push_back( new TH1D(Form("hOuTot_%d",iRun),"",npt,0,tmax) ); // output
  }
  
  
  //-----------------------
  // fill static histograms
  // laser response
  double sigma = 0.021; // 21 ps sigma for our laser
  double mu = 5.; // arbitrary
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
  myLoadData(runsVec,dataFolder,baseName,psName); // fill graphs from files
  
  myFitPulse(bool(doMinimization),
             nPars_amp, nPars_t0, nPars_Rq, nPars_Cd, nPars_BW, nPars_BW2);
  
  myDrawFit(plotFolder,commonLabel); // draw the fit result
  
  return 0;
}
