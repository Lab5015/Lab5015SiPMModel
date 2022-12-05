#include "interface/Convolution.h"



// **********************************************************
// Pulse convolution: hFFT = hFun1 * hFun2
// - uses histograms with the same number of bins (samples)
// - and the same bin-width (reciprocal of the sampling freq)
// - iflag = 0 - convolution 
// - iflag = 1 - deconvolution
// 
// Example: 
//  TH1D *hPulse = new TH1D("hPulse","hPulse",npt,0,tmax);
//  hConvol(hFun1,hFun2,hPulse);       
//  TH1D * hDeOut =(TH1D*)hPulse->Clone("Detector output");
//
// Or (overwriting hFun1): 
//   hConvol(hFun1,hFun2,hFun1);
//        
// **********************************************************

void hConvol(TH1D* h1, TH1D* h2, TH1D* hFFT)
{
  // number of samples and sampling frequency
  int n = h1 -> GetNbinsX(); 
  double freq = 1./h1->GetBinWidth(1);
  
  TVirtualFFT::SetTransform(0);

  // Do the DFT of h1 
  TH1D* h1_freq = new TH1D(Form("%s_DFT_MAG",h1->GetName()),"",n,0.,freq);
  h1 -> FFT(h1_freq,"MAG R2C M"); 
  TVirtualFFT* FFT1 = TVirtualFFT::GetCurrentTransform();
  // Get complex coefficients 
  double* re1 = new double[n];
  double* im1 = new double[n];
  FFT1 -> GetPointsComplex(re1,im1);

  //  Do the DFT of h2 
  TH1D* h2_freq = new TH1D(Form("%s_DFT_MAG",h2->GetName()),"",n,0.,freq);
  h2 -> FFT(h2_freq,"MAG R2C M");
  TVirtualFFT* FFT2 = TVirtualFFT::GetCurrentTransform();
  // Get complex coefficients 
  double* re2 = new double[n];
  double* im2 = new double[n];
  FFT2 -> GetPointsComplex(re2,im2);
  
  // COMPLEX PRODUCT - (de)convolution in the conjugate domain - (including the DFT normalization term) 
  double* re_back = new double[n];
  double* im_back = new double[n];
  for(int i = 0; i < n; ++i)
  {
    re_back[i] = (re1[i]*re2[i]-im1[i]*im2[i]) / n ;
    im_back[i] = (im1[i]*re2[i]+re1[i]*im2[i]) / n ;
  }
  
  //Now let's make a backward transform:
  TVirtualFFT* FFT_back = TVirtualFFT::FFT(1, &n, "C2R M K");
  FFT_back -> SetPointsComplex(re_back,im_back);
  FFT_back -> Transform();
  TH1::TransformHisto(FFT_back,hFFT,"RE");
  
  delete h1_freq;
  delete h2_freq;
  
  return;
}



void hConvol(TH1D* h1, TF1* fBand_re, TF1* fBand_im, TH1D* hFFT)
{
  // number of samples and sampling frequency
  int n = h1 -> GetNbinsX(); 
  double freq = 1./h1->GetBinWidth(1);

  TVirtualFFT::SetTransform(0);

  // Do the DFT of h1 
  TH1D* h1_freq = new TH1D(Form("%s_DFT_MAG",h1->GetName()),"",n,0.,freq);
  h1 -> FFT(h1_freq,"MAG R2C M"); 
  TVirtualFFT* FFT1 = TVirtualFFT::GetCurrentTransform();
  // Get complex coefficients 
  double* re1 = new double[n];
  double* im1 = new double[n];
  FFT1 -> GetPointsComplex(re1,im1);
  
  double* re2 = new double[n];
  double* im2 = new double[n];
  
  // COMPLEX PRODUCT - (de)convolution in the conjugate domain - (including the DFT normalization term) 
  double* re_back = new double[n];
  double* im_back = new double[n];
  for(int i = 0; i < n; ++i)
  {
    double currFreq = freq / n * i;
    if( currFreq > freq/2. )
      currFreq = freq - currFreq;
    
    re2[i] = fBand_re -> Eval(currFreq);
    im2[i] = fBand_im -> Eval(currFreq);
    re_back[i] = (re1[i]*re2[i]-im1[i]*im2[i]) / n ;
    im_back[i] = (im1[i]*re2[i]+re1[i]*im2[i]) / n ;
  }
  
  //Now let's make a backward transform:
  TVirtualFFT* FFT_back = TVirtualFFT::FFT(1, &n, "C2R M K");
  FFT_back -> SetPointsComplex(re_back,im_back);
  FFT_back -> Transform();
  TH1::TransformHisto(FFT_back,hFFT,"RE");
  
  delete h1_freq;
  
  return;
}



void hConvol(TH1D* h1, TH1D* hBand_re, TH1D* hBand_im, TH1D* hFFT)
{
  // number of samples and sampling frequency
  int n = h1 -> GetNbinsX(); 
  double freq = 1./h1->GetBinWidth(1);

  TVirtualFFT::SetTransform(0);

  // Do the DFT of h1 
  TH1D* h1_freq = new TH1D(Form("%s_DFT_MAG",h1->GetName()),"",n,0.,freq);
  h1 -> FFT(h1_freq,"MAG R2C M"); 
  TVirtualFFT* FFT1 = TVirtualFFT::GetCurrentTransform();
  // Get complex coefficients 
  double* re1 = new double[n];
  double* im1 = new double[n];
  FFT1 -> GetPointsComplex(re1,im1);
  
  double* re2 = new double[n];
  double* im2 = new double[n];
  
  // COMPLEX PRODUCT - (de)convolution in the conjugate domain - (including the DFT normalization term) 
  double* re_back = new double[n];
  double* im_back = new double[n];
  for(int i = 0; i < n; ++i)
  {
    double currFreq = freq / n * i;
    
    re2[i] = hBand_re -> GetBinContent(hBand_re->FindBin(currFreq));
    im2[i] = hBand_im -> GetBinContent(hBand_im->FindBin(currFreq));
    re_back[i] = (re1[i]*re2[i]-im1[i]*im2[i]) / n ;
    im_back[i] = (im1[i]*re2[i]+re1[i]*im2[i]) / n ;
  }
  
  //Now let's make a backward transform:
  TVirtualFFT* FFT_back = TVirtualFFT::FFT(1, &n, "C2R M K");
  FFT_back -> SetPointsComplex(re_back,im_back);
  FFT_back -> Transform();
  TH1::TransformHisto(FFT_back,hFFT,"RE");
  
  delete h1_freq;
  
  return;
}



void hConvol(TH1D* h1, TF1* fBand_mag, TH1F* hBand_ph, TH1D* hFFT)
{
  // number of samples and sampling frequency
  int n = h1 -> GetNbinsX(); 
  double freq = 1./h1->GetBinWidth(1);

  TVirtualFFT::SetTransform(0);

  // Do the DFT of h1 
  TH1D* h1_freq = new TH1D(Form("%s_DFT_MAG",h1->GetName()),"",n,0.,freq);
  h1 -> FFT(h1_freq,"MAG R2C M"); 
  TVirtualFFT* FFT1 = TVirtualFFT::GetCurrentTransform();
  // Get complex coefficients 
  double* re1 = new double[n];
  double* im1 = new double[n];
  FFT1 -> GetPointsComplex(re1,im1);
  
  double* mag2 = new double[n];
  double* ph2 = new double[n];
  double* re2 = new double[n];
  double* im2 = new double[n];
  
  // COMPLEX PRODUCT - (de)convolution in the conjugate domain - (including the DFT normalization term) 
  double* re_back = new double[n];
  double* im_back = new double[n];
  for(int i = 0; i < n; ++i)
  {
    double currFreq = freq / n * i;
    
    mag2[i] = fBand_mag -> Eval(currFreq);
    ph2[i] = hBand_ph -> GetBinContent(hBand_ph->FindBin(currFreq));
    re2[i] = mag2[i]*cos(ph2[i]);
    im2[i] = mag2[i]*sin(ph2[i]);
    re_back[i] = (re1[i]*re2[i]-im1[i]*im2[i]) / n ;
    im_back[i] = (im1[i]*re2[i]+re1[i]*im2[i]) / n ;
  }
  
  //Now let's make a backward transform:
  TVirtualFFT* FFT_back = TVirtualFFT::FFT(1, &n, "C2R M K");
  FFT_back -> SetPointsComplex(re_back,im_back);
  FFT_back -> Transform();
  TH1::TransformHisto(FFT_back,hFFT,"RE");
  
  delete h1_freq;
  
  return;
}
