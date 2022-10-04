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




void hConvolBode(TH1D *hFun1, TProfile *hFun2_RE, TProfile *hFun2_IM, TH1D *hFFT, int iflag)
{
  // number of samples and sampling frequency
  int n = hFun1->GetNbinsX(); 
  double freq = 1./hFun1->GetBinWidth(1);
  double freq_RE = hFun2_RE->GetBinWidth(1);
  int n_RE = hFun2_RE->GetNbinsX();
  double freq_IM = hFun2_IM->GetBinWidth(1);
  int n_IM = hFun2_IM->GetNbinsX();
    
  // auxiliary histos and profiles
  TH1 *hP = 0;

  // Do the DFT of hFun1 
  TVirtualFFT::SetTransform(0);
  hP = hFun1->FFT(hP,"MAG1");
  TVirtualFFT *ffP = TVirtualFFT::GetCurrentTransform();
  // Get complex coefficients 
  double *reFun1 = new double[n];
  double *imFun1 = new double[n];
  ffP->GetPointsComplex(reFun1,imFun1);
  
  //get the real and imaginary part
  double *reFun2 = new double[n_RE];
  double *imFun2 = new double[n_IM];
  
  for(int ff = 0; ff < n_RE; ff++)
  {
     	reFun2[ff] = hFun2_RE -> GetBinContent(1+ff);
     	imFun2[ff] = hFun2_IM -> GetBinContent(1+ff);
  	//reFun2[ff] = 1;
  	//imFun2[ff] = 0;	
  }
  
  // COMPLEX PRODUCT - (de)convolution in the conjugate domain - (including the DFT normalization term) 
  double *re_back = new double[n];
  double *im_back = new double[n];
  for (int i=0;i<n;i++){
      double ifreq = i*freq/n;
      int iRE = int(ifreq/(freq_RE));
      int iIM = int(ifreq/(freq_IM));
      //std::cout  <<"i = " << i << ", " << "ifreq = " << ifreq << ", " << "iRE = " << iRE << std::endl;
      if(iRE >= n_RE) iRE = n_RE -1;
      if(iIM >= n_IM) iIM = n_IM -1;
      if (iflag==0) {
        re_back[i] = (reFun1[i]*reFun2[iRE]-imFun1[i]*imFun2[iIM]) / n ;
        im_back[i] = (imFun1[i]*reFun2[iRE]+reFun1[i]*imFun2[iIM]) / n ;
      }
      else {
        double norm = reFun2[i]*reFun2[iRE]+imFun2[i]*imFun2[iIM]; 
        re_back[i] = (reFun1[i]*reFun2[iRE]+imFun1[i]*imFun2[iIM]) / norm / n ;
        im_back[i] = (imFun1[i]*reFun2[iRE]-reFun1[i]*imFun2[iIM]) / norm / n ;
      }
    }
  
  //Now let's make a backward transform
  TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &n, "C2R M K");
  fft_back->SetPointsComplex(re_back,im_back);
  fft_back->Transform();
  TH1 *hb = 0;
  // gets the output and normalized to the bin-width (discrete transform norm)
  hb = TH1::TransformHisto(fft_back,hb,"Re");
  for (int ipt = 0; ipt<=n; ipt++){
    if (iflag==0) 
      hFFT -> SetBinContent(ipt+1, hb->GetBinContent(ipt+1) / freq );
    else
      hFFT -> SetBinContent(ipt+1, hb->GetBinContent(ipt+1) * freq );

  }

  delete hP;
  delete hb;

  return;
}



void hConvolBode(TH1D* hFun1, TF1* fBandwidth, TH1D* hFFT, int iflag)
{
  // number of samples and sampling frequency
  int n = hFun1->GetNbinsX(); 
  double freq = 1./hFun1->GetBinWidth(1);

  // auxiliary histos
  TH1 *hP = 0; 
  
  // Do the DFT of hFun1 
  TVirtualFFT::SetTransform(0);
  hP = hFun1->FFT(hP,"MAG1"); 
  TVirtualFFT *ffP = TVirtualFFT::GetCurrentTransform();
  // Get complex coefficients 
  double *reFun1 = new double[n];
  double *imFun1 = new double[n];
  ffP->GetPointsComplex(reFun1,imFun1);
  
  // COMPLEX PRODUCT - (de)convolution in the conjugate domain - (including the DFT normalization term) 
  double *re_back = new double[n];
  double *im_back = new double[n];
  for (int i=0;i<n;i++){

    double currFreq = freq / n * i;
    if( currFreq > freq/2. )
      currFreq = freq - currFreq;
    std::cout << "i: " << i << " currFreq = " << currFreq << "   func = " << fBandwidth->Eval(currFreq) << std::endl;
    
    // re_back[i] = reFun1[i] / n;
    // im_back[i] = imFun1[i] / n ;
    re_back[i] = reFun1[i]*fBandwidth->Eval(currFreq) / fBandwidth->Eval(0.) / n;
    im_back[i] = imFun1[i]*fBandwidth->Eval(currFreq) / fBandwidth->Eval(0.) / n;
  }
  
  //Now let's make a backward transform:
  TVirtualFFT *fft_back = TVirtualFFT::FFT(1, &n, "C2R M K");
  fft_back->SetPointsComplex(re_back,im_back);
  fft_back->Transform();
  TH1 *hb = 0;
  // gets the output and normalized to the bin-width (discrete transform norm)
  hb = TH1::TransformHisto(fft_back,hb,"Re");
  for (int ipt = 0; ipt<=n; ipt++){
    if (iflag==0) 
      hFFT -> SetBinContent(ipt+1, hb->GetBinContent(ipt+1) / freq );
    else
      hFFT -> SetBinContent(ipt+1, hb->GetBinContent(ipt+1) * freq );

  }
  delete hb;
  delete hP;
  
  return;
}
