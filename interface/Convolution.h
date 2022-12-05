#ifndef CONVOLUTION_H
#define CONVOLUTION_H

#include "interface/Functions.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TF1.h"
#include "TVirtualFFT.h"



void hConvol(TH1D* hFun1, TH1D* hFun2, TH1D* hFFT);
void hConvol(TH1D* hFun1, TF1* fBand_re, TF1* fBand_im, TH1D* hFFT);
void hConvol(TH1D* hFun1, TH1D* hBand_re, TH1D* hBand_im, TH1D* hFFT);
void hConvol(TH1D* h1, TF1* fBand_mag, TH1F* hBand_ph, TH1D* hFFT);

#endif
