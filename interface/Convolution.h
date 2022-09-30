#ifndef CONVOLUTION_H
#define CONVOLUTION_H

#include "interface/Functions.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TF1.h"
#include "TVirtualFFT.h"



void hConvol(TH1D* hFun1, TH1D* hFun2, TH1D* hFFT, int iflag=0);
void hConvolBode(TH1D* hFun1, TProfile* hFun2_RE, TProfile* hFun2_IM, TH1D* hFFT, int iflag=0);
void hConvolBode(TH1D* hFun1, TF1* fBandwidth, TH1D* hFFT, int iflag=0);

#endif
