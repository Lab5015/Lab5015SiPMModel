#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include <iostream>
#include <cmath>
#include <string>

#include "TString.h"
#include "TGraphAsymmErrors.h"
#include "TH1D.h"

// SiPM parameters
struct SiPMParams {
  float OV;
  float Npe;
  float Rq;
  float Cq;
  float Rd;
  float Cd;
  float Nc;
  float Cg;
  float RL;
  float RF;
  float RqErr;
  float CqErr;
  float RdErr;
  float CdErr;
  float CgErr;
  float RLErr;
  float RFErr;
  std::string title;
  std::string label;
};

int CountUnique(const std::vector<int>& v);

std::map<int,int> RemapPars(const std::vector<int>& params);

void GetSiPMParsFromCfg(char* cfg, std::vector<SiPMParams>& vec, std::vector<std::string>* runs = NULL);
void GetFitParsFromCfg(char* cfg, const int& nRuns, std::map<int,int>& parIndex,
                       int& nPars_amp, int& nPars_t0, int& nPars_Rq, int& nPars_Cd, int& nPars_BW, int& nPars_L);

double funcLP(const double& xx, const double& tau);
double funcScint(const double& xx, const double& tau_r, const double& tau_d);

double SiPMPulseShape(const double& x, const SiPMParams& sipmPars, const double& amp, const double& x0, const int& whichComponent = -1);

#endif
