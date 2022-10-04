#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include <iostream>
#include <cmath>
#include <string>

#include "TString.h"

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
  std::string title;
  std::string label;
};

int CountUnique(const std::vector<int>& v);

void GetSiPMParsFromCfg(char* cfg, std::vector<SiPMParams>& vec, std::vector<int>& runs);
void GetFitParsFromCfg(char* cfg, const int& nRuns, std::map<int,int>& parIndex,
                       int& nPars_amp, int& nPars_t0, int& nPars_Rq, int& nPars_BW);

double funcLP(const double& xx, const double& tau);

double SiPMPulseShape(const double& x, const SiPMParams& sipmPars, const double& amp, const double& x0);

#endif
