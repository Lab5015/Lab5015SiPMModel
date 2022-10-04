#include "interface/Functions.h"



int CountUnique(const std::vector<int>& v)
{
  std::vector<int> v_sorted(v);
  std::sort(v_sorted.begin(),v_sorted.end());

  auto last = std::unique(v_sorted.begin(),v_sorted.end());
  v_sorted.erase(last,v_sorted.end());

  int count = 0;
  for(auto val : v_sorted)
    if( val >= 0 ) count += 1;

  return count;
}



void GetSiPMParsFromCfg(char* cfg, std::vector<SiPMParams>& vec, std::vector<int>& runs)
{  
  // parse the config file
  CfgManager opts;
  opts.ParseConfigFile(cfg);
  
  std::vector<std::string> SiPMList = opts.GetOpt<std::vector<std::string> >("Input.SiPMList");
  
  for(auto iSiPM : SiPMList)
  {
    SiPMParams pars;
    pars.Npe   = opts.GetOpt<float>(Form("%s.Npe",iSiPM.c_str()));
    pars.OV    = opts.GetOpt<float>(Form("%s.OV",iSiPM.c_str()));
    pars.Rq    = opts.GetOpt<float>(Form("%s.Rq",iSiPM.c_str()));
    pars.Cq    = opts.GetOpt<float>(Form("%s.Cq",iSiPM.c_str()));
    pars.Rd    = opts.GetOpt<float>(Form("%s.Rd",iSiPM.c_str()));
    pars.Cd    = opts.GetOpt<float>(Form("%s.Cd",iSiPM.c_str()));
    pars.Nc    = opts.GetOpt<float>(Form("%s.Nc",iSiPM.c_str()));
    pars.Cg    = opts.GetOpt<float>(Form("%s.Cg",iSiPM.c_str()));
    pars.RL    = opts.GetOpt<float>(Form("%s.RL",iSiPM.c_str()));
    std::vector<std::string> tokens = opts.GetOpt<std::vector<std::string> >(Form("%s.title",iSiPM.c_str()));
    std::string title = "";
    for(auto token : tokens) title += token + " ";
    pars.title = title;
    pars.label = opts.GetOpt<std::string>(Form("%s.label",iSiPM.c_str()));
    vec.push_back(pars);

    int run = opts.GetOpt<int>(Form("%s.run",iSiPM.c_str()));
    runs.push_back(run);
  }

  return;
}



void GetFitParsFromCfg(char* cfg, const int& nRuns, std::map<int,int>& parIndex,
                       int& nPars_amp, int& nPars_t0, int& nPars_Rq, int& nPars_BW)
{
  // parse the config file
  CfgManager opts;
  opts.ParseConfigFile(cfg);
  
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
  
  nPars_amp = CountUnique(params_amp);
  nPars_t0 = CountUnique(params_t0);
  nPars_Rq = CountUnique(params_Rq);
  nPars_BW = CountUnique(params_BW);
  
  std::cout << "-------- FIT PARAMETERS --------" << std::endl;
  for(unsigned int iRun = 0; iRun < nRuns; ++iRun)
  {
    if( params_amp.at(iRun) >= 0 ) parIndex[0+4*iRun] = params_amp.at(iRun);                                  else parIndex[0+4*iRun] = -1;
    if( params_t0.at(iRun)  >= 0 ) parIndex[1+4*iRun] = nPars_amp + params_t0.at(iRun);                       else parIndex[1+4*iRun] = -1;
    if( params_Rq.at(iRun)  >= 0 ) parIndex[2+4*iRun] = nPars_amp + nPars_t0 + params_Rq.at(iRun);            else parIndex[2+4*iRun] = -1;
    if( params_BW.at(iRun)  >= 0 ) parIndex[3+4*iRun] = nPars_amp + nPars_t0 + nPars_Rq + params_BW.at(iRun); else parIndex[3+4*iRun] = -1;
    std::cout << "iRun: " << iRun << "  parIndex[amp] = " << parIndex[0+4*iRun] << "   parIndex[t0] = " << parIndex[1+4*iRun] << "   parIndex[Rq] = " << parIndex[2+4*iRun] << "   parIndex[BW] = " << parIndex[3+4*iRun] << std::endl;
  }
  std::cout << "--------------------------------" << std::endl;
  
  return;
}



// **************************** **************************** ****************************
// Exponential or generic low-pass filter (time domain corresponding to L^-1[Ï„/(sÏ„+1)] 
// - can be used to set the kernel of a high-pass filter L^-1[sÏ„/(sÏ„+1)], 
// - which is given by a Î´(0) - the LP kernel
// **************************** **************************** ****************************
double funcLP(const double& xx, const double& tau)
{
  // Generic low-pass filter with cut-off time constant tau 
  // good for LYSO scintillation, RL filter, and delay line with RC-net
  double ff = 1/tau * exp(-xx/tau); 
  return ff;
}



// **************************** **************************** ****************************
// Exponential solution of the SiPM circuit from Abhinav K. Jha et.al., 2013
// IL = Ge * ( a1*exp(-xx/tcd1) + a2*exp(-xx/tcd2) + a3*exp(-xx/t_decay) + a4*exp(-xx/t_rise))
// **************************** **************************** ****************************
double SiPMPulseShape(const double& xx, const SiPMParams& sipmPars, const double& amp, const double& x0)
{
  double OV = sipmPars.OV;
  
  double Ge = (OV+0.25) * (sipmPars.Cq+sipmPars.Cd);
  
  double tmr = sipmPars.Rd * (sipmPars.Cq+sipmPars.Cd) * 1e9; // rise time, ns
  double tmd = sipmPars.Cd * (sipmPars.Rd*sipmPars.Rq) / (sipmPars.Rd + sipmPars.Rq) * 1e9; // decay time, ns
  
  double t1 = sipmPars.Rq * (sipmPars.Cq + sipmPars.Cd) * 1e9;
  double t2 = sipmPars.Rq * sipmPars.Cq * 1e9;
  double t3 = sipmPars.RL * sipmPars.Nc * sipmPars.Cd * 1e9;
  double t4 = sipmPars.RL * sipmPars.Cg * 1e9;
  double tcd1 = 0.5 * ( t1+t3+t4 + sqrt((t1+t3+t4)*(t1+t3+t4)-4*(t1*t4+t2*t3)) );
  double tcd2 = 0.5 * ( t1+t3+t4 - sqrt((t1+t3+t4)*(t1+t3+t4)-4*(t1*t4+t2*t3)) );
  
  double a1 = ( tcd1*(tcd1 - t2) ) / ( (tcd1-tcd2)*(tcd1-tmd)*(tcd1-tmr) );
  double a2 = ( tcd2*(tcd2 - t2) ) / ( (tcd2-tcd1)*(tcd2-tmd)*(tcd2-tmr) );
  double a3 = ( tmd*(tmd - t2) ) / ( (tmd-tcd1)*(tmd-tcd2)*(tmd-tmr) );
  double a4 = ( tmr*(tmr - t2) ) / ( (tmr-tcd1)*(tmr-tmd)*(tmr-tcd2) );

  double IL = 1E-15;
  if( xx >= x0 ) 
    IL = amp * Ge * 1e9 * ( a1*exp(-(xx-x0)/tcd1) + a2*exp(-(xx-x0)/tcd2) + a3*exp(-(xx-x0)/tmd) + a4*exp(-(xx-x0)/tmr));
  return IL;
}
