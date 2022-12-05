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



std::map<int,int> RemapPars(const std::vector<int>& params)
{
  // initialize remapping
  std::map<int,int> ret;
  for(unsigned int it = 0; it < 100; ++it)
    {
      ret[it] = -1;
    }
  
  int itPar = 0;
  for(unsigned int it = 0; it < params.size(); ++it)
    {
      int currPar = params.at(it);
      
      if( currPar == -1 ) ret[currPar] = -1;
      
      else
	{
	  if( ret[currPar] == -1 )
	    {
	      ret[currPar] = itPar;
	      ++itPar;
	    }
	  else
	    {
	      continue;
	    }
	}
    }
  
  return ret;
}



void GetSiPMParsFromCfg(char* cfg, std::vector<SiPMParams>& vec, std::vector<std::string>* runs)
{
  // parse the config file
  CfgManager opts;
  opts.ParseConfigFile(cfg);
  
  std::vector<std::string> SiPMList = opts.GetOpt<std::vector<std::string> >("Input.SiPMList");
  
  for(auto iSiPM : SiPMList)
  {
    SiPMParams pars;
    std::vector<float> vec_Npe = opts.GetOpt<std::vector<float> >(Form("%s.Npe",iSiPM.c_str()));
    std::vector<float> vec_OV  = opts.GetOpt<std::vector<float> >(Form("%s.OV",iSiPM.c_str()));
    std::vector<float> vec_Rq  = opts.GetOpt<std::vector<float> >(Form("%s.Rq",iSiPM.c_str()));
    std::vector<float> vec_Cq  = opts.GetOpt<std::vector<float> >(Form("%s.Cq",iSiPM.c_str()));
    std::vector<float> vec_Rd  = opts.GetOpt<std::vector<float> >(Form("%s.Rd",iSiPM.c_str()));
    std::vector<float> vec_Cd  = opts.GetOpt<std::vector<float> >(Form("%s.Cd",iSiPM.c_str()));
    std::vector<float> vec_Nc  = opts.GetOpt<std::vector<float> >(Form("%s.Nc",iSiPM.c_str()));
    std::vector<float> vec_Cg  = opts.GetOpt<std::vector<float> >(Form("%s.Cg",iSiPM.c_str()));
    std::vector<float> vec_RL  = opts.GetOpt<std::vector<float> >(Form("%s.RL",iSiPM.c_str()));
    std::vector<float> vec_RF  = opts.GetOpt<std::vector<float> >(Form("%s.RF",iSiPM.c_str()));
    pars.Npe   = vec_Npe.at(0);
    pars.OV    = vec_OV.at(0);
    pars.Rq    = vec_Rq.at(0);
    pars.Cq    = vec_Cq.at(0);
    pars.Rd    = vec_Rd.at(0);
    pars.Cd    = vec_Cd.at(0);
    pars.Nc    = vec_Nc.at(0);
    pars.Cg    = vec_Cg.at(0);
    pars.RL    = vec_RL.at(0);
    pars.RF    = vec_RF.at(0);
    if( vec_Rq.size() > 1 ) pars.RqErr = vec_Rq.at(1); else pars.RqErr = 0.;
    if( vec_Cq.size() > 1 ) pars.CqErr = vec_Cq.at(1); else pars.CqErr = 0.;
    if( vec_Rd.size() > 1 ) pars.RdErr = vec_Rd.at(1); else pars.RdErr = 0.;
    if( vec_Cd.size() > 1 ) pars.CdErr = vec_Cd.at(1); else pars.CdErr = 0.;
    if( vec_Cg.size() > 1 ) pars.CgErr = vec_Cg.at(1); else pars.CgErr = 0.;
    if( vec_RL.size() > 1 ) pars.RLErr = vec_RL.at(1); else pars.RLErr = 0.;
    if( vec_RF.size() > 1 ) pars.RFErr = vec_RF.at(1); else pars.RFErr = 0.;
    std::vector<std::string> tokens = opts.GetOpt<std::vector<std::string> >(Form("%s.title",iSiPM.c_str()));
    std::string title = "";
    for(auto token : tokens) title += token + " ";
    pars.title = title;
    pars.label = opts.GetOpt<std::string>(Form("%s.label",iSiPM.c_str()));
    vec.push_back(pars);
    
    if( runs )
      {
	std::string run = opts.GetOpt<std::string>(Form("%s.run",iSiPM.c_str()));
	(*runs).push_back(run);
      }
  }

  return;
}



void GetFitParsFromCfg(char* cfg, const int& nRuns, std::map<int,int>& parIndex,
                       int& nPars_amp, int& nPars_t0, int& nPars_Rq, int& nPars_Cd, int& nPars_BW, int& nPars_BW2)
{
  // parse the config file
  CfgManager opts;
  opts.ParseConfigFile(cfg);
  
  std::vector<int> params_amp;
  std::vector<int> params_t0;
  std::vector<int> params_Rq;
  std::vector<int> params_Cd;
  std::vector<int> params_BW;
  std::vector<int> params_BW2;
  
  std::vector<std::string> SiPMList = opts.GetOpt<std::vector<std::string> >("Input.SiPMList");
  for(auto iSiPM : SiPMList)
  {
    int param_amp = opts.GetOpt<int>(Form("%s.paramAmp",iSiPM.c_str())); params_amp.push_back( param_amp );
    int param_t0  = opts.GetOpt<int>(Form("%s.paramT0", iSiPM.c_str())); params_t0.push_back( param_t0 );
    int param_Rq  = opts.GetOpt<int>(Form("%s.paramRq", iSiPM.c_str())); params_Rq.push_back( param_Rq );
    int param_Cd  = opts.GetOpt<int>(Form("%s.paramCd", iSiPM.c_str())); params_Cd.push_back( param_Cd );
    int param_BW  = opts.GetOpt<int>(Form("%s.paramBW", iSiPM.c_str())); params_BW.push_back( param_BW );
    int param_BW2 = opts.GetOpt<int>(Form("%s.paramBW2",iSiPM.c_str())); params_BW2.push_back( param_BW2 );
  }
  
  nPars_amp = CountUnique(params_amp);
  nPars_t0  = CountUnique(params_t0);
  nPars_Rq  = CountUnique(params_Rq);
  nPars_Cd  = CountUnique(params_Cd);
  nPars_BW  = CountUnique(params_BW);
  nPars_BW2 = CountUnique(params_BW2);
  
  std::map<int,int> remap_amp = RemapPars(params_amp);
  std::map<int,int> remap_t0  = RemapPars(params_t0);
  std::map<int,int> remap_Rq  = RemapPars(params_Rq);
  std::map<int,int> remap_Cd  = RemapPars(params_Cd);
  std::map<int,int> remap_BW  = RemapPars(params_BW);
  std::map<int,int> remap_BW2 = RemapPars(params_BW2);
  
  std::cout << "-------- FIT PARAMETERS --------" << std::endl;
  std::cout << "nPars_amp: " << nPars_amp << std::endl;
  std::cout << "nPars_t0:  " << nPars_t0  << std::endl;
  std::cout << "nPars_Rq:  " << nPars_Rq  << std::endl;
  std::cout << "nPars_Cd:  " << nPars_Cd  << std::endl;
  std::cout << "nPars_BW:  " << nPars_BW  << std::endl;
  std::cout << "nPars_BW2: " << nPars_BW2 << std::endl;
  for(unsigned int iRun = 0; iRun < nRuns; ++iRun)
  {
    if( params_amp.at(iRun) >= 0 ) parIndex[0+6*iRun] = remap_amp[params_amp.at(iRun)];                                                         else parIndex[0+6*iRun] = -1;
    if( params_t0.at(iRun)  >= 0 ) parIndex[1+6*iRun] = nPars_amp + remap_t0[params_t0.at(iRun)];                                               else parIndex[1+6*iRun] = -1;
    if( params_Rq.at(iRun)  >= 0 ) parIndex[2+6*iRun] = nPars_amp + nPars_t0 + remap_Rq[params_Rq.at(iRun)];                                    else parIndex[2+6*iRun] = -1;
    if( params_Cd.at(iRun)  >= 0 ) parIndex[3+6*iRun] = nPars_amp + nPars_t0 + nPars_Rq + remap_Cd[params_Cd.at(iRun)];                         else parIndex[3+6*iRun] = -1;
    if( params_BW.at(iRun)  >= 0 ) parIndex[4+6*iRun] = nPars_amp + nPars_t0 + nPars_Rq + nPars_Cd + remap_BW[params_BW.at(iRun)];              else parIndex[4+6*iRun] = -1;
    if( params_BW2.at(iRun) >= 0 ) parIndex[5+6*iRun] = nPars_amp + nPars_t0 + nPars_Rq + nPars_Cd + nPars_BW + remap_BW2[params_BW2.at(iRun)]; else parIndex[5+6*iRun] = -1;
    std::cout << "iRun: " << iRun 
	      << "  parIndex[amp]  = " << parIndex[0+6*iRun] 
	      << "   parIndex[t0]  = " << parIndex[1+6*iRun] 
	      << "   parIndex[Rq]  = " << parIndex[2+6*iRun] 
	      << "   parIndex[Cd]  = " << parIndex[3+6*iRun] 
	      << "   parIndex[BW]  = " << parIndex[4+6*iRun] 
	      << "   parIndex[BW2] = " << parIndex[5+6*iRun] 
	      << std::endl;
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
// Double-exponential function to model scintillation (with rise and decay time)
// **************************** **************************** ****************************
double funcScint(const double& xx, const double& tau_r, const double& tau_d)
{
  double ff = ( exp(-1.*xx/tau_d) - exp(-1.*xx/tau_r) ) / (tau_d-tau_r);
  return ff;
}



// **************************** **************************** ****************************
// Exponential solution of the SiPM circuit from Abhinav K. Jha et.al., 2013
// IL = Ge * ( a1*exp(-xx/tcd1) + a2*exp(-xx/tcd2) + a3*exp(-xx/t_decay) + a4*exp(-xx/t_rise))
// **************************** **************************** ****************************
double SiPMPulseShape(const double& xx, const SiPMParams& sipmPars, const double& amp, const double& x0, const int& whichComponent)
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
    {
      IL = amp * Ge * 1e9 * ( a1*exp(-(xx-x0)/tcd1) + a2*exp(-(xx-x0)/tcd2) + a3*exp(-(xx-x0)/tmd) + a4*exp(-(xx-x0)/tmr) );
      
      if( whichComponent == 1 )
	IL = amp * Ge * 1e9 * a1*exp(-(xx-x0)/tcd1);
      
      if( whichComponent == 2 )
	IL = amp * Ge * 1e9 * a2*exp(-(xx-x0)/tcd2);
      
      if( whichComponent == 3 )
	IL = amp * Ge * 1e9 * a3*exp(-(xx-x0)/tmd);
      
      if( whichComponent == 4 )
	IL = amp * Ge * 1e9 * a4*exp(-(xx-x0)/tmr);
    }
  
  return IL;
}
