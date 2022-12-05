#include "interface/FitUtils.h"

#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;



std::pair<double,double> FindSmallestInterval(std::vector<double>* vals,
					      const double fraction, const bool verbosity)
{
  if( verbosity )
    std::cout << ">>>>>> FindSmallestInterval" << std::endl;
  
  
  std::sort(vals->begin(),vals->end());
  
  unsigned int nPoints = vals->size();
  unsigned int maxPoints = (unsigned int)(fraction * nPoints);
  
  double min = 0.;
  double max = 0.;
  unsigned int minPoint = 0;
  unsigned int maxPoint = 0;
  double delta = 999999.;
  for(unsigned int point = 0; point < nPoints-maxPoints; point++)
    {
      double tmpMin = vals -> at(point);
      double tmpMax;
      if( point+maxPoints > 0 ) tmpMax = vals -> at(point+maxPoints-1);
      else tmpMax = 1;
      
      if( tmpMax-tmpMin < delta )
	{
	  delta = tmpMax - tmpMin;
	  min = tmpMin;
	  max = tmpMax;
	  minPoint = point;
	  maxPoint = point + maxPoints - 1;
	}
    }
  
  
  std::pair<double,double> ret(min,max);
  
  return ret;
}
