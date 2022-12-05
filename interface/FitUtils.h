#ifndef FITUTILS_H
#define FITUTILS_H

#include <vector>
#include <iostream>


using namespace std;

std::pair<double,double> FindSmallestInterval(std::vector<double>* vals,                        
					      const double fraction = 0.68, const bool verbosity = false);


#endif
