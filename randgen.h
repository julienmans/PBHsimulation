#ifndef __RANDGEN_H
#define __RANDGEN_H

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <math.h>
#include <vector>
#include <map>
#include <random>

class Input;

using namespace std;

class randgen {

 public:
  randgen();
  randgen(Input*, double, vector<double>, vector<double>, map<double, vector<double> >);
  ~randgen();
  double NFW(double);
  double constantDistribution(double);
  //double luminosityFunction(double);
  double MaxwellBoltzmann(double,double);
  double random_arbitrary_function(string, double, double,double,double) ;
  double genericVDistribution(double, double);

 protected:
  Input* inp;
  double rho0;
  vector<double> radiusVector;
  vector<double> velocityVector;
  map<double, vector<double> > velocityDistributions;
  unsigned short seed1; 
  int tid;
  // using mersene twister algorithm
  std::mt19937_64 generator;
  // Using a uniform random dist from standard library
  std::uniform_real_distribution<double> unif;


};


#endif
