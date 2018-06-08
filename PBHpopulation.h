#ifndef __PBHpopulation_H
#define __PBHpopulation_H

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#include <math.h>
#include <vector>
#include <map>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h> 
#include <gsl/gsl_vector.h>

#include "cgs.h"

class Input;
class randgen;

using namespace std;

// **********************************************************************************************************************

class PBHpopulation {

 public:

   PBHpopulation();
   PBHpopulation(Input*, double , double);
   PBHpopulation(Input*, double , double, double, int);
   vector<double> runSimulation(Input*);
   double gasDensity(double, double, double); 
   void columnDensityPrint(Input* );
   double gasDensitySpheres(double x, double y, double z);
   double gasDensityEllipses(double x, double y, double z);
   double nH2_Gal(double r, double z);
   double nHI_Gal(double Rkpc, double Zkpc);
   double nHII_Gal(double Rkpc,double Zkpc);
   ~PBHpopulation();
   void WriteMassToFile();
   void WriteToFile();
   void RadioWriteToFile();
   void pix2ang_ring(long,long,double*,double*);
   void ang2pix_ring(long,double,double,long*); 

   vector<double> retrievevDetRadio() {
        return vDetRadio;
   };
   vector<double> retrievevDetXray() {
        return vDetXray;
   };
  
   vector<vector<double>> retrieveVelRings() {
	return saveVVector;
   };
 
   inline long index(int il, int ib) { 
       return il*nb + ib;
    };
   
   inline long index_gasMap(int il, int ib) { 
       return il*nb_gas + ib;
    };
   
   inline long index_gasCube(int ix, int iy, int iz) { 
       return (ix*ny_gas + iy)*nz_gas + iz;
    };

   
   inline double sech(const double A) {return 2.*exp(A)/(exp(2.*A)+1.);}
   
   inline double constantMassUpToR(double R) {
	   
	   if (R > xBulge) R = xBulge;
	   double density = bulgeMass / ( 4./3. * M_PI * xBulge*xBulge*xBulge );	   
	   return density * ( 4./3. * M_PI * R*R*R );
   }
   inline double constantRingMass(double Rmin, double Rmax) {
	   
	   double m1 = constantMassUpToR(Rmax);
	   double m2 = constantMassUpToR(Rmin);
	   return m1-m2;
   }
   
   inline double NFWMassUpToR(double R, double rho0) {
	   
	   return 4.*M_PI*rho0*Rs*Rs*Rs*( log((Rs + R)/Rs) - R/(Rs + R) );	
   }
   inline double NFWRingMass(double Rmin, double Rmax, double rho0) {
	   
	   double m1 = NFWMassUpToR(Rmax,rho0);
	   double m2 = NFWMassUpToR(Rmin,rho0);
	   return m1-m2;
   }
   
   inline double BurkertMassUpToR(double R, double rho0) {
	   
	   return M_PI*rho0*Rs*Rs*Rs*( -2.*atan(R/Rs) + 2.*log(1. + R/Rs) + log(1. + R*R/(Rs*Rs)) );
   }
   inline double BurkertRingMass(double Rmin, double Rmax, double rho0) {
	   
	   double m1 = BurkertMassUpToR(Rmax,rho0);
	   double m2 = BurkertMassUpToR(Rmin,rho0);
	   return m1-m2;
   }
   
   
   inline double EddingtonLuminosity(double M) {
	   
	   return 1.26 * 1.e31 * 1.e7 * (M/mass_sun); //erg/s
   }
   
   inline double MdotEddington(double M) {
	   
	   return 1.26 * 1.e31 * 1.e7 * (M/mass_sun) /c_light/c_light; //g/s
   }

   
   //inline double XrayEfficiency_old(double Mdot) {  
	   //return efficiency_SgrA * (Mdot / Mdot_SgrA);	   
   //}

 protected:
   
   double lambda;
   double eta0;
   double PBHmass_parameter;
   double PBHmasses[100100]; // + 100 to be certain not to exceed array
   double radioFluxes[100100];
   double xrayFluxes[100100];
   vector<double> vDetRadio;
   vector<double> vDetXray;
   vector<vector<double>> saveVVector;    

   vector<double> l_vec;
   vector<double> b_vec;
   double deltal;
   double deltab;
   vector<double> l_vec_gas;
   vector<double> b_vec_gas;
   double deltal_gas;
   double deltab_gas;
   
   double myRho0;
   vector<double> R_vec;
   vector<double> M_vec;
   
   vector<double> PBHmap;
   vector<double> PBHRadioMap;
   
   vector<double> gasCube;
   vector<double> gasMap;
   
   vector<double> brightSourceNumber;
   
   //vector<double> R_vec;
   //vector<double> M_vec;
   double deltar;
    
   long resolution;    
   long nside;
   long npix;
        
   randgen* myrandom;
   
   int nPBH;
   int tid;   
   int it;

   vector<double> radiusVector;
   vector<double> velocityVector;
   map<double, vector<double> > velocityDistributions;
   //vector<double> velEddVector;
   vector< vector<double> > velEddMatrix;
      
};

// **********************************************************************************************************************


#endif
