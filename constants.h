#ifndef __CONSTANTS_H
#define __CONSTANTS_H

#include "math.h"
#include "cgs.h"

// For radio+X-ray: define Chandra and required MD

#define MAXWELLIANS
//#define DISTRIBUTIONS
//#define PESSIMISTIC
//#define EDDPHASESPACE

//#define SAVEVELIRING

//#define CONSERVATIVE
//#define NUSTAR
#define CHANDRA
//#define SKA

//#define DELTA
#define LOGNORMAL
//#define POWERLAW

//#define NEWLAMBDA

//#define BURKERT

//#define StromgrenCheck

//
//const std::string mainFile = "prova.dat";
//const std::string mainFile = "test_200kms.dat";

//const std::string mainFile = "NsourcesVLA_f1_MAXWELLIANS_lambda0.01_30Msun.dat";
//const std::string mainFile = "NsourcesVLA_f1_MAXWELLIANS_lambda0.01_FINAL.dat";
//const std::string mainFile = "NsourcesVLA_f0.1_MAXWELLIANS_lambda1_FINAL.dat";
//const std::string mainFile = "NsourcesVLA_f1_MAXWELLIANS_lambda0.1_LARGE.dat";
const std::string mainFile = "NsourcesSKA_PhotonFlux_.txt";
const double lambdaRicotti = .02;//.01;
const int    nIterations = 1000;// 30;

//const std::string mainFile = "NsourcesNuStar_f1_MAXWELLIANS_PESSIMISTIC_lambda1_PROVA.dat";

const double gasV = 50.*km/s;

const std::string histFileBase = "prova_";

const std::string vDistributionFilename = "FvAsFctOfRadius1t39.txt"; 
//const std::string velocitiesMaxwellianFilename = "Maxwellians_MODIFIED.txt";
const std::string velocitiesMaxwellianFilename = "Maxwellians.txt";
const std::string velocitiesEddingtonFilename = "fvrEdd_grid_FINAL.txt";

const int file_downsampling = 1000;
const int Ndebug_bigcounter = 1000000;
const int Ndebug = 1.e10;
const int NringDebug = 100;

const double DMfraction = .1;

const double vmin = 0.*km/s;
const double vmax = 1000.*km/s;
const double vminEdd = 0.*km/s;
const double vmaxEdd = 775.*km/s;
const int nvel = 1000;
const double dvEdd = (vmaxEdd - vminEdd)/(nvel-1); // -1 to make sure it includes
//const int saveV_iRing = 0;

const double soundSpeed = 1.0*km/s;
const double soundSpeedIonized = 10.0*km/s;
const double Delta_T = 0.5*pow(soundSpeedIonized/soundSpeed,2.);
const double machD = sqrt(2*Delta_T) * ( 1 - sqrt(1 - 0.5/Delta_T) );
const double machR = sqrt(2*Delta_T) * ( 1 + sqrt(1 - 0.5/Delta_T) );


const double milkyWayMass = 1.e12*mass_sun;
const double bulgeMass    = 0.4*1.84e10*mass_sun; //arXiv:1502.00633 //0.4 is the portion of the mass due to DM

const double xBulge = 2.2*kpc;
const double yBulge = 1.4*kpc;
const double zBulge = 1.2*kpc;
const double rBulge = 2.*kpc;

const double Rs = 16.*kpc; //NFW scale radius
const double rhoSun = 1.4e-2 * mass_sun/pc3; //http://arxiv.org/pdf/1304.5127.pdf

const double minPBHMass = 20.*mass_sun;//8.*mass_sun;
//const double minPBHMass =   10.*mass_sun;
const double maxPBHMass = 100.*mass_sun;// 200.*mass_sun;
const int numMass = 10;// 10;
//const int numMass = 3;

const double sigma_LN = 1.00;
const double Mstar = pow(10.,1.);

#ifdef POWERLAW
const double minPBHMassPL = 0.1*Mstar*mass_sun;
const double maxPBHMassPL = Mstar*mass_sun;
#endif

const double etaFender = 0.1; 

const double Xslope = 1.6;

const double bandWidth = 18.e6; //Hz
const double obsFrequency = 1.4e9; //Hz //run 12.09.2016
const double Jy_erg_conversion =  1.e-23; //erg/cm2/s/Hz/Jy

const double l = 0.;//  0.113*deg;
const double b = 0.;//  -1.424*deg;

#ifdef NUSTAR
	const double lmin =  -0.9*deg; //l - angRadius;
	const double lmax =  0.3*deg; //l + angRadius;
	const double bmin =  -0.1*deg; //b - angRadius;
	const double bmax = 0.4*deg; //b + angRadius;
#endif
	
#ifndef NUSTAR
	const double lmin = -.9*deg;//  -0.5*deg; //-0.5 //l - angRadius;
	const double lmax = .7*deg;// 0.5*deg; //0.5 //l + angRadius;
	const double bmin = -0.3*deg; //b - angRadius;
	const double bmax =  0.3*deg; //b + angRadius;
#endif

const double nl = 300; // 80;
const double nb = 300; //  32;

const double lmin_gas = -10.*deg; //l - angRadius;
const double lmax_gas =  10.*deg; //l + angRadius;
const double bmin_gas = -10.*deg; //b - angRadius;
const double bmax_gas =  10.*deg; //b + angRadius;
const int nl_gas = 200;
const int nb_gas = 200;
const double xmax_gas = 5.*kpc;
const double ymax_gas = 5.*kpc;
const double zmax_gas = 5.*kpc;
const int nx_gas = 500;
const int ny_gas = 500;
const int nz_gas = 500;

const double xsun = 8.3 * kpc; //kpc
const double rsun = xsun;
const double ysun = 0.;
const double zsun = 0.;

const double rho_s = rhoSun *  (1 + pow(rsun/Rs, 2.)) * (1. + rsun/Rs);

const double XPhotonFluxChandraThreshold = 2.e-6/cm/cm/s;  //1.e-26*erg/cm/cm/s;
const double Nustar_threshold =  8.e32*erg/s;
const double Chandra_threshold =  2.2e32*erg/s;
const double RadioThreshold5GHz_mJy = 2.47; //mJy
const double RadioThreshold1GHz_mJy = 1.60;// 1.60;// 0.012; //1.60;// 0.012; // 1.60; //mJy
const double RadioThresholdSKA_mJy = 85.e-6; // 85.e-6;// 1.e-3;//1.e-3; //mJy
//1000 h -> 85 nJ  = 85e-6 mJ
//   1 h ->  2.7muJ = 2.7e-3 mJ

const double photonEnergy = 3.*keV; //CHECK

const double rmin = 0.*kpc;
const double rmax = 5.*kpc;
const int nr = 500;

const double X_CO = 1.9;

const int randgen_numsteps = 10;

/*
const double GeV_erg = 0.00160217646;
const double erg_GeV = 624.15; 
const double cLight = 2.99792458e10;               // speed of light in vacuum cm*s^-1
const double yr = 31556925.2;             // tropical year in seconds
const double sec_yr = 1./yr;          // to pass from seconds to year
const double pc = 3.0856775807e18 ;        // parsec in cm
const double kpc = 3.0856775807e21;        // kpc in cm
*/

//const int nE = 31 ;
//const double Emin = 1.*keV; //GeV
//const double Emax = 1000.*keV; //GeV
//const double E_factor = pow((Emax/Emin),1./(nE-1));

//const double angRadius = 2.56*arcmin; // not used

//const double XluminosityThreshold = 2.e-6/cm/cm/s; // photons/cm2/s not used
//const double radioLuminosityThreshold = 2.e-6/cm/cm/s; // photons/cm2/s not used

//const double rho0 = 1.4e-2 * mass_sun/pc3; //http://arxiv.org/pdf/1304.5127.pdf

//const double vAverage = 200.*km/s * (2.*sqrt(2./M_PI));  //12.09.2016
//const double sigmav   = 200.*km/s ; // vAverage / (2.*sqrt(2./M_PI));

//const double LIGOmass = 30.*mass_sun;

//const double gasDensity = 10.*mass_proton/cm3;

//const double efficiency_SgrA = 1.e-8; //Lx = 10^33 erg/s
//const double Lx_SgrA = 10.*1.e33; //erg/s //FLARING
//const double MdotBondi_SgrA = 1.e-6*mass_sun/year; //10^41 erg/s
//const double XemissionEfficiency = 1.e-5;// 1.e-5; // L/(Mdot c^2)

//const double measuredFlux = 0.;// 8.6e-11 / squareDeg; //erg/cm2/s/srs

#endif
