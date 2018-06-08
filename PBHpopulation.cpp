#include "constants.h"
#include "PBHpopulation.h"
#include "input.h"
#include "randgen.h"
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h> 
#include <gsl/gsl_vector.h>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <ctime>

#include <omp.h> // For parallelization (to retrieve the thread number)
#include <chrono> // For seed random variable generator
#include <random>
//#include <limits> // For max precision

#define DEBUG 0
#define NUMTHREADS 4

using namespace std;


PBHpopulation::PBHpopulation() { 
	cout << "I am the default costructor. I do nothing" << endl;
}

double PBHpopulation::nH2_Gal(double r, double z) {

	int i;
	double nH2_ = 0.0, fR,fZ0,fZh;                                              // [B88]/Table 3
	double R[18] ={ 0.00, 2.25, 2.75, 3.25, 3.75, 4.25, 4.75, 5.25, 5.75,       // kpc, col.1
			6.25, 6.75, 7.25, 7.75, 8.25, 8.75, 9.25, 9.75,10.25},
			Y[18] ={ 0.00,  1.5,  3.3,  5.8,  5.5,  8.4,  9.0,  9.6,  8.6,       // CO, K km s^-1
					9.1,  7.9,  9.2,  7.7,  5.0,  3.6,  4.8,  1.7,  0.0},// (col.4)
					Z0[18]={0.039,0.039,0.036,0.000,-.008,0.001,-.010,-.001,-.004,       // kpc, col.7
							-.019,-.022,-.014,-.009,-.004,0.013,-.004,-.020,-.020},
							Zh[18]={0.077,0.077,0.080,0.061,0.065,0.071,0.072,0.082,0.083,       // kpc, col.10
									0.073,0.063,0.058,0.072,0.080,0.066,0.023,0.147,0.147};

	double H2toCO = 1.e20;              // [SM96]

	if(r > R[17]) return nH2_*mass_proton;
	for (i=0; i<17; i++)  if(R[i] <= r && r <= R[i+1])  break;

	fR =  Y[i] + ( Y[i+1] - Y[i])/(R[i+1] - R[i])*(r - R[i]); 
	fZ0= Z0[i] + (Z0[i+1] -Z0[i])/(R[i+1] - R[i])*(r - R[i]);
	fZh= Zh[i] + (Zh[i+1] -Zh[i])/(R[i+1] - R[i])*(r - R[i]);
	nH2_ =  fR * exp( -log(2.)*pow( (z-fZ0)/fZh, 2 ) )  * H2toCO/kpc * 2.*mass_proton;

	return nH2_< 0. ? 0.: nH2_;
}


double PBHpopulation::nHI_Gal(double Rkpc, double Zkpc)
{
	int i;                                                             // Table 1 [GB76]
	double R[30] ={ 0.0, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0,  // kpc, col.1
			6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5,10.0,10.5,11.0,
			11.5,12.0,12.5,13.0,13.5,14.0,14.5,15.0,15.5,16.0},
			Y[30] ={ .10, .13, .14, .16, .19, .25, .30, .33, .32, .31,  // nHI, cm^-3
					.30, .37, .38, .36, .32, .29, .38, .40, .25, .23,  // (col.3)
					.32, .36, .32, .25, .16, .10, .09, .08, .06, .00};
	double fR, fZ,fZ1=0.,fZ2=0., R1,R2=R[29], Y1,Y2=Y[29];
	double nGB =0.33, nDL =0.57;       // cm^-3, disk density @ 4-8 kpc; [GB76], [DL90]
	double A1=0.395,     z1=0.212/2.,  // cm^-3, kpc; Z-distribution parameters from [DL90]
			A2=0.107,     z2=0.530/2.,
			B =0.064,     zh=0.403;

	for (i=0; i<29; i++)  if(R[i] <= Rkpc && Rkpc <= R[i+1])  break;

	R1 = (R[i]+R[i+1])/2;   Y1 = Y[i];
	if(Rkpc < R1)
	{  
		if(i> 0)    { R2 = (R[i-1]+R[i])/2;   Y2 = Y[i-1]; }
		else        { R2 = R[0];              Y2 = Y[0];   }
	}
	else  if(i<28) { R2 = (R[i+1]+R[i+2])/2; Y2 = Y[i+1]; }

	fR = Y1 +(Y2 -Y1)/(R2 -R1)*(Rkpc -R1);                             // interpolation in R

	R2 = (R[28] +R[29]) /2;
	if(Rkpc > R2) fR = Y[28]*exp(-(Rkpc-R2)/3);                        // extrapolation in R

	// calculation of Z-dependence
	if(Rkpc <10.)                                                      // [DL90]
		fZ1 =A1*exp(-log(2.)*pow(Zkpc/z1,2))+A2*exp(-log(2.)*pow(Zkpc/z2,2))+B*exp(-fabs(Zkpc)/zh);
	if(Rkpc > 8.) 
		fZ2=nDL*exp(-pow(Zkpc /(0.0523*exp(0.11*Rkpc)), 2)); // [C86] 

	if(Rkpc <= 8.) 
		fZ = fZ1;
	else {   
		if(Rkpc >=10.) 
			fZ = fZ2;
		else 
			fZ = fZ1 +(fZ2 -fZ1)/2.*(Rkpc -8.);                       // interp. [DL90] & [C86]
	}
	return mass_proton*(fZ *fR/nGB);
}

double PBHpopulation::nHII_Gal(double Rkpc,double Zkpc)
{

	double fne1=0.025, H1=1.00, A1=20.0;
	double fne2=0.200, H2=0.15, A2= 2.0;
	double R2=4.0;
	double ne1 = fne1 * exp(-fabs(Zkpc)/H1) * exp (-pow( Rkpc    /A1, 2));
	double ne2 = fne2 * exp(-fabs(Zkpc)/H2) * exp (-pow((Rkpc-R2)/A2, 2));
	return mass_proton*(ne1+ne2);
}

double PBHpopulation::gasDensitySpheres(double x, double y, double z) {

	double r = sqrt(x*x+y*y+z*z);
	double H2_CMZ = 0.;
	double HI_CMZ  = 0.;

	if (r > 0.*kpc && r <= 1.3*kpc) {

		if (r > 0.*kpc && r <= 0.450*kpc) {
			double volume = 4/3 * 3.1415926 * (pow(0.450, 3));
			double mass_H2_CMZ = 2.4e7; //Msun mass H2
			H2_CMZ += mass_H2_CMZ/volume * mass_sun /(1e9 * pc * pc * pc);
		}
		else if (r > 0.45*kpc && r <= 1.3*kpc) {
			double volume = 4/3 * 3.1415926 * (pow(1.3, 3) - pow(0.450, 3));
			double mass_H2_CMZ = 2.9e7; //Msun mass H2
			H2_CMZ += mass_H2_CMZ/volume * mass_sun /(1e9 * pc * pc * pc);
		}
		double volume = 4/3 * 3.1415926 * (pow(1.3, 3));
		double mass_HI_CMZ = 1.1e7;
		HI_CMZ += mass_HI_CMZ/volume * mass_sun /(1e9 * pc * pc * pc);
		return (H2_CMZ + HI_CMZ);
	}
	else if (r > 1.3*kpc && r <= 3*kpc) {
		double volume = 4/3 * 3.1415926 * (pow(3., 3) - pow(1.3, 3));
		double mass_H2_CMZ = 1.3e7; //Msun mass H2
		H2_CMZ += mass_H2_CMZ/volume * mass_sun /(1e9 * pc * pc * pc);
		return H2_CMZ;					   
	}
	else {
		return 1.*mass_proton/cm3;
	}
	return 0.; 

}

double PBHpopulation::gasDensityEllipses(double x, double y, double z) {

	double r = sqrt(x*x + y*y);
	double H2_CMZ = 0.;
	double HI_CMZ  = 0.;
	double zheight = 0.100; // kpc
	if (r > 0.*kpc && r <= 1.3*kpc && z <= zheight*kpc) {

		if (r > 0.*kpc && r <= 0.450*kpc) {
			double volume = 3.1415926 * zheight * (0.450 * 0.450);
			double mass_H2_CMZ = 2.4e7 ;
			H2_CMZ += mass_H2_CMZ/volume * mass_sun /(1e9 * pc * pc * pc);
		}
		else if (r > 0.450*kpc && r <= 1.3*kpc) {
			double volume = 3.1415926 * zheight * (1.3 * 1.3 - 0.450 * 0.450);;
			double mass_H2_CMZ = 2.9e7;
			H2_CMZ += mass_H2_CMZ/volume * mass_sun /(1e9 * pc * pc * pc);
		}
		double volume = 3.1415926 * zheight * (1.3 * 1.3);; // Volume HI, z = 0.050 pc
		double mass_HI_CMZ = 1.1e7;
		HI_CMZ += mass_HI_CMZ/volume * mass_sun /(1e9 * pc * pc * pc);
		return (H2_CMZ + HI_CMZ);
	}
	else if (r > 1.3*kpc && r <= 3.*kpc && z <= zheight*kpc) {
		double volume = 3.1415926 * zheight * (3. * 3. - 1.3 * 1.3);;
		double mass_H2_CMZ = 1.3e7 ;
		H2_CMZ += mass_H2_CMZ/volume * mass_sun /(1e9 * pc * pc * pc);
		return H2_CMZ;					   
	}

	else {
		return 1.*mass_proton/cm3;
	}
	return 0.; 

}

double PBHpopulation::gasDensity(double x, double y, double z) {

	//cout << "gasDensity" << endl;
	//cout << "x = " << x/kpc << " y = " << y/kpc << " z = " << z/kpc << " ";
	double r = sqrt(x*x+y*y+z*z);
	double H2_CMZ, HI_CMZ, HII_CMZ, H2_GBdisk, HI_GBdisk;
	double xc, yc, thetac, xb, yb, x2, y2, Xc, Lc, Hc2, HcI;
	double alpha1, beta1, thetad1, Xd, Ld, Hd, HdI, x3, y3, z3;
	if (r < 2.0*kpc) {
		//CMZ contribution ********************************************************************************
		xc = -50.*pc; yc = 50.*pc;
		thetac = 70. * degreeToRad;
		xb = x - xc; yb = y - yc;
		x2 = (xb)*cos(thetac) + (yb)*sin(thetac);
		y2 = - (xb)*sin(thetac) + (yb)*cos(thetac);
		Xc = 125.*pc;
		Lc = 137.*pc;
		Hc2 = 18.*pc;
		HcI = 54.*pc;       
		H2_CMZ = 150.*mass_proton/cm3 * exp( - pow( ((sqrt( pow(x2,2) + pow((2.5*y2),2) ) - Xc) / (Lc)) ,4) ) * exp ( - pow((z/Hc2),2)  );
		HI_CMZ = 8.8 *mass_proton/cm3 * exp( -  pow( ((sqrt( pow(x2,2) + pow((2.5*y2),2) ) - Xc) / (Lc)) ,4) ) * exp ( - pow((z/HcI),2)  ); 
		double y_3 = -10.*pc; //pc
		double z_3 = -20.*pc;
		double L3 = 145.*pc;
		double H3 = 26.*pc;
		double L2 = 3700.*pc;
		double H2 = 140.*pc;
		double L1 = 17000.*pc;
		double H1 = 950.*pc;    
		double temp = 8.*0.005*mass_proton/cm3*cos(M_PI*r/(2*L1))*sech(z/H1)*sech(z/H1);        
		HII_CMZ = 8.0*mass_proton/cm3 * (exp(-(x*x+(y-y_3)*(y-y_3))/(L3*L3)) * exp(-((z-z_3)*(z-z_3))/(H3*H3))) + 8.*0.009*mass_proton/cm3 * exp(-((r-L2)*(r-L2))/(L2*L2/4.))*sech(z/H2)*sech(z/H2) + temp;
		//GB disk contribution *****************************************************************************       
		alpha1 = 13.5 * degreeToRad;
		beta1 = 20. * degreeToRad;
		thetad1 = 48.5 * degreeToRad;        
		Xd = 1200.*pc; //pc
		Ld =  438.*pc;
		Hd =   42.*pc;
		HdI = 120.*pc;        
		x3 = x*cos(beta1)*cos(thetad1) - y*(sin(alpha1)*sin(beta1)*cos(thetad1) - cos(alpha1)*sin(thetad1)) - z*(cos(alpha1)*sin(beta1)*cos(thetad1) + sin(alpha1)*sin(thetad1));        
		y3 = -x*cos(beta1)*cos(thetad1) + y*(sin(alpha1)*sin(beta1)*sin(thetad1) + cos(alpha1)*cos(thetad1)) + z*(cos(alpha1)*sin(beta1)*sin(thetad1) - sin(alpha1)*cos(thetad1));        
		z3 = x*sin(beta1) + y*sin(alpha1)*cos(beta1) + z*cos(alpha1)*cos(beta1);        
		H2_GBdisk = 4.8*mass_proton/cm3  * exp(- pow( ((sqrt( pow(x3,2) + pow((3.1*y3),2) ) - Xd ) / (Ld)) ,4) ) * exp(- pow((z3/Hd),2));        
		HI_GBdisk = 0.34*mass_proton/cm3 * exp(- pow( ((sqrt( pow(x3,2) + pow((3.1*y3),2) ) - Xd ) / (Ld)) ,4) ) * exp(-pow((z3/HdI),2));
		//cout << (2.*H2_CMZ + HI_CMZ + HII_CMZ + H2_GBdisk + HI_GBdisk)/mass_proton << endl;
		return 2.*H2_CMZ + HI_CMZ + HII_CMZ + H2_GBdisk + HI_GBdisk;
	}
	else {
		//cout << (2.*X_CO*nH2_Gal(r/kpc,z/kpc) + nHI_Gal(r/kpc,z/kpc) + nHII_Gal(r/kpc,z/kpc))/mass_proton << endl;
		return  2.*X_CO*nH2_Gal(r/kpc,z/kpc) + nHI_Gal(r/kpc,z/kpc) + nHII_Gal(r/kpc,z/kpc);
	}
}


void PBHpopulation::columnDensityPrint(Input* inp)
{

	deltal_gas = (lmax_gas - lmin_gas)/(nl_gas-1);
	deltab_gas = (bmax_gas - bmin_gas)/(nb_gas-1);
	for (int il=0; il<nl_gas; il++)
		l_vec_gas.push_back(lmin_gas + il*deltal_gas);
	for (int ib=0; ib<nb_gas; ib++)
		b_vec_gas.push_back(bmin_gas + ib*deltab_gas);
	for (int il=0; il<nl_gas; il++)
		for (int ib=0; ib<nb_gas; ib++) 
			gasMap.push_back(0.);

	double deltax_gas = (2.*xmax_gas)/(nx_gas-1);
	double deltay_gas = (2.*ymax_gas)/(ny_gas-1);
	double deltaz_gas = (2.*zmax_gas)/(nz_gas-1);
	for (int ix = 0; ix < nx_gas; ix++) {
		if (ix%100==0)
			cout << 100.*double(ix)/double(nx_gas)<< "%" << endl;
		for (int iy = 0; iy < ny_gas; iy++) {
			for (int iz = 0; iz < nz_gas; iz++) {
				//cout << ix << " " << iy << " " << iz << endl;
				double x_ = (-xmax_gas) + deltax_gas*ix;
				double y_ = (-ymax_gas) + deltay_gas*iy;
				double z_ = (-zmax_gas) + deltaz_gas*iz;
				double value = this->gasDensity(x_,y_,z_);
				gasCube.push_back( value );
				//cout << x_/kpc << " " << y_/kpc << " " << z_/kpc << " " << value/mass_proton << endl << endl;
			}
		}
	}

	double limitR = 8.5*kpc;
	double limitZ = 1.*kpc;
	double ds = 0.010*kpc;

	for (int il=0; il<nl_gas; il++)
		for (int ib=0; ib<nb_gas; ib++) 
			gasMap.push_back(0.);

	for (int il=0; il<nl_gas; il++) {
		if (il%10==0)
			cout << 100.*double(il)/double(nl_gas)<< "%" << endl;
		for (int ib=0; ib<nb_gas; ib++) { 

			double l = l_vec_gas[il];
			double b = b_vec_gas[ib];
			double sinl = sin(l);
			double sinb = sin(b);
			double cosl = cos(l);
			double cosb = cos(b);
			double currentRadius, currentX, currentY, currentZ;
			double d = 0.*kpc;
			currentRadius = 0.*kpc;
			currentZ = 0.*kpc;

			//line-of-sight integral
			if ((il==nl_gas/2) && (ib==nb/2)) {
				cout << "Line of sight integral; l = " << l*180./M_PI << " b =  " << b*180./M_PI << endl; 
			}
			
			int counter = 0;

			while ( (currentRadius < limitR) && (fabs(currentZ) < limitZ) ) {
				
				d += ds;

				currentZ=d*sinb;
				currentRadius=(xsun>0) ? xsun*xsun+pow(d*cosb,2)-2.0*xsun*d*cosb*cosl : xsun*xsun+pow(d*cosb,2)+2.0*xsun*d*cosb*cosl; // Galactocentric distance of point
				if(currentRadius<1.e-100*kpc) 
					currentRadius = 1.e-100*kpc;
				currentRadius = sqrt(currentRadius);
				double phi = atan2(d*cosb*sinl, xsun - d*cosb*cosl);
				currentX = currentRadius*cos(phi); 
				currentY = currentRadius*sin(phi); 

				int ixgas = int(((currentX-(-xmax_gas))/deltax_gas));
				if (ixgas > nx_gas-1) ixgas = nx_gas-1;
				if (ixgas < 0) ixgas = 0;
				//
				int iygas = int(((currentY-(-ymax_gas))/deltay_gas));
				if (iygas > ny_gas-1) iygas = ny_gas-1;
				if (iygas < 0) iygas = 0;
				//
				int izgas = int(((currentZ-(-zmax_gas))/deltaz_gas));
				if (izgas > nz_gas-1) izgas = nz_gas-1;
				if (izgas < 0) izgas = 0;
				
				double value = gasCube[index_gasCube(ixgas,iygas,izgas)];
				if ((il==nl_gas/2) && (ib==nb/2) && (counter%100==0)) {
					cout << "d = " << d/kpc << endl;
					cout << "x = " << currentX/kpc << " y = " << currentY/kpc << " z = " << currentZ/kpc << " " << value/mass_proton << endl;
				}
				
				gasMap[index_gasMap(il,ib)] += value*ds;
				counter++;

			}
		}
	}

	ofstream gasTable;
	gasTable.open("gasTable.txt");
	for (int il=0; il<nl_gas; il++)
		for (int ib=0; ib<nb_gas; ib++) 
			gasTable << l_vec_gas[il]*180./M_PI << "\t" << b_vec_gas[ib]*180./M_PI << "\t" << gasMap[index_gasMap(il,ib)] << "\n";
	gasTable.close();

}


// constructor --------------------------------------------------------------------------------------------------------------------------------

PBHpopulation::PBHpopulation(Input* inp, double lambda_, double eta_)
{

	/*cout << "-------------------------------------" << endl;
  cout << "This is the PBHpopulation constructor" << endl;
  cout << "lambda = " << lambda_ << endl;
  cout << "eta = " << eta_ << endl;
  cout << "-------------------------------------" << endl << endl;
  lambda = lambda_;
  eta0 = eta_;

  nPBH = 0;

  myRho0 = bulgeMass/( 4.*M_PI*Rs*Rs*Rs*( log((Rs + rBulge)/Rs) - rBulge/(Rs + rBulge) ) );
  cout << "rho0 [Msun/pc3] = " << myRho0/(mass_sun/(pc3)) << endl;

  myrandom  = new randgen(inp,myRho0,radiusVector,velocityVector,velocityDistributions);

  deltar = (rmax-rmin)/(nr-1);
  for (int ir=0; ir < nr; ir++) {
	  R_vec.push_back(rmin + ir*deltar);
	  M_vec.push_back(0.);
  } 

  resolution = inp->resolution;
  nside = pow(2,resolution);
  npix = 12*nside*nside;
  deltal = (lmax-lmin)/(nl-1);
  deltab = (bmax-bmin)/(nb-1);
  cout << "l_vec " << endl;
  for (int il=0; il<nl; il++){
	  l_vec.push_back(lmin + il*deltal);
	  brightSourceNumber.push_back(0.);
	  cout << l_vec.back()/deg << " " ;
  }
  cout << endl;
  cout << "b_vec " << endl;
  for (int ib=0; ib<nb; ib++){
	  b_vec.push_back(bmin + ib*deltab);
	  cout << b_vec.back()/deg << " " ;
  }
  cout << endl;
  cout << "PBHmap " << endl;
  for (int il=0; il<nl; il++)
	  for (int ib=0; ib<nb; ib++) {
		  PBHmap.push_back(0.);
		  PBHRadioMap.push_back(0.);
	  }

  cout << "Grid correctly initialized" << endl;*/

}


PBHpopulation::PBHpopulation(Input* inp, double lambda_, double eta_, double mass_, int it_)
{
	tid = omp_get_thread_num();
	if (tid == 0) {
	  cout << "----------------------------------------------------" << endl;
	  cout << "This is the PBHpopulation constructor; mass = " << mass_/mass_sun << " Msun " << endl;
	  cout << "----------------------------------------------------" << endl;
	}
	lambda = lambda_;
	eta0 = eta_;
	PBHmass_parameter = mass_/mass_sun;
	it = it_;

	nPBH = 0;

	myRho0 = bulgeMass/( 4.*M_PI*Rs*Rs*Rs*( log((Rs + rBulge)/Rs) - rBulge/(Rs + rBulge) ) );

#ifdef BURKERT
	myRho0 = rho_s;// bulgeMass/(  M_PI*Rs*Rs*Rs*( -2.*atan(rBulge/Rs) + 2.*log(1. + rBulge/Rs) + log(1. + rBulge*rBulge/(Rs*Rs)) ) );
#endif	
	
	//cout << "rho0 [Msun/pc3] = " << myRho0/(mass_sun/(pc3)) << endl;

	deltar = (rmax-rmin)/(nr-1);
	for (int ir=0; ir < nr; ir++) {
		R_vec.push_back(rmin + ir*deltar);
		M_vec.push_back(0.);
	} 

	resolution = inp->resolution;
	nside = pow(2,resolution);
	npix = 12*nside*nside;
	deltal = (lmax-lmin)/(nl-1);
	deltab = (bmax-bmin)/(nb-1);
	//cout << "l_vec " << endl;
	for (int il=0; il<nl; il++){
		l_vec.push_back(lmin + il*deltal);
		brightSourceNumber.push_back(0.);
		//cout << l_vec.back()/deg << " " ;
	}
	//cout << endl;
	//cout << "b_vec " << endl;
	for (int ib=0; ib<nb; ib++){
		b_vec.push_back(bmin + ib*deltab);
		//cout << b_vec.back()/deg << " " ;
	}
	//cout << endl;
	//cout << "PBHmap " << endl;
	for (int il=0; il<nl; il++)
		for (int ib=0; ib<nb; ib++) {
			PBHmap.push_back(0.);
			PBHRadioMap.push_back(0.);
		}
	if (tid == 0) {
	cout << "Grid correctly initialized!" << endl;
	}
#ifdef DISTRIBUTIONS    
	if (tid == 0) {
	cout << "Reading velocity distribution... " << endl;
	}
	std::ifstream vFile(vDistributionFilename);
	std::string line;
	int counter = 0;
	while (std::getline(vFile, line))
	{
		std::istringstream iss(line);
		if (counter == 0) {
			double temp;
			while (iss >> temp){
				velocityVector.push_back(temp*km/s);
				//cout << temp << " km/s " << endl;
			}
		}
		else {
			int counter2 = 0;
			double temp;
			while (iss >> temp) {
				if (counter2==0) {
					//cout << temp << " kpc " << endl;
					radiusVector.push_back(temp*kpc);
				}
				else {
					//cout << " v = " << velocityVector[counter2-1]/(km/s) << " " << temp << endl;
					velocityDistributions[counter-1].push_back(temp);
				}
				counter2++;
			}
			//cout << endl;
		}
		counter++;
	}
	if (tid == 0) {
	cout << "...done" << endl;
	}
#endif

#ifdef MAXWELLIANS  
	if (tid == 0) {
	cout << "Reading velocity distribution... " << endl;
	}
	std::ifstream vFile(velocitiesMaxwellianFilename.c_str());
	std::string line;
	int counter = 0;
	while (std::getline(vFile, line))
	{
		std::istringstream iss(line);
		double temp1, temp2, temp3;
		iss >> temp1 >> temp2 >> temp3;
		radiusVector.push_back(temp2*kpc);
		velocityVector.push_back(temp3*km/s);
		velocityDistributions[0].push_back(0.);
		counter++;
	}
#endif   

#ifdef EDDPHASESPACE
	if (tid == 0) {
		cout << "Reading velocity distribution... " << endl;
	}

	std::ifstream vEddFile(velocitiesEddingtonFilename.c_str());
	std::string lineEdd;
	int nolines = 0;
	vector<double> normVector;
	while (std::getline(vEddFile, lineEdd)) {
		std::istringstream issEdd(lineEdd);
		double norm, velEdd;
		issEdd >> norm; // Read normalization constant
		normVector.push_back(norm); 
		vector<double> velEddVector;
		while(issEdd >> velEdd) { // Read columns and put in vector
			velEddVector.push_back(velEdd);
		}
		velEddMatrix.push_back(velEddVector); // Put vector of column in matrix
		nolines++;
	}
	//if (tid != 0) { cout << "NOT tid 0!!! ===== I've done the velocity stuff .. tid #" << tid  << endl;}
	
	if (tid == 0) {
	        // Uncomment this to verify reading in went well
		//cout << endl;
		//for (int i=0; i < 100; i++) { cout << scientific << velEddMatrix[i][0] << '\t'; }
		//cout << endl;	

		cout << "...done " << endl;
	} 
#endif

//cout << "Initializing random generator... " << endl;
	myrandom  = new randgen(inp,myRho0,radiusVector,velocityVector,velocityDistributions);
	//cout << "...done" << endl;

}


vector<double> PBHpopulation::runSimulation(Input* inp) {

	//ofstream outfile;
	//outfile.open("PBH_list.dat");

	//ofstream outfileRadio;
	//outfileRadio.open("../bigData/temp.dat");
	//outfileRadio << "R[kpc]" << "\t" << "v[km/s]" << "\t" << "Mdot[Msun/y]" << "\t" << "Lx[erg/s]" << "\t" << "Lradio[erg/s]" << "\t" << "RadioFlux[mJy]" << endl;
	//outfileRadio << "R" << "\t" << "v" << "\t" << "Mdot" << "\t" << "Lx" << "\t" << "Lradio" << "\t" << "RadioFlux" << endl;
	
	vector<double> outputVector;
	outputVector.push_back(0.);
	outputVector.push_back(0.);

	double cumulatedMass = 0.;
	double cumulatedMassInTheBulge = 0.;
	long int iPBH = 0;
	long int iPBHlarger = 0;
	long int iPBHsmaller = 0;
	int i = 0;

	double integratedFlux = 0.;
	double integratedRadioFlux = 0.;

	int numberRadioBrightSources = 0;
	int numberRadioBrightSourcesIonized = 0;
	int numberSKABrightSources = 0;
	int numberSKABrightSourcesIonized = 0;
	int numberSources = 0;
	int numberNuStarSources = 0;
	int numberNuStarSourcesIonized = 0;
	int numberChandraSources = 0.;
	int numberChandraSourcesIonized = 0.;

	ofstream sourceTable, sourceTableX;
	std::stringstream iss, iss2;
	iss << "sourceTableComplete_" <<  PBHmass_parameter << "_" << lambdaRicotti << "_LARGE.dat";
	string tableFileName = iss.str();
	iss2 << "XsourceTableComplete_" <<  PBHmass_parameter << "_" << lambdaRicotti << ".dat";
	string tableFileNameX = iss2.str();

	if (it == 0) {
		sourceTable.open(tableFileName.c_str());
		sourceTableX.open(tableFileNameX.c_str());
	}

	//ofstream bin1;
	//bin1.open("bin1.dat");

	if (it == 0) {
		sourceTable << "r[kpc]" << "\t" << "l" << "\t" << "b" << "\t " << "v[km/s]" << "\t" << "radioFlux[mJy]" << "\t" << endl;
		sourceTableX << "r[kpc]" << "\t" << "l" << "\t" << "b" << "\t " << "v[km/s]" << "\t" << "TotalXFlux[erg/s/cm2]" << "\t" << endl;
	}

#ifdef MAXWELLIANS
	gsl_spline *splineVelocity;
	gsl_interp_accel *acc;
  
	if (tid == 0) {
	  cout << "Initializing spline..." << endl;
	}
	int Nr = radiusVector.size();
	double r_vec_spline[Nr];
	double v_vec_spline[Nr];  
	for (int i=0; i < Nr; i++) {
		r_vec_spline[i] = radiusVector[i];
		v_vec_spline[i] = velocityVector[i];
	}
	//cout << "I'm here.. oh and btw: " << Nr << endl;
	splineVelocity = gsl_spline_alloc(gsl_interp_cspline, Nr);
	//cout << "Still alive!" << endl;
	gsl_spline_init(splineVelocity, r_vec_spline, v_vec_spline, Nr);
	acc = gsl_interp_accel_alloc ();
#endif

#ifdef EDDPHASESPACE
	// Creating velocity vector
	vector<double> velVector;
	for (int i=0; i < nvel; i++) {
		velVector.push_back(vminEdd + i*dvEdd);
	}
#endif

	if (tid == 0) {
	  cout << "Generating PBHs..." << endl;
	}

	// Initializing the seed for generating random values (EMD, radius, iAngle) + (velocity ifdef EDD) 
        unsigned short timedep = std::chrono::system_clock::now().time_since_epoch().count();
	unsigned short seed = timedep * (1 + omp_get_thread_num());

        // using mersene twister algorithm
        std::mt19937_64 generator(seed);
	// Using a uniform random dist from standard library
	std::uniform_real_distribution<double> unif;

#ifdef LOGNORMAL
	// Initializing the EMD distribution (lognormal)
        //std::default_random_engine generator (seed); 
	//std::mt19937_64 generator(seed); // using mersene twister algorithm
        std::lognormal_distribution<double> distribution (log(PBHmass_parameter),sigma_LN); // Parameters (m = log(mu), s = sigma)
	if (tid == 0) {
	cout << "The lognormal distribution is initialized with (mu;sigma) = (" << PBHmass_parameter << ";" << sigma_LN << ")" << endl;
	}
#endif

#ifdef POWERLAW
        double Mmin = PBHmass_parameter; // In units of Msun
	// Calculating the y_star value for the PL distribution
	double y_star = log(Mstar/Mmin)/(log(Mstar/Mmin)+2);
#endif

	for (int iRing=0; iRing < nr-1; iRing++) {
	//	cout << "I am looping over rings.. " << iRing << endl;

		//cout << "R = " << 0.5*(R_vec[iRing+1] + R_vec[iRing])/kpc << " kpc; vMean = " << vMean/(km/s) << " km/s" << endl;
		//cout << 100.*double(iRing)/double(nr) << "% " << endl;

		/*std::istringstream iss;
	  iss << "histFile_simulated_" << iRing << ".dat";
	  ofstream histFile;
	  histFile.open(iss.c_str());

	  std::istringstream iss2;
	  iss2 << "histFile_predicted_" << iRing << ".dat";
	  ofstream histFileMark;
	  histFileMark.open(iss2.c_str());*/

		double cumulatedMassInThisRing = 0.;
		int blackHolesInThisRing = 0;

		double predictedCumulatedMass, ringMass;

		predictedCumulatedMass = DMfraction*NFWMassUpToR(R_vec[iRing+1], myRho0);
		ringMass = DMfraction*NFWRingMass(R_vec[iRing], R_vec[iRing+1], myRho0);
		
#ifdef BURKERT
		predictedCumulatedMass = DMfraction*BurkertMassUpToR(R_vec[iRing+1], myRho0);
		ringMass = DMfraction*BurkertRingMass(R_vec[iRing], R_vec[iRing+1], myRho0);		
#endif		

#ifdef CONSERVATIVE
		predictedCumulatedMass = DMfraction*constantMassUpToR(R_vec[iRing+1]);
		ringMass = DMfraction*constantRingMass(R_vec[iRing], R_vec[iRing+1]);
#endif	  

		double deltaR = R_vec[iRing+1] - R_vec[iRing];
	        if (tid == 0) {
		  if (iRing<5) {
	    	    cout << endl << "*** r =  [" << R_vec[iRing]/kpc << "; " << R_vec[iRing+1]/kpc << "] ***"<< endl;
	    	    cout << "ring mass should be [Msun] = " << ringMass/mass_sun << endl;
	      	    cout << "mass up to this radius should be [Msun] = " << predictedCumulatedMass/mass_sun << endl;
	    	  }
		}

#ifdef EDDPHASESPACE
		// Create a spline for the function f(v) for each ring/shell
                gsl_spline *splineVelocityEdd;
                gsl_interp_accel *accEdd;

                //if (tid == 0) { cout << iRing+1 << ") Initializing velocity spline..." << endl; }
		vector<double> fvVector = velEddMatrix[iRing];
                int NvEdd = fvVector.size();
                double fv_vec_spline[NvEdd];
                double vel_vec_spline[NvEdd];
                for (int i=0; i < NvEdd; i++) {
			//if (tid==0) { if (iRing==0) { cout << velVector[i]/km*s << ' ' << fvVector[i] << endl; } }
                	fv_vec_spline[i] = fvVector[i];
                	vel_vec_spline[i] = velVector[i];
                }

                splineVelocityEdd = gsl_spline_alloc(gsl_interp_cspline, NvEdd);

		gsl_spline_init(splineVelocityEdd, vel_vec_spline, fv_vec_spline, NvEdd);
                accEdd = gsl_interp_accel_alloc ();
		
		/* // Tested if gsl_spline works correctly
		if (tid==0) { if (iRing==0) {
			for (int i=0; i < NvEdd; i++) { 
				cout << "vel,splvec " << velVector[i]/km*s << ' ' << vel_vec_spline[i]/km*s << " fv,splvec,spl "  << fvVector[i] << ' ' <<fv_vec_spline[i] << ' ' << gsl_spline_eval(splineVelocityEdd, velVector[i], accEdd) << endl;
			}
		} }*/

                // Doing it manually: 1) find box, 2) do hit or miss till true in while loop
                double yminE = 0;
                double ymaxE = 0;
                for (int i=0; i < nvel; i++) {
                        double fvtemp = fvVector[i];
                        if (fvtemp > ymaxE) { ymaxE = fvtemp; } // Finding the maximum of f(v)
                }
		if (tid == 0) { if (iRing < 5) { cout << iRing << ") ymin: " << yminE << " ymax: " << ymaxE << endl; } }
                // Box given by vminEdd < v < vmaxEdd on the x-axis and yminE < f(v) ymaxE on the y-axis
#endif

#ifdef SAVEVELIRING
		// Creating vector to store velocities of this ring in
	        vector<double> saveVeliRing;
#endif



		while (cumulatedMassInThisRing < ringMass)    {

			//extract the radius
			//double radius = myrandom->random_arbitrary_function("NFW", R_vec[iRing], R_vec[iRing+1]); //kpc
			double radius = R_vec[iRing] + unif(generator) * deltaR; 
			
			//if (iPBH%file_downsampling==0)
			//  outfile << radius/kpc << " ";
			blackHolesInThisRing++;
			iPBH++;
			
			double PBHmass;
#ifdef DELTA
			PBHmass = mass_sun * PBHmass_parameter; // Delta function
#endif

#ifdef LOGNORMAL
			PBHmass = mass_sun * distribution(generator); // Lognormal distribution
#endif

#ifdef POWERLAW
			// Broken power-law distribution, Deng & Vilenkin 2017:  http://adsabs.harvard.edu/abs/2017JCAP...12..044D
			// Using inverse CDF to draw from distribution
			double y_uniform = unif(generator); // Random number drawn from uniform dist (0,1)
			if (y_uniform < y_star) PBHmass = mass_sun * Mmin*pow((Mstar/Mmin),(y_uniform/y_star));
			else PBHmass = mass_sun * (4*Mstar/log(Mstar/Mmin)/log(Mstar/Mmin)) * y_star*y_star/(1-y_uniform)/(1-y_uniform);
#endif

			if (iPBH < 100010) PBHmasses[iPBH-1] = PBHmass; // Fill array with PBHmasses

			M_vec[iRing] += PBHmass;     
			if (radius < rBulge)
				cumulatedMassInTheBulge +=PBHmass;
			cumulatedMassInThisRing += PBHmass;
			cumulatedMass += PBHmass;
			//if ((iRing%NringDebug==0) && iPBH%Ndebug_bigcounter==0){
			//  cout << "\t PBH n. " << iPBH << endl;
			//  cout << "\t ring mass [Msun] =  " << cumulatedMassInThisRing/mass_sun << endl;			  
			//  cout << "\t cumulated mass [Msun] =  " << cumulatedMass/mass_sun << endl;
			//}

			// extract an angle on the sphere
			double iAngle_ = unif(generator);
			iAngle_ *= npix;
			int iAngle = (int) iAngle_;      
			double theta, phi;
			pix2ang_ring(nside, iAngle, &theta, &phi);
			double x = radius*sin(theta)*cos(phi);
			double y = radius*sin(theta)*sin(phi);
			double z = radius*cos(theta);
			double distance = sqrt( pow(x-xsun,2.) + pow(y-ysun,2.) + pow(z-zsun,2.)  ); 
			double l = atan2( (ysun-y), (xsun-x) ); // -pi... pi
			double b = asin(z/(xsun-x)); // -pi/2... pi/2, 0 = GP
			int il = lower_bound(l_vec.begin(), l_vec.end(), l) - l_vec.begin() - 1;
			int ib = lower_bound(b_vec.begin(), b_vec.end(), b) - b_vec.begin() - 1;

			//if ((iRing%NringDebug==0) && iPBH%Ndebug_bigcounter==0){
			//		cout << "\t radius [kpc] =  " << radius/kpc << " ir = " << iRing << endl;
			//			cout << "\t distance [kpc] =  " << distance/(kpc) << endl;
			//			cout << "\t l = " << l/deg <<	" b = " << b/deg << endl;
			//cout << il << " " << ib << endl;
			//}

			if (l < lmin || l > lmax)
				continue;
			if (b < bmin || b > bmax)
				continue;
			if (il < 0 || il > nl)
				continue;
			if (ib < 0 || ib > nb)
				continue;

			i++;
			if (tid == 0) {
			  if (i%Ndebug==0) cout << "\t\t"<< "*** PBH in the selected window n. " << i << endl;
			  if (i%Ndebug==0) cout << "\t\t"<<"radius [kpc] =  " << radius/kpc << " ir = " << iRing << endl;
			  if (i%Ndebug==0) cout << "\t\t"<<"distance [kpc] =  " << distance/(kpc) << endl;
			  if (i%Ndebug==0) cout << "\t\t"<<"l, b =  " << l/deg << "°; " << b/deg << "° " << endl;
			}
			//double vAverage = gsl_spline_eval(splineVelocity, radius, acc);
			//double sigmav = vAverage/ (2.*sqrt(2./M_PI));

			double BHvelocity = 200.*km/s;
			double vMean = 0.;
			double sigmaV = 0.;
			
			if (radius < 2.*kpc) {
#ifdef DISTRIBUTIONS		  
				BHvelocity = myrandom->random_arbitrary_function("generic", 0.*km/s, 380.*km/s, 0., radius);
				//if (velocity < soundSpeed)
					//velocity = soundSpeed;
#endif		  
#ifdef MAXWELLIANS
				vMean    = gsl_spline_eval(splineVelocity, 0.5*(R_vec[iRing+1] + R_vec[iRing]), acc);
				sigmaV   = vMean / sqrt(8./M_PI);
				BHvelocity = myrandom->random_arbitrary_function("MaxwellBoltzmann", 0.*km/s, 380.*km/s, sigmaV,  radius);
				//if (velocity < soundSpeed)
					//velocity = soundSpeed;
#endif	
#ifdef PESSIMISTIC
				BHvelocity = myrandom->random_arbitrary_function("MaxwellBoltzmann", 0.*km/s, 500.*km/s, 200.*km/s/sqrt(8./M_PI), radius);
#endif			
			}	
#ifdef EDDPHASESPACE
			// Why radius < 2 kpc???
			// Doing it manually: 1) find box, 2) do hit or miss till true in while loop
			// Box given by vminEdd < v < vmaxEdd on the x-axis and yminE < f(v) ymaxE on the y-axis

			bool need_rand_vel = true; // Using boolean to exit loop once value has been found
			while (need_rand_vel) { // Hit or miss method to generate random number
                		double xr = unif(generator);
                		double yr = unif(generator);
                		double xrand = vminEdd + (vmaxEdd - vminEdd) * xr;
                		double yrand = yminE + (ymaxE - yminE) * yr;
                		double f_of_x = gsl_spline_eval(splineVelocityEdd, xrand, accEdd); 

  		                if (yrand <= f_of_x) { BHvelocity = xrand; need_rand_vel = false; }
  			}
#endif
			
			double velocity = BHvelocity;

#ifdef SAVEVELIRING		
			// Storing velocities of this ring in iRing vector
			saveVeliRing.push_back(velocity);
#endif
			if (tid == 0) {
			  if (i%Ndebug==0) cout << "\t\t"<<"vMean [km/s] =  " << vMean/(km/s) << endl;
			  if (i%Ndebug==0) cout << "\t\t"<<"v [km/s] =  " << velocity/(km/s) << endl;		  
			}
			// Full Ferriere model
			double localGasDensity = gasDensity(x,y,z);
			// Gas approximation spheres
			//double localGasDensity = gasDensitySpheres(x,y,z);
			// Gas approximation ellipses
			//double localGasDensity = gasDensityEllipses(x,y,z);
			if (tid == 0) {
			  if (i%Ndebug==0) cout << "\t\t"<<"local gas density [cm^-3] = " << localGasDensity/(mass_proton/cm3) << endl;
			}
#ifndef NEWLAMBDA
			double MdotBondi = 4.*M_PI*localGasDensity*G_Newton*G_Newton*PBHmass*PBHmass / pow((velocity*velocity+soundSpeed*soundSpeed), 3./2.); //Bondi accretion rate
			if (tid == 0) {
			  if (i%Ndebug==0) cout << "\t\t"<<"MdotBondi [Msun/yr] =  " << MdotBondi/(mass_sun/year) << endl;
			}
						
			MdotBondi = 4.*M_PI*localGasDensity*G_Newton*G_Newton*PBHmass*PBHmass / pow((velocity*velocity+soundSpeed*soundSpeed), 3./2.); //Bondi accretion rate
			
			double Mdot = lambda * MdotBondi;
			//cout << "lambda = " << lambda << " | MdotBondi = " << MdotBondi << " | Mdot = " << Mdot << endl;
#endif
			

#ifdef NEWLAMBDA
			// =========== Implementing new lambda ============ Park & Ricotti 2013
			double Delta_rho;
			double mach = velocity/float(soundSpeed);
			double DD = (1 + mach*mach)*(1 + mach*mach) - 8*mach*mach*Delta_T;
			if (mach < machD) {
				Delta_rho = (1 + mach*mach + sqrt(DD))/(4*Delta_T);
			} else if (mach < machR) {
				Delta_rho = (1 + mach*mach)/(4*Delta_T);
			} else {
				Delta_rho = (1 + mach*mach - sqrt(DD))/(4*Delta_T);
			}
			double newlambda = Delta_rho * pow(2*Delta_T + mach*mach/(Delta_rho*Delta_rho),-3./2.);
			double MdotB = exp(1.5)*M_PI*localGasDensity*G_Newton*G_Newton*PBHmass*PBHmass / (soundSpeed*soundSpeed*soundSpeed);
			double Mdot = MdotB * newlambda;
			//cout << "new lambda = " << newlambda << " | MdotBondi = " << MdotB << " | Mdot = " << Mdot << endl; 
#endif

			if (tid == 0) {
			  if (i%Ndebug==0) cout << "\t\t"<<"Mdot [Msun/yr] =  " << Mdot/(mass_sun/year) << endl;
			}
			//double MdotCritical = 0.01 * EddingtonLuminosity(7.*mass_sun) / ( 0.1 * c_light * c_light); //Fender et al.
			
			double MdotCritical = 0.01 * EddingtonLuminosity(PBHmass) / ( 0.1 * c_light * c_light); //Fender et al.
			if (tid == 0) {
			  if (i%Ndebug==0) cout << "\t\t"<<"->MdotCritical  [Msun/yr] =  " << MdotCritical/(mass_sun/year) << endl; 
			}

			double mdot = Mdot / MdotEddington(PBHmass);
			if (tid == 0) {
			  if (i%Ndebug==0) cout << "\t\t"<<"->MdotEddington [Msun/yr] =  " << MdotEddington(PBHmass)/(mass_sun/year) << endl;
			  if (i%Ndebug==0) cout << "\t\t"<<"mdot [dimensionless] =  " << mdot << endl;
			}

			//double dimensionlessBolometricL = 0.;
			//if (mdot > 1) 
			//  dimensionlessBolometricL = min(0.1*mdot, 1.); // THIN DISK
			//else
			//  dimensionlessBolometricL = eta0*mdot*mdot;   // QUASI-SPHERICAL accretion
			//if (i%Ndebug==0) cout << "\t\t"<<"dimensionlessBolometricL =  " << dimensionlessBolometricL << endl;
			//double XLuminosity = 0.3 * dimensionlessBolometricL * EddingtonLuminosity(PBHmass);

			double eta =  eta0 * (Mdot/MdotCritical); //radiative efficiency according to Fender

			double XLuminosity = 0.3 * eta * Mdot * c_light * c_light;
						
			double rBondi = 2.*G_Newton * PBHmass * pow(velocity, -2.); //cm
			
			
			double photonFlux = 1.46e9 * (XLuminosity/(erg/s)); //s^-1
			double sigmaRecombination = 2.59e-13; //cm3/s

                        bool IONIZED = false; // Uncommented this bit to avoid errors (undefined)
#ifndef NEWLAMBDA
			// Implementing new lambda without Stromgren - Julien
			double rStromgren = pow( 3.* photonFlux / ( 4.*M_PI* pow( localGasDensity/(mass_proton) , 2.) * sigmaRecombination)  , 1./3);
			double rStromgrenMoving = sqrt( photonFlux / ( 4.*M_PI * ( localGasDensity/(mass_proton) ) * velocity  ) ); //cm
			
			double t_recombination = 1./(sigmaRecombination * localGasDensity/mass_proton); //http://w.astro.berkeley.edu/~echiang/rad/ps6ans.pdf
			
			double t_ionization = localGasDensity/mass_proton * 4/3. * M_PI * pow(rStromgrenMoving, 3.) / photonFlux;
			double t_crossing   = rBondi/velocity;
			
			if (tid == 0) {
			  if (i%Ndebug==0) cout << "ionization time "     << t_ionization/year << " y " << endl;
			  if (i%Ndebug==0) cout << "Bondi crossing time " << t_crossing/year << " y " << endl;
			}
						
			if (t_ionization < t_crossing) {
				IONIZED = true;
				if (tid == 0) {
				  if (i%Ndebug==0) cout << "Ionization time smaller than Bondi time! " << endl;
				}
				MdotBondi = 4.*M_PI*localGasDensity*G_Newton*G_Newton*PBHmass*PBHmass / pow((velocity*velocity+soundSpeedIonized*soundSpeedIonized), 3./2.); //Bondi accretion rate
				Mdot = lambda * MdotBondi;
				mdot = Mdot / MdotEddington(PBHmass);
				eta =  eta0 * (Mdot/MdotCritical);
				XLuminosity = 0.3 * eta * Mdot * c_light * c_light;
				////if (i%Ndebug==0) cout << "\t\t"<<"MdotBondi for IONIZED gas [Msun/yr] =  " << MdotBondi/(mass_sun/year) << endl;
			}
#endif
						
			double sigmaPhoto = 6.3e-18; //cm2 //sigma = sigmaPhoto * (E / 13,6eV)^(-3)
			
			double Eion = 13.6*eV;
			//double N = 0.7 * referenceLuminosity / pow(Eion, -0.7); //S = \int{ N E^-1.7 dE }
			//double t_ionization = 1./  ( sigmaPhoto*N/(4.*M_PI*pow(rStromgren,2.)) * (1./3.7) * pow(Eion,3.) * pow(Eion,-3.7) );
			
			double EX1 =  2.*keV;
			double EX2 = 10.*keV;
			//double E_nustar = 80.*keV;
			double A = (1. - Xslope)*XLuminosity/keV / ( pow(EX2/keV, 1.-Xslope)  -  pow(EX1/keV, 1.-Xslope)  );
			//double XLuminosity_Nustar = A * pow(E_nustar/keV, -Xslope);
			double EX3 = 10*keV; //10*keV;
			double EX4 = 40.*keV; //40*keV;
#ifdef CHANDRA			
			EX3 = 0.5*keV; //10*keV;
			EX4 = 8.*keV; //40*keV;			
#endif
			double XLuminosity_Nustar = A*keV/(1. - Xslope) * ( pow(EX4/keV, 1.-Xslope)  -  pow(EX3/keV, 1.-Xslope)  );
			//double XFlux_Nustar =  XLuminosity_Nustar / (4.*M_PI*distance*distance);
			if (tid == 0) {
		 	  if (i%Ndebug==0) cout << "\t\t"<<"L_X_2_10keV [erg/s] =  " << XLuminosity << "; L_X_10_40keV [erg/s]  " << XLuminosity_Nustar << endl;
			}
			//if (i%Ndebug==0) cout << "\t\t"<<"eta0*lambda^2/Mcritical =  " << XLuminosity/(MdotBondi*MdotBondi*c_light*c_light) << endl;
			//if (i%Ndebug==0) cout << "\t\t"<<"eta0*lambda^2/Mcritical [SgrA*] =  " << Lx_SgrA/(MdotBondi_SgrA*MdotBondi_SgrA*c_light*c_light) << endl;

			double totalXFlux = XLuminosity_Nustar / (4.*M_PI*distance*distance); // erg/s/cm^2
			if (tid == 0) {
			  if (i%Ndebug==0) cout << "\t\t"<<"X-ray flux [erg/s/cm^2] =  " << totalXFlux << endl;
			}

			//double XphotonFlux = totalXFlux / photonEnergy;
			double XphotonFlux = A/(-Xslope) * ( pow(EX4/keV, -Xslope)  -  pow(EX3/keV, -Xslope)  ) / (4.*M_PI*distance*distance);
			if (tid == 0) {
			  if (i%Ndebug==0) cout << "\t\t"<<"X-ray photon flux [photons/s/cm^2] =  " << XphotonFlux << endl;
			}
			// according to http://arxiv.org/pdf/1105.3211v2.pdf
			// log(Lx) = 1.45 log(Lr) - 0.88log(Mbh) - 6.07

			double logLx  = log10(XLuminosity);
			double logMbh = log10(PBHmass/mass_sun);
			double radioLuminosity = pow(10., (logLx + 0.88*logMbh + 6.07)/1.45 );

			if (tid == 0) {
			  if (i%Ndebug==0) cout << "\t\t" << "Lradio [erg/s] =  " << radioLuminosity << endl;
			  if (i%Ndebug==0) cout << "\t\t" << "log(Lradio) log[erg/s] =  " << log10(radioLuminosity) << endl;
			  if (i%Ndebug==0) cout << "\t\t" << "log(Lx) - xi log(Mbh) log[erg/s] =  " << logLx + 0.88*logMbh  << endl;
			}

			double radioFlux = radioLuminosity / (4.*M_PI*distance*distance); //erg/cm2/s
			double radioFlux_mJy = 1000. * radioFlux / obsFrequency / Jy_erg_conversion;

			if (tid == 0) {
			  if (i%Ndebug==0) cout << "\t\t" << "FluxRadio [mJy] =  " << 1000. * radioFlux / obsFrequency / Jy_erg_conversion << endl;
			}

			if (iPBH < 100010){
			  radioFluxes[iPBH-1] = radioFlux_mJy; // Fill array with PBHmasses
			  xrayFluxes[iPBH-1] = XphotonFlux;
			}

			//if (XLuminosity_Nustar >= Chandra_threshold){
			if (XphotonFlux >= XPhotonFluxChandraThreshold){
				numberChandraSources += 1;
				vDetXray.push_back(velocity); // TURN BACK ON!
				if (IONIZED == true)
					numberChandraSourcesIonized += 1;
				//cout << "\t\t\t" << "WOW I found a bright X-ray source at l = " << l_vec[il]/deg << "°; v [km/s] = " << velocity/(km/s) << endl;
				//if (t_ionization < t_crossing) cout << "\t\t\t t_Ionizaton < t_Crossing! Gas may be ionized " << endl;
				//if (rStromgrenMoving < rBondi) cout << "\t\t\t very small Strongren sphere! Gas ionization is negligible " << endl;
				//if ( (t_ionization < t_crossing) && (rStromgrenMoving > rBondi)  ) 
				//s	cout << "\t\t\t Conflicting conditions " << endl;
			}

			//if (R_vec[iRing] <= 200.*pc)
			//  continue;

			//if ((radius > 0.000*kpc) && (radius < 0.050*kpc))
			//  bin1 << velocity/(km/s) << endl;
			//cout << radioFlux_mJy << " >=??? " << RadioThreshold1GHz_mJy << endl;
			if (radioFlux_mJy >= RadioThreshold1GHz_mJy) {			
				numberRadioBrightSources += 1;
				vDetRadio.push_back(velocity); // TURN BACK ON!
				//cout << "Radio source found!" << endl;
				if (IONIZED == true)
					numberRadioBrightSourcesIonized += 1;
			}

			if (radioFlux_mJy >= RadioThresholdSKA_mJy) {
				numberSKABrightSources += 1;
				if (IONIZED == true)
					numberSKABrightSourcesIonized += 1;
				sourceTable << radius/kpc << "\t" << l << "\t" << b << "\t " << velocity/(km/s) << "\t" << radioFlux_mJy << "\t" << endl;
			}

#ifdef NUSTAR			
			if (XLuminosity_Nustar >= Nustar_threshold) {
				if (tid == 0) {					
				  cout << endl << "bright X-ray source found! " << endl;
				  cout << "r = " << radius/kpc << " kpc" << endl;
				  cout << "n = " << localGasDensity/mass_proton << " cm^(-3)" << endl;
				  cout << "v = " << velocity/(km/s) << " km/s" << endl;
				  cout << "Mdot = " << Mdot << endl;
				  cout << "MdotCritical = " << MdotCritical << endl;
				  cout << "ionization time "     << t_ionization/year << " y " << endl;
				  cout << "Bondi crossing time " << t_crossing/year << " y " << endl;
				  cout << "Xflux Chandra = " << XphotonFlux << endl;
				}
				numberNuStarSources += 1;
				if (IONIZED == true)
					numberNuStarSourcesIonized += 1;
				//if ((it == 0) && (PBHmass > 10.*mass_sun))
					//sourceTableX << radius/kpc << "\t" << l << "\t" << b << "\t " << velocity/(km/s) << "\t" << totalXFlux << "\t" << endl;
				
			}
#endif			

			//outfileRadio << radius/kpc << "\t" << velocity/(km/s) << "\t" << Mdot/(mass_sun/year) << "\t" << XLuminosity << "\t" << radioLuminosity << "\t" << 1000. * radioFlux / bandWidth / Jy_erg_conversion << endl;

			PBHmap[index(il,ib)] += totalXFlux;
			PBHRadioMap[index(il,ib)] += radioFlux;
			integratedFlux += totalXFlux;
			integratedRadioFlux += radioFlux;
		}

		#ifdef EDDPHASESPACE                  
        	gsl_spline_free (splineVelocityEdd);
        	gsl_interp_accel_free (accEdd);
        	//if (tid == 0) {
          	//	cout << "Spline deleted" << endl;
        	//}
		#endif



		/*if (iRing%NringDebug==0){

		  cout << endl << "PBHs in this ring: " << blackHolesInThisRing << endl;
		  cout << "actual mass in this ring [Msun]: " << blackHolesInThisRing*PBHmass/mass_sun << endl << endl;
	  }*/

		//histFile.close();
		//histFileMark.close();

#ifdef SAVEVELIRING 
		saveVVector.push_back(saveVeliRing);
#endif

	}
	if (it == 0) {
		sourceTable.close();
		sourceTableX.close();
	}
	//bin1.close();

	//cout << endl;
	if (tid == 0) {
	  cout << "...done!" << endl;
	  cout << "total PBHs within Rmax = " << iPBH << endl;
	  cout << "total PBHs within Rmax in the selected window = " << i << endl;
	  cout << "bulge mass from the sum of all PBHs [Msun] = " << cumulatedMassInTheBulge/mass_sun << endl;

#ifndef BURKERT	
	  cout << "bulge mass from NFW integration times DMfraction [Msun] " << DMfraction*NFWMassUpToR(rBulge,myRho0 )/mass_sun << endl;
#endif
#ifdef BURKERT	
	  cout << "bulge mass from NFW integration times DMfraction [Msun] " << DMfraction*BurkertMassUpToR(rBulge,myRho0 )/mass_sun << endl;
#endif
	}
	/*
  cout << "predicted integrated X-ray flux from PBHs from the selected window [erg/s/cm^2]: " << integratedFlux << endl;
  cout << "predicted integrated radio flux from PBHs from the selected window [erg/s/cm^2]: " << integratedRadioFlux << endl;
  double area = 4. * angRadius * angRadius; //srad
  cout << "area of the selected window [square deg] " << 4. * angRadius/deg * angRadius/deg << endl;
  integratedFlux /= area;
  integratedRadioFlux /= area;
  outfile << endl;
  //outfileRadio << endl;
  outfile.close();
  outfileRadio.close();
  cout << "predicted X-ray flux density from PBHs from the selected window [erg/s/cm^2/sr]: " << integratedFlux << endl;
  cout << "predicted radio flux density from PBHs from the selected window [erg/s/cm^2/sr]: " << integratedRadioFlux << endl;
  cout << "measured X-ray flux density from the selected window [erg/s/cm^2/sr] TO BE CHECKED: " << measuredFlux << endl;

  ofstream outfile_binned;
  outfile_binned.open("PBH_list_binned.dat");
  for (int ir = 0; ir < nr; ir++) {
	  outfile_binned << R_vec[ir]/kpc << "\t" <<  M_vec[ir]/mass_sun << endl;
	  if (ir % 2000 == 0)
		  cout << R_vec[ir]/kpc << "\t" <<  M_vec[ir]/mass_sun << endl;
  }
  outfile_binned.close();

  cout << "*** Number of bright sources (F > 10^-6/cm2/s) " << endl;
  double number = 0;
  ofstream outfile_bright;
  outfile_bright.open("bright_sources_binned.dat");
  cout << "bmax-bmin = " << (bmax-bmin)/arcmin << " arcmin " << endl;
  cout << "deltal = " << deltal/arcmin << " arcmin " << endl;
  double area_arcmin2 = (bmax-bmin)*deltal/arcmin/arcmin;
  for (int il=0; il<l_vec.size(); il++){
	  //cout << l_vec[il]/deg << "°\t" <<  brightSourceNumber.at(il) << endl;
	  number += brightSourceNumber.at(il);
	  outfile_bright << l_vec[il]/deg << "\t" << brightSourceNumber.at(il)/area_arcmin2 << "\n" ;
	  cout << l_vec[il]/deg << "°\t" <<  brightSourceNumber.at(il)/area_arcmin2 << "/arcmin2 " << endl;
  }
  outfile_bright.close();
  cout << "TOTAL NUMBER OF SOURCES: " << number << endl;*/

	//cout << RAND_MAX << endl; 
	if (tid == 0) {
	  cout << "****** Number of bright VLA sources (F > 2.47 mJy): " << numberRadioBrightSources << " error: " << sqrt(numberRadioBrightSources) << endl;
	  cout << "****** Number of bright SKA sources (F > 1 muJy): " << numberSKABrightSources << endl;
	  cout << "****** Number of bright NuStar sources: " << numberNuStarSources << " error: " << sqrt(numberNuStarSources) << endl;
	  cout << "****** Number of sources: " << numberSources << endl;
	  cout << "****** Number of PBHs with velocity larger and smaller resp.: " << iPBHlarger << " (" << 100*iPBHlarger/float(iPBH) << "%); " << iPBHsmaller << " (" << 100*iPBHsmaller/float(iPBH) << "%)" << endl;
	  cout << "****** Total number of PBHs generated: " << iPBH << " " << iPBHlarger+iPBHsmaller << endl;
	}

#ifdef MAXWELLIANS		    
	gsl_spline_free (splineVelocity);
	gsl_interp_accel_free (acc);
	if (tid == 0) {
	  cout << "Spline deleted" << endl;
	}
#endif  

	delete myrandom;
	if (tid == 0) {
	  cout << "Randgen object deleted" << endl;
	}

	// I am using ifndef NUSTAR and ifdef CHANDRA
#ifdef NUSTAR
	outputVector[0] = numberNuStarSources;
	outputVector[1] = numberNuStarSourcesIonized;
#endif	

	// VLA radio sources
#ifndef NUSTAR
	outputVector[0] = numberRadioBrightSources;
//	outputVector[1] = numberRadioBrightSourcesIonized;
#endif
	
#ifdef SKA
	outputVector[0] = numberSKABrightSources;
	outputVector[1] = numberSKABrightSourcesIonized;
#endif
	
#ifdef CHANDRA
	outputVector[1] = numberChandraSources;
//	outputVector[1] = numberChandraSourcesIonized;
#endif	


	
	return outputVector;
}

void PBHpopulation::WriteMassToFile() {

        stringstream tempString;
#ifdef DELTA
	tempString << "PBHmassesD_" << PBHmass_parameter << ".dat";
#endif
#ifdef LOGNORMAL
        tempString << "PBHmassesLN_" << PBHmass_parameter << "_" << sigma_LN << ".dat";
#endif
#ifdef POWERLAW
	tempString << "PBHmassesPL_" << PBHmass_parameter << "_" << Mstar << ".dat"; //Parameters Mmin=PBHmass_parameter _ Mstar
#endif
        string outFileName = tempString.str();
        cout << "Writing masses to " << outFileName << endl;

        ofstream outfile;
	outfile.open(outFileName.c_str());

        for (int i=0; i<100000; i++)
        	outfile << PBHmasses[i]/mass_sun << "\t" << radioFluxes[i] << "\t" << xrayFluxes[i] << endl;

        outfile.close();

}


void PBHpopulation::WriteToFile() {

	ofstream outfile;
	outfile.open("PBH_map.dat");

	outfile << "l" << "\t" << "b" << "\t" << "flux" << endl;
	for (int il=0; il<nl; il++)
		for (int ib=0; ib<nb; ib++)
			outfile << l_vec[il]/deg << "\t" << b_vec[ib]/deg << "\t" << PBHmap[index(il,ib)] << endl;

	outfile.close();    
}

void PBHpopulation::RadioWriteToFile() {

	ofstream outfile;
	outfile.open("PBH_radio_map.dat");

	outfile << "l" << "\t" << "b" << "\t" << "flux" << endl;
	for (int il=0; il<nl; il++)
		for (int ib=0; ib<nb; ib++)
			outfile << l_vec[il]/deg << "\t" << b_vec[ib]/deg << "\t" << PBHRadioMap[index(il,ib)] << endl;

	outfile.close();    
}


// --------------------------------------------------------------------------------------------------------------------------------------------
// destructor ---------------------------------------------------------------------------------------------------------------------------------

PBHpopulation::~PBHpopulation()
{
	if (tid == 0) {
	  cout << "I am the destructor " << endl << endl;
	}
}

void PBHpopulation::pix2ang_ring(long nside, long ipix, double *theta, double *phi) {
	/*
     c=======================================================================
     c     gives theta and phi corresponding to pixel ipix (RING) 
     c     for a parameter nside
     c=======================================================================
	 */

	int nl2, nl4, npix, ncap, iring, iphi, ip, ipix1;
	double  fact1, fact2, fodd, hip, fihip;
	double PI=M_PI;
	//      PARAMETER (pi     = 3.1415926535897932384626434d0)
	//      parameter (ns_max = 8192) ! 2^13 : largest nside available

	int ns_max=8192;

	if( nside<1 || nside>ns_max ) {
		fprintf(stderr, "%s (%d): nside out of range: %ld\n", __FILE__, __LINE__, nside);
		exit(0);
	}
	npix = 12*nside*nside;      // ! total number of points
	if( ipix<0 || ipix>npix ) { // I removed -1 in "ipix>npix-1", because RAND/RAND_MAX includes both 0 and 1.
		fprintf(stderr, "%s (%d): ipix out of range: %ld\n", __FILE__, __LINE__, ipix);
		exit(0);
	}

	ipix1 = ipix + 1; // in {1, npix}
	nl2 = 2*nside;
	nl4 = 4*nside;
	ncap = 2*nside*(nside-1);// ! points in each polar cap, =0 for nside =1
	fact1 = 1.5*nside;
	fact2 = 3.0*nside*nside;

	if( ipix1 <= ncap ) {  //! North Polar cap -------------

		hip   = ipix1/2.;
		fihip = floor(hip);
		iring = (int)floor( sqrt( hip - sqrt(fihip) ) ) + 1;// ! counted from North pole
		iphi  = ipix1 - 2*iring*(iring - 1);

		*theta = acos( 1. - iring*iring / fact2 );
		*phi   = (1.*iphi - 0.5) * PI/(2.*iring);
	}
	else if( ipix1 <= nl2*(5*nside+1) ) {//then ! Equatorial region ------

		ip    = ipix1 - ncap - 1;
		iring = (int)floor( ip / nl4 ) + nside;// ! counted from North pole
		iphi  = (int)fmod(double(ip),double(nl4)) + 1;

		fodd  = 0.5 * (1 + fmod((double)(iring+nside),2));//  ! 1 if iring+nside is odd, 1/2 otherwise
		*theta = acos( (nl2 - iring) / fact1 );
		*phi   = (1.*iphi - fodd) * PI /(2.*nside);
	}
	else {//! South Polar cap -----------------------------------

		ip    = npix - ipix1 + 1;
		hip   = ip/2.;
		/* bug corrige floor instead of 1.* */
		fihip = floor(hip);
		iring = (int)floor( sqrt( hip - sqrt(fihip) ) ) + 1;//     ! counted from South pole
		iphi  = (int)(4.*iring + 1 - (ip - 2.*iring*(iring-1)));

		*theta = acos( -1. + iring*iring / fact2 );
		*phi   = (1.*iphi - 0.5) * PI/(2.*iring);
	}
}

void PBHpopulation::ang2pix_ring( const long nside, double theta, double phi, long *ipix) {
	/*
    c=======================================================================
    c     gives the pixel number ipix (RING) 
    c     corresponding to angles theta and phi
    c=======================================================================
	 */

	int nl2, nl4, ncap, npix, jp, jm, ipix1;
	double  z, za, tt, tp, tmp;
	int ir, ip, kshift;

	double piover2 = 0.5*M_PI;
	double twopi=2.0*M_PI;
	double z0=2.0/3.0; 

	z = cos(theta);
	za = fabs(z);
	if( phi >= twopi)  phi = phi - twopi;
	if (phi < 0.)     phi = phi + twopi;
	tt = phi / piover2;//  ! in [0,4)

	nl2 = 2*nside;
	nl4 = 4*nside;
	ncap  = nl2*(nside-1);// ! number of pixels in the north polar cap
	npix  = 12*nside*nside;

	if( za <= z0 ) {

		jp = (int)floor(nside*(0.5 + tt - z*0.75)); /*index of ascending edge line*/
		jm = (int)floor(nside*(0.5 + tt + z*0.75)); /*index of descending edge line*/

		ir = nside + 1 + jp - jm;// ! in {1,2n+1} (ring number counted from z=2/3)
		kshift = 0;
		if (fmod(ir,2)==0.) kshift = 1;// ! kshift=1 if ir even, 0 otherwise

		ip = (int)floor( ( jp+jm - nside + kshift + 1 ) / 2 ) + 1;// ! in {1,4n}
		if( ip>nl4 ) ip = ip - nl4;

		ipix1 = ncap + nl4*(ir-1) + ip ;
	}
	else {

		tp = tt - floor(tt);//      !MOD(tt,1.d0)
		tmp = sqrt( 3.*(1. - za) );

		jp = (int)floor( nside * tp * tmp );// ! increasing edge line index
		jm = (int)floor( nside * (1. - tp) * tmp );// ! decreasing edge line index

		ir = jp + jm + 1;//        ! ring number counted from the closest pole
		ip = (int)floor( tt * ir ) + 1;// ! in {1,4*ir}
		if( ip>4*ir ) ip = ip - 4*ir;

		ipix1 = 2*ir*(ir-1) + ip;
		if( z<=0. ) {
			ipix1 = npix - 2*ir*(ir+1) + ip;
		}
	}
	*ipix = ipix1 - 1;// ! in {0, npix-1}

}

