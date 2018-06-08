#include "constants.h"
#include "randgen.h"
#include "input.h"
#include <string>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h> 
#include <gsl/gsl_vector.h>

#include <omp.h> // For parallelization (to retrieve the thread number)
#include <chrono> // For seed random variable generator
#include <random> // For random generator

using namespace std;

randgen::randgen() { 

}

randgen::randgen(Input* inp_, double rho0_, vector<double> radiusVector_, vector<double> velocityVector_, map<double, vector<double> > velocityDistributions_) {
 inp = inp_;
 rho0 = rho0_;
 radiusVector          = radiusVector_;
 velocityVector        = velocityVector_;
 velocityDistributions = velocityDistributions_;

 // Initializing the seed for the two random variables
 unsigned short timedep = std::chrono::system_clock::now().time_since_epoch().count();
 tid = omp_get_thread_num();
 seed1 = timedep * (1 + tid);
 // using mersene twister algorithm
 generator.seed(seed1);
}

randgen::~randgen() {
	if (tid == 0) {	
	  cout << "I am the destructor of randgen class" << endl;
	}
}

// *****************************************

double randgen::NFW(double radius){

	double radius_ = radius;
	if (radius < 10.*pc)
		radius_ = 10.*pc;
	
    double denominator = (radius_/Rs)*(1+radius_/Rs)*(1+radius_/Rs);
	return rho0/(denominator);
}

double randgen::constantDistribution(double radius){

	return 1.;
}

/*double randgen::luminosityFunction(double L){

    double result = pow(L, -inp->slopeLum);
    return result;
}*/

double randgen::MaxwellBoltzmann(double v, double sigmav){

    double result = sqrt(2./M_PI)*v*v/(sigmav*sigmav*sigmav)*exp(-(v*v)/(2.*sigmav*sigmav));
    
    return result;
}

//double randgen::logNormal(double M, double mu, double sigma){
//	
//	double dPhidM = exp(-pow(log(M/mu),2)) / (sqrt(2*M_PI)*sigma*M)	
//
//	return dPhidM;
//}

double randgen::genericVDistribution(double radius, double v){

	// radialBins contains the bin edges
	// 0  1 2 3 4 5
	// radiusVector contains the bin centers
	//   0 1 2 3 4
	
    double result = 0.;
    
    int nR = radiusVector.size();
    int nV = velocityVector.size();  
    int whichBin = 0;
    
    vector<double> radialBins;
    radialBins.push_back(0.);
    for (int ir = 0; ir < nR-1; ir++)
    	radialBins.push_back(0.5*(radiusVector[ir]+radiusVector[ir+1]));
    radialBins.push_back( radiusVector[nR-1] + 0.5*(radiusVector[nR-1] - radiusVector[nR-2]) );
    
    /*cout << "Radial bins" << endl;
    for (int ir = 0; ir < nR+1; ir++)
    	cout << radialBins[ir]/kpc << "\t";
    cout << endl;*/
    
    for (int ir=0; ir<nR+1; ir++)     	
    	if ((radius >= radialBins[ir]) && (radius < radialBins[ir+1]))
    		whichBin = ir;
    
    //cout << "radius = " << radius/kpc << endl;
    //cout << "radial bin n. " << whichBin << endl; 
    
    double  v_vec_spline[nV+1];
    double fv_vec_spline[nV+1];
    v_vec_spline[0] = 0.;
    fv_vec_spline[0] = 1.e-7;
    for (int i=0; i < nV; i++) {
  	  v_vec_spline[i+1]  = velocityVector[i];
  	  fv_vec_spline[i+1] = velocityDistributions[whichBin][i];
  	  //cout << v_vec_spline[i]/(km/s) << "\t" << fv_vec_spline[i] << endl; 
    }
    gsl_spline *mySplineVelocity;
    gsl_interp_accel *myAcc;
    myAcc = gsl_interp_accel_alloc ();
    mySplineVelocity = gsl_spline_alloc(gsl_interp_cspline, nV+1);
    gsl_spline_init(mySplineVelocity, v_vec_spline, fv_vec_spline, nV+1);
    //cout << "velocity = " << v/(km/s) << endl;
    //cout << "f(v) = " << gsl_spline_eval(splineVelocity, v, acc) << endl;
    result = gsl_spline_eval(mySplineVelocity, v, myAcc);
    gsl_spline_free (mySplineVelocity);
    gsl_interp_accel_free (myAcc);
    
    return result;
    
    
/*    int iDown = 0;
    int iUp = 0;
    int ivDown = 0;
    int ivUp = 0;
    double deltaR = radiusVector[0];
    double deltaV = velocityVector[0];
    double downR  = 1.;
    double upR = 0.;
    double downV = 1.;
    double upV = 0.;
    
	for (int i = 0; i < radiusVector.size()-1; i++) {
    	if ( (radius >= radiusVector[i]) && (radius < radiusVector[i+1]) ){
    		iDown = i;
    		iUp   = i+1;
    		deltaR = radiusVector[i+1] - radiusVector[i];
    		downR = (radiusVector[iUp] - radius)/deltaR;
    		upR   = (radius - radiusVector[iDown])/deltaR;
    	}
    }
	if (radius >= radiusVector.back()) {
		iDown = nR-1;
		iUp   = nR-1;
		downR = 1.;
		upR   = 0.;
		//cout << "warning" << endl;
	}
	
	for (int iv = 0; iv < velocityVector.size()-1; iv++) {
    	if ( (v >= velocityVector[iv]) && (v < velocityVector[iv+1]) ){
    		ivDown = iv;
    		ivUp   = iv+1;
    		deltaV = velocityVector[iv+1] - velocityVector[iv];
    		downV = (velocityVector[ivUp] - v)/deltaV;
    		upV   = (v - velocityVector[ivDown])/deltaV;
    	}
    }
	if (v >= velocityVector.back()) {
		ivDown = nV-1;
		ivUp   = nV-1;
		downV = 1.;
		upV   = 0.;
		//cout << "warning" << endl;

	}
	
	result = velocityDistributions[iDown][ivDown]*downR*downV + velocityDistributions[iDown][ivUp]*downR*upV 
			+  velocityDistributions[iUp][ivDown]*upR*downV +   velocityDistributions[iUp][ivUp]*upR*upV;
    
    return result; */
}


// *****************************************

double randgen::random_arbitrary_function(string func_name, double xmin, double xmax, double xmean, double R) {

	// randgen_numsteps is in constants.h

	double ymin = 0.;
	double ymax;

	// Initialize ymin and ymax to create box 
	if (func_name == "NFW")
		ymin = NFW(xmin);
	else if (func_name == "constant")
		ymin = constantDistribution(xmin);
	else if (func_name == "MaxwellBoltzmann")
		ymin = MaxwellBoltzmann(xmin,xmean);
	else if (func_name == "generic")
		ymin = genericVDistribution(R,xmin);

	ymax = ymin;

	// Find lowest and highest y-values in xrange
	for (int i=0; i<randgen_numsteps; i++) {
    		double x = xmin + (xmax - xmin) * i / randgen_numsteps;    
		double y = 0.	;	
	
		if (func_name == "NFW")
			y = NFW(x);
		else if (func_name == "constant")
			y = constantDistribution(x);		
        	else if (func_name == "MaxwellBoltzmann")
			y = MaxwellBoltzmann(x,xmean);
    		else if (func_name == "generic")
    			y = genericVDistribution(R,x);
        
		if (y < ymin) 
			ymin = y;
	        if (y > ymax) 
			ymax = y;
	}
	// Box is created


	while (true) {
		double xr = unif(generator);
	        double yr = unif(generator);
	        double xrand = xmin + (xmax - xmin) * xr;
        	double yrand = ymin + (ymax - ymin) * yr;
		double f_of_x = 0.;
		
		if (func_name == "NFW")
			f_of_x = NFW(xrand);
		else if (func_name == "constant")
			f_of_x = constantDistribution(xrand);		
        	else if (func_name == "MaxwellBoltzmann")
			f_of_x = MaxwellBoltzmann(xrand,xmean);
	        else if (func_name == "generic")
			f_of_x = genericVDistribution(R,xrand);
        
        	if (yrand <= f_of_x)
            		return xrand;            		
	}
}
