#include <iostream>
#include <stdio.h>
#include <cstdlib> 
#include <ctime>

#include "constants.h"
#include "errorcode.h"
#include "PBHpopulation.h"
#include "input.h"
#include "randgen.h"

using namespace std;

int main(int argc, char** argv) {
    
  cout << "PBH template builder" << endl;
    
  if (argc != 2) {
        cerr << "Usage: ./PBH <xml-file>" << endl;
        exit(-1);
  }

  time_t time_s;
  time_t time_e;

  Input* inp = new Input();
  inp->LoadFile( argv[1] );

  cout << "Parameter file read successfully" << endl;
  inp->Print();
    
  srand(time(NULL));    

  cout << "Constructing the template... " << endl;
  
  ofstream outfile_constraints;
  outfile_constraints.open(mainFile.c_str());
  
  vector<double> lambdaVec;
  vector<double> etaVec;
  double logLambdaFactor = log10(lambdaMax/lambdaMin) / (numLambda-1);
  double logEtaFactor    = log10(etaMax/etaMin) / (numEta-1);
  
  cout << "lambdaVec" << endl;
  for (int i=0; i<numLambda; i++) {
	  lambdaVec.push_back( pow(10., log10(lambdaMin) + i*logLambdaFactor  ) );
	  cout << lambdaVec.back() << "\t";
  }
  cout << endl;
  cout << "etaVec" << endl;
  for (int i=0; i<numEta; i++) {
	  etaVec.push_back( pow(10., log10(etaMin) + i*logEtaFactor  ) );
	  cout << etaVec.back() << "\t";
  }
  cout << endl;
  
  /*for (int iLambda = 0; iLambda < numLambda; iLambda++)
	  for (int iEta = 0; iEta < numEta; iEta++) {
 
		  cout << endl;
		  time(&time_s); // fix initial time
		  PBHpopulation* PBH = new PBHpopulation(inp, lambdaVec[iLambda], etaVec[iEta]);
		  double value = PBH->runSimulation(inp);
		  time(&time_e); // fix final time
		  cout << "Simulation computed in " << (double)(time_e-time_s) << " s." << endl;
		  outfile_constraints << lambdaVec[iLambda] << "\t" << etaVec[iEta] << "\t" << value << endl;
		  cout << endl;
		  cout <<  "****** " << endl;
		  cout <<  "****** " << lambdaVec[iLambda] << "\t" << etaVec[iEta] << "\t" << value << endl;
		  cout <<  "****** " << endl;
		  cout << endl;
		  delete PBH;

	  }*/
  
  
  /*cout << endl;
  time(&time_s); // fix initial time
  PBHpopulation* PBH = new PBHpopulation(inp, 1., 0.01);
  double value = PBH->runSimulation(inp);
  time(&time_e); // fix final time
  cout << "Simulation computed in " << (double)(time_e-time_s) << " s." << endl;
  outfile_constraints << 1. << "\t" << 0.01 << "\t" << value << endl;
  cout << endl;
  cout <<  "****** " << endl;
  cout <<  "****** " << 1. << "\t" << 0.01 << "\t" << value << endl;
  cout <<  "****** " << endl;
  cout << endl;
  delete PBH;*/
  
  
  
  outfile_constraints.close();
  
  //cout << "Saving the template" << endl;    
  //PBH->WriteToFile();
  //PBH->RadioWriteToFile();

  delete inp;

  return 0;

}
