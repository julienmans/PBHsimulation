#include "input.h"
#include "constants.h"
#include "tinyxml.h"
#include "errorcode.h"
#include "math.h"

#include <string>
#include <stdlib.h>
#include <algorithm>
#include <sstream>
#include <vector>

using namespace std;

int Input::QueryIntAttribute(string obj, TiXmlElement* el) {
  TiXmlElement* el1 = el->FirstChildElement(obj.c_str());
  int result = -1;
  if (el1) el1->QueryIntAttribute("value", &result);
  return result;
}

double Input::QueryDoubleAttribute(string obj, TiXmlElement* el) {
  TiXmlElement* el1 = el->FirstChildElement(obj.c_str());
  double result = -1.0;
  if (el1) el1->QueryDoubleAttribute("value", &result);
  return result;
}

double Input::QueryDoubleAttributeWithDefault(string obj, TiXmlElement* el, double def) {
  TiXmlElement* el1 = el->FirstChildElement(obj.c_str());
  double result = def;
  if (el1) el1->QueryDoubleAttribute("value", &result);
  return result;
}


int Input::LoadFile(const string inputfilename) {
   
  TiXmlDocument* doc = new   TiXmlDocument(inputfilename.c_str());
   
  TiXmlElement* el = NULL;
  //TiXmlElement* el1 = NULL;
  //TiXmlElement* el2 = NULL;
   
  if (doc->LoadFile(inputfilename.c_str(), TIXML_DEFAULT_ENCODING)) {
      
    string gridstr = "Healpix";     
    el = doc->FirstChildElement(gridstr.c_str());
 
    if (el){
	  string resstr = "resolution";	
      resolution = QueryDoubleAttribute(resstr, el);
	}
	cout << endl;      
      
    /*string PBHstr = "XraySpectrum";     
    el = doc->FirstChildElement(PBHstr.c_str());
 
    if (el){
	  string E0str = "E0";	
      E0 = QueryDoubleAttribute(E0str, el);
	  string Gammastr = "Gamma";	
      Gamma = QueryDoubleAttribute(Gammastr, el);
	  string Ecutstr = "Ecut";	
      Ecut = QueryDoubleAttribute(Ecutstr, el);
    }
	cout << endl;*/   

  }

  delete doc;
   
  cout << "Input file read successfully!" << endl;
   
  return 0;
}


void Input::Print() {

  //cout << endl;
  //cout << "************************************** " << endl;
  //cout << "You selected the following parameters: " << endl;
  //cout << "************************************** " << endl;
  //cout << endl;

}
