#ifndef _INPUT_H
#define _INPUT_H

#include <iostream>
#include <string>
#include <vector>

#include "constants.h"

class TiXmlElement;

using namespace std;

class Input {
   
private:
   int QueryIntAttribute(string, TiXmlElement*);
   double QueryDoubleAttribute(string, TiXmlElement*);
   double QueryDoubleAttributeWithDefault(string, TiXmlElement*, double);
   
public:
   Input() { }  /**< The default constructor. */

   int LoadFile(const string);
   void Print();
   
   ~Input() { } /**< Destructor. Does nothing. */
   
   int resolution;
    
   //double E0;
   //double Gamma;
   //double Ecut;

};

#endif
