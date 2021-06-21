/*

 Atz_XML_Helper_ParseData.h
 
 Paul J. Atzberger
 http://atzberger.org/
  
 Collection of routines used for help
 in parsing the XML data sets.
 
*/

#ifndef ATZ_XML_HELPER_PARSEDATA_H_
#define ATZ_XML_HELPER_PARSEDATA_H_


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

#include "Atz_XML_Package.h"
#include "SELM_Package.h"

using namespace std; /* ensures standard template library in namespace */
using namespace Atz_XML; /* ensures standard template library in namespace */
using namespace LAMMPS_NS;

class Atz_XML_Helper_ParseData  {

public:

  /* get data from attributes */
  static int     getIntFromAttr(Atz_XML::AttributesType *attributes);
  static int     getIntFromAttr(const char *fieldName, Atz_XML::AttributesType *attributes);

  static double  getDoubleFromAttr(Atz_XML::AttributesType *attributes);
  static double  getDoubleFromAttr(const char *fieldName, Atz_XML::AttributesType *attributes);

  static char *getCStringFromAttr(Atz_XML::AttributesType *attributes);

  static string *getStringFromAttr(Atz_XML::AttributesType *attributes);
  static string *getStringFromAttr(const char *fieldName, Atz_XML::AttributesType *attributes);

  static void getDoubleArrayFromAttr(Atz_XML::AttributesType *attributes, double **doubleArray_ptr, int *numDoubleArray_ptr);
  static void getDoubleArrayFromAttr(const char *fieldName, Atz_XML::AttributesType *attributes, double **doubleArray_ptr, int *numDoubleArray_ptr);

  static void parseDoubleArrayFromString(string dataStr, double **doubleArray_ptr, int *numDoubleArray_ptr);
  static void parseIntArrayFromString(string dataStr, int **intArray_ptr, int *numIntArray_ptr);

protected:

  public:
  Atz_XML_Helper_ParseData();
  virtual ~Atz_XML_Helper_ParseData();

  protected:

  public:


};

#endif /* ATZ_XML_HELPER_PARSEDATA_H_ */
