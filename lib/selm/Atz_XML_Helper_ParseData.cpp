/*

 Atz_XML_Helper_ParseData.cpp
 
 Paul J. Atzberger
 http://atzberger.org/
 
*/

#include "Atz_XML_Helper_ParseData.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


Atz_XML_Helper_ParseData::Atz_XML_Helper_ParseData() {
  // TODO Auto-generated constructor stub
}


Atz_XML_Helper_ParseData::~Atz_XML_Helper_ParseData() {
  // TODO Auto-generated destructor stub
}


int Atz_XML_Helper_ParseData::getIntFromAttr(Atz_XML::AttributesType *attributes) {
  return getIntFromAttr("value", attributes);
}


int Atz_XML_Helper_ParseData::getIntFromAttr(const char *fieldName, Atz_XML::AttributesType *attributes) {
  /* get the values from the attributes */
  /* convert from valueStr to integer value */
  string valueStr;
  int    value;

  valueStr = attributes->find(fieldName)->second;
  value    = atoi(valueStr.c_str());

  return value;
}


double Atz_XML_Helper_ParseData::getDoubleFromAttr(Atz_XML::AttributesType *attributes) {
  return getDoubleFromAttr("value", attributes);
}


double Atz_XML_Helper_ParseData::getDoubleFromAttr(const char *fieldName, Atz_XML::AttributesType *attributes) {

  /* get the values from the attributes */
  /* convert from valueStr to integer value */
  string valueStr;
  double value;

  valueStr = attributes->find(fieldName)->second;
  value    = atof(valueStr.c_str());

  return value;
}

char *Atz_XML_Helper_ParseData::getCStringFromAttr(Atz_XML::AttributesType *attributes) {

  string *cpp_str = getStringFromAttr("value", attributes);
  char   *c_return;

  const char *c_str = cpp_str->c_str();
  int len           = strlen(c_str);
  c_return          = (char *)malloc(sizeof(char)*(len + 1));
  strcpy(c_return, c_str);

  return c_return;
}

string *Atz_XML_Helper_ParseData::getStringFromAttr(Atz_XML::AttributesType *attributes) {
  return getStringFromAttr("value", attributes);
}


string *Atz_XML_Helper_ParseData::getStringFromAttr(const char *fieldName, Atz_XML::AttributesType *attributes) {
  return &attributes->find(fieldName)->second;
}


void Atz_XML_Helper_ParseData::getDoubleArrayFromAttr(Atz_XML::AttributesType *attributes, double **doubleArray_ptr, int *numDoubleArray_ptr) {
  getDoubleArrayFromAttr("value", attributes, doubleArray_ptr, numDoubleArray_ptr);
}

void Atz_XML_Helper_ParseData::getDoubleArrayFromAttr(const char *fieldName, Atz_XML::AttributesType *attributes, double **doubleArray_ptr, int *numDoubleArray_ptr) {

  /* get the values from the attributes */
  /* convert from valueStr to integer value */
  string  valueStr;
  double *doubleArray    = *doubleArray_ptr;
  int     numDoubleArray = *numDoubleArray_ptr;

  valueStr = attributes->find(fieldName)->second;
  parseDoubleArrayFromString(valueStr, &doubleArray, &numDoubleArray);

  /* return values */
  (*doubleArray_ptr)    = doubleArray;
  (*numDoubleArray_ptr) = numDoubleArray;

}


void Atz_XML_Helper_ParseData::parseDoubleArrayFromString(string dataStr, double **doubleArray_ptr, int *numDoubleArray_ptr) {

  const char *error_str_code = "Atz_XML_Helper_ParseData.cpp";
  const char *error_str_func = "parseDoubleArrayFromString()";

  /* WARNING: Be sure when using this routine to set
   * *doubleArray_ptr = NULL if you want the routine
   * to allocated the array.  Otherwise, the routine
   * uses the pre-allocated array provided by
   * *doubleArray_ptr.
   */

  /* create string stream... */
  stringstream       dataStream(dataStr);
  string             token;

  vector<double>     doubleList;
  double             value;
  double            *doubleArray;

  int                N;

  /* parse the doubles */
  while (dataStream.good()) {
    token.clear();
    dataStream >> token;
    if (token != "") {
      value = atof(token.c_str());
      doubleList.push_back(value);
    }
  }

  /* construct the double array */
  N        = doubleList.size();
  if ((*doubleArray_ptr) == NULL) {
    if (N != 0) {
      doubleArray = (double *)malloc(sizeof(double)*N);
    } else {
      doubleArray = NULL;
    }
  } else {
    doubleArray = (*doubleArray_ptr);

    if ((*numDoubleArray_ptr) < N) {
      /* issue ERROR array too small!!! */
      stringstream message;
      message << "doubleArray_ptr pre-allocated by user but array not large enough." << endl;
      message << "(*numDoubleArray_ptr) = " << (*numDoubleArray_ptr) << " < " << N << endl;
      message << endl;
      SELM_Package::packageError(error_str_code, error_str_func, message);
    }

  }

  for (int i = 0; i < N; i++) {
    doubleArray[i] = doubleList[i];
  }

  /* return the results */
  (*doubleArray_ptr)    = doubleArray;
  (*numDoubleArray_ptr) = N;

}




void Atz_XML_Helper_ParseData::parseIntArrayFromString(string dataStr, int **intArray_ptr, int *numIntArray_ptr) {

  const char *error_str_code = "Atz_XML_Helper_ParseData.cpp";
  const char *error_str_func = "parseIntArrayFromString()";

  /* WARNING: Be sure when using this routine to set
   * *intArray_ptr = NULL if you want the routine
   * to allocated the array.  Otherwise, the routine
   * uses the pre-allocated array provided by
   * *intArray_ptr.
   */

  /* create string stream... */
  stringstream    dataStream(dataStr);
  string          token;

  vector<int>     intList;
  int             value;
  int            *intArray;

  int             N;

  /* parse the ints */
  while (dataStream.good()) {
    token.clear();
    dataStream >> token;
    if (token != "") {
      value = atoi(token.c_str());
      intList.push_back(value);
    }
  }

  /* construct the int array */
  N        = intList.size();
  if ((*intArray_ptr) == NULL) {
    if (N != 0) {
      intArray = (int *)malloc(sizeof(int)*N);
    } else {
      intArray = NULL;
    }
  } else {
    intArray = (*intArray_ptr);

    if ((*numIntArray_ptr) < N) {
      /* issue ERROR array too small!!! */

      stringstream message;
      message << "ERROR: Atz_XML_Helper_ParseData::parseIntArrayFromString" << endl;
      message << "intArray_ptr pre-allocated by user but array not large enough." << endl;
      message << "(*numIntArray_ptr) = " << (*numIntArray_ptr) << " < " << N << endl;
      message << endl;
      SELM_Package::packageError(error_str_code, error_str_func, message);

    }

  }

  for (int i = 0; i < N; i++) {
    intArray[i] = intList[i];
  }

  /* return the results */
  (*intArray_ptr)    = intArray;
  (*numIntArray_ptr) = N;

}



