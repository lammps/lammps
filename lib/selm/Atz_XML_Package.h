/* 
 * XML Handling
 * 
 * Paul J. Atzberger
 * http://atzberger.org/
 *
*/

#ifndef ATZ_XML_PACKAGE_H
#define ATZ_XML_PACKAGE_H

#include "Atz_XML.h"

#include "Atz_XML_Parser.h"
#include "Atz_XML_SAX_DataHandler.h"
#include "Atz_XML_SAX_Handler_Multilevel.h"
#include "Atz_XML_SAX_Handler_PrintToScreen.h"
#include "Atz_XML_Helper_ParseData.h"
#include "Atz_XML_Helper_DataHandler_List.h"
#include "Atz_XML_Helper_Handler_SkipNextTag.h"

#include "SELM_Package.h"

#include <sstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>

using namespace std; /* ensures standard template library names used */

class Atz_XML_Package {

 public:

  /* ================ Constants ================= */

  /* ================ Data structure type definitions ================= */

  /* ================ Variables ================= */

  /* ================ Function prototypes ================= */
  Atz_XML_Package();
  virtual ~Atz_XML_Package();

  static void writeXMLHeader(FILE *fid);
  static void writeXMLHeader(ofstream &fid);

  static void writeTagStart(FILE *fid, string tagName);
  static void writeTagStart(FILE *fid, string tagName, string extras);

  static void writeTagStart(FILE *fid, const char *tagName);
  static void writeTagStart(FILE *fid, const char *tagName, const char *extras);

  static void writeTagStart(ofstream &fid, string tagName);
  static void writeTagStart(ofstream &fid, string tagName, string extras);

  static void writeTagStart(ofstream &fid, const char *tagName);
  static void writeTagStart(ofstream &fid, const char *tagName, const char *extras);

  static void writeTagEnd(FILE *fid, string tagName);
  static void writeTagEnd(FILE *fid, const char *tagName);

  static void writeTagEnd(ofstream &fid, string tagName);
  static void writeTagEnd(ofstream &fid, const char *tagName);

  static void writeTagValueDouble(FILE *fid, const char *tagName, double value);
  static void writeTagValueDoubleArray(FILE *fid, const char *tagName, int N, double *values);

  static void writeTagValueDouble(ofstream &fid, const char *tagName, double value);
  static void writeTagValueDoubleArray(ofstream &fid, const char *tagName, int N, double *values);

  static void packageError(const char *error_str_code, const char *error_str_func, stringstream &message);
  static void packageError(const char *error_str_code, const char *error_str_func, string &message);
  static void packageError(const char *error_str_code, const char *error_str_func, const char *message);

};

#endif
