/*

  Atz_XML_Parser.h
 
  Paul J. Atzberger
  http://atzberger.org/
 
*/

#ifndef ATZ_XML_PARSER_H_
#define ATZ_XML_PARSER_H_

using namespace std; /* ensures standard template library in namespace */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>

#include "Atz_XML.h"
#include "Atz_XML_Package.h"
#include "Atz_XML_SAX_DataHandler.h"

class Atz_XML_Parser {

public:
  Atz_XML_Parser();
  virtual
  ~Atz_XML_Parser();

public:

  static const char *error_str_code;

  static void   parse(string fileName, Atz_XML_SAX_DataHandler *dataHandler);
  static void   parse(const char *fileName, Atz_XML_SAX_DataHandler *dataHandler);
  static void   parse(ifstream *fileStream, Atz_XML_SAX_DataHandler *dataHandler);

protected:
  static void   processTag(string tagStr);
  static void   removeLeadingWhiteSpace(ifstream* strStream);
  static void   getTagName(ifstream *tagStream, string& tagName);
  static void   parseUntilEndComment(ifstream *tagStream);
  static void   getTagAttributes(ifstream *tagStream, Atz_XML::AttributesType *attributes);
  static void   getAttrName(ifstream* attrStream, string& attrName);
  static void   getAttrValue(ifstream* attrStream, string& attrValue);
  static string getPassableName(string tagName);

};

#endif /* ATZ_XML_PARSER_H_ */
