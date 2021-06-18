/*
 * Atz_XML_SAX_DataHandler.h
 *
 * Paul J. Atzberger
 * http://atzberger.org/
 *
 */

#ifndef ATZ_XML_HELPER_HANDLER_SKIPNEXTTAG_H_
#define ATZ_XML_HELPER_HANDLER_SKIPNEXTTAG_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>

using namespace std; /* ensures standard template library in namespace */

#include "Atz_XML_SAX_DataHandler.h"
#include "Atz_XML.h"

class Atz_XML_Helper_Handler_SkipNextTag : public Atz_XML_SAX_DataHandler {

  public:
  Atz_XML_Helper_Handler_SkipNextTag();
  virtual ~Atz_XML_Helper_Handler_SkipNextTag();

  public:
    virtual void XML_startDocument(Atz_XML_SAX_DataHandler* sourceHandler);

    virtual void XML_endDocument(Atz_XML_SAX_DataHandler* sourceHandler);

    virtual void XML_startElement(string qName, Atz_XML::AttributesType *attributes, Atz_XML_SAX_DataHandler* sourceHandler);

    virtual void XML_characters(string xmlString, Atz_XML_SAX_DataHandler* sourceHandler);

    virtual void XML_endElement(string qName, Atz_XML_SAX_DataHandler* sourceHandler);

    virtual void *XML_getData(); /* gets data from parsing the XML */

};

#endif /* ATZ_XML_HELPER_HANDLER_SKIPNEXTTAG_H_ */
