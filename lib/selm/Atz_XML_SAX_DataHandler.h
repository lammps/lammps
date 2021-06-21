/*
 * Atz_XML_SAX_DataHandler.h
 *
 * Paul J. Atzberger
 * http://atzberger.org/
 *
 * Abstract class defining the interface for data handling.
 *
 */

#ifndef ATZ_XML_SAX_DATAHANDLER_H_
#define ATZ_XML_SAX_DATAHANDLER_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>

using namespace std; /* ensures standard template library in namespace */

#include "Atz_XML.h"

class Atz_XML_SAX_DataHandler {

  public:
    string  DataHandlerName;
    string  DataHandlerType;

    void   *extras;  /* extra data used for parsing */

  public:
    Atz_XML_SAX_DataHandler();
    virtual ~Atz_XML_SAX_DataHandler();

  public:
    virtual void XML_startDocument(Atz_XML_SAX_DataHandler* sourceHandler) = 0;

    virtual void XML_endDocument(Atz_XML_SAX_DataHandler* sourceHandler) = 0;

    virtual void XML_startElement(string qName, Atz_XML::AttributesType *attributes, Atz_XML_SAX_DataHandler* sourceHandler) = 0;

    virtual void XML_characters(string xmlString, Atz_XML_SAX_DataHandler* sourceHandler) = 0;

    virtual void XML_endElement(string qName, Atz_XML_SAX_DataHandler* sourceHandler) = 0;

    virtual void *XML_getData() = 0; /* gets data from parsing the XML */

};

#endif /* ATZ_XML_SAX_DATAHANDLER_H_ */
