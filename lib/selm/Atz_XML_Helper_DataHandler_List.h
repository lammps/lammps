/*
 * Atz_XML_Helper_DataHandler_List.h
 *
 * Paul J. Atzberger
 * http://atzberger.org/
 *
 */

#ifndef ATZ_XML_HELPER_DATAHANDLER_LIST_H_
#define ATZ_XML_HELPER_DATAHANDLER_LIST_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>

#include "Atz_XML_Package.h"

using namespace std; /* ensures standard template library in namespace */
using namespace Atz_XML; /* ensures standard template library in namespace */

class Atz_XML_Helper_DataHandler_List : public Atz_XML_SAX_Handler_Multilevel {

public:
  typedef vector<void*>             ObjListType;
  typedef map<string,ObjListType*>  DataListMapType;
  typedef pair<string,ObjListType*> DataListMapPairType;

protected:
  //Atz_XML_SAX_DataHandlerInterface dataHandler  = null;
  DataListMapType    *tagDataLists;
  string              xmlString;
  AttributesType      xmlAttributes;

  public:
  Atz_XML_Helper_DataHandler_List();
  Atz_XML_Helper_DataHandler_List(Atz_XML_SAX_DataHandler* dataHandler_in);
  virtual ~Atz_XML_Helper_DataHandler_List();

  protected:
    void setupGeneric();

  public:
    DataListMapType *getTagDataLists();

    virtual void XML_startDocument(Atz_XML_SAX_DataHandler* sourceHandler);

    virtual void XML_endDocument(Atz_XML_SAX_DataHandler* sourceHandler);

    virtual void XML_startElement(string qName, Atz_XML::AttributesType *attributes, Atz_XML_SAX_DataHandler* sourceHandler);

    virtual void XML_characters(string xmlString, Atz_XML_SAX_DataHandler* sourceHandler);

    virtual void XML_endElement(string qName, Atz_XML_SAX_DataHandler* sourceHandler);

    virtual void *XML_getData(); /* gets data from parsing the XML */

};

#endif /* ATZ_XML_HELPER_DATAHANDLER_LIST_H_ */
