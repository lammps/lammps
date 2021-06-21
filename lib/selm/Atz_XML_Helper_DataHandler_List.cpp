/*
 * Atz_XML_Helper_DataHandler_List.cpp
 *
 * Paul J. Atzberger
 * http://atzberger.org/
 *
 */

#include "Atz_XML_Helper_DataHandler_List.h"

Atz_XML_Helper_DataHandler_List::Atz_XML_Helper_DataHandler_List() {
  setupGeneric();
}

Atz_XML_Helper_DataHandler_List::Atz_XML_Helper_DataHandler_List(Atz_XML_SAX_DataHandler *dataHandler_in) {
  setupGeneric();
  setDataHandler(dataHandler_in); /* set up current data handler */
  scopeDepthCount = 0;
}

Atz_XML_Helper_DataHandler_List::~Atz_XML_Helper_DataHandler_List() {
  // TODO Auto-generated destructor stub
}

void Atz_XML_Helper_DataHandler_List::setupGeneric() {
  DataHandlerName = "ListHandler";
  DataHandlerType = "Atz_XML_Helper_DataHandler_List";

  tagDataLists = new DataListMapType();
}

Atz_XML_Helper_DataHandler_List::DataListMapType *Atz_XML_Helper_DataHandler_List::getTagDataLists() {
  return tagDataLists;
}

void Atz_XML_Helper_DataHandler_List::XML_startDocument(Atz_XML_SAX_DataHandler* sourceHandler) {
  this->Atz_XML_SAX_Handler_Multilevel::XML_startDocument(sourceHandler);
}

void Atz_XML_Helper_DataHandler_List::XML_endDocument(Atz_XML_SAX_DataHandler* sourceHandler) {
  this->Atz_XML_SAX_Handler_Multilevel::XML_endDocument(sourceHandler);
}

void Atz_XML_Helper_DataHandler_List::XML_startElement(string qName, Atz_XML::AttributesType *attributes, Atz_XML_SAX_DataHandler* sourceHandler) {

  if (scopeDepthCount == 0) { /* for zero scope setup the parsing behavior for this next tag */

    if (dataHandlerStack.empty()) { /* this indicates tag was recently processed and handler popped */
      /* reset to use the last parser to process the next tag */
      this->Atz_XML_SAX_Handler_Multilevel::parseNextTagWithDataHandler(getLastUsedDataHandler());
    } else {
      /* otherwise current data handler which should already be set will be used */
      /* this case will usually occur for the first tag */
      int a = 0;
    }

  } /* for scope of greater depth nothing special needs to be done */

  this->Atz_XML_SAX_Handler_Multilevel::XML_startElement(qName, attributes, sourceHandler);

}

void Atz_XML_Helper_DataHandler_List::XML_characters(string xmlString, Atz_XML_SAX_DataHandler* sourceHandler) {
  this->Atz_XML_SAX_Handler_Multilevel::XML_characters(xmlString, sourceHandler);
}

void Atz_XML_Helper_DataHandler_List::XML_endElement(string qName, Atz_XML_SAX_DataHandler* sourceHandler) {

  this->Atz_XML_SAX_Handler_Multilevel::XML_endElement(qName, sourceHandler);

  /* only record lists when the scope is zero (relative to this listDataHandler parser) */
  if (scopeDepthCount == 0) {
    ObjListType               *dataHandlerList = NULL;
    DataListMapType::iterator  dataIterator;
    void                      *XMLData;
    Atz_XML_SAX_DataHandler   *lastDataHandler;

    /* check if current tag type was already encountered */
    dataIterator = tagDataLists->find(qName);
    if (dataIterator == tagDataLists->end()) { /* indicates key not found, new data type */
      dataHandlerList = new ObjListType();    /* create new list */
    } else { /* otherwise add to the existing list */
      dataHandlerList = dataIterator->second; /* dereference value */
      //tagDataLists->remove(qName); /* remove the list, (added back below) */
    }

    lastDataHandler = getLastUsedDataHandler();
    XMLData = lastDataHandler->XML_getData();
    dataHandlerList->push_back(XMLData); /* add parsed data to the list */
    tagDataLists->insert(DataListMapPairType(qName, dataHandlerList));  /* save the list for later */
  }

  if (scopeDepthCount == -1) {
    int a = 1;
  }

}

void *Atz_XML_Helper_DataHandler_List::XML_getData() { /* gets data from parsing the XML */
  return tagDataLists; /* TODO: Need to make as clone!!!! */ /* really make a copy of the tagDataLists */
}
