/*
 * Atz_XML_SAX_DataHandler.cpp
 *
 * Paul J. Atzberger
 * http://atzberger.org/
 *
 */

#include "Atz_XML_Helper_Handler_SkipNextTag.h"


Atz_XML_Helper_Handler_SkipNextTag::Atz_XML_Helper_Handler_SkipNextTag() {
  DataHandlerName = "SkipNextTag";
  DataHandlerType = "Atz_XML_Helper_Handler_SkipNextTag";
}

Atz_XML_Helper_Handler_SkipNextTag::~Atz_XML_Helper_Handler_SkipNextTag() {
  // TODO Auto-generated destructor stub
}

void Atz_XML_Helper_Handler_SkipNextTag::XML_startDocument(Atz_XML_SAX_DataHandler* sourceHandler) {

}

void Atz_XML_Helper_Handler_SkipNextTag::XML_endDocument(Atz_XML_SAX_DataHandler* sourceHandler) {
}

void Atz_XML_Helper_Handler_SkipNextTag::XML_startElement(string qName, Atz_XML::AttributesType *attributes, Atz_XML_SAX_DataHandler* sourceHandler) {
}

void Atz_XML_Helper_Handler_SkipNextTag::XML_characters(string xmlString, Atz_XML_SAX_DataHandler* sourceHandler) {
}

void Atz_XML_Helper_Handler_SkipNextTag::XML_endElement(string qName, Atz_XML_SAX_DataHandler* sourceHandler) {
}

void *Atz_XML_Helper_Handler_SkipNextTag::XML_getData() { /* gets data from parsing the XML */
  return NULL;
}


