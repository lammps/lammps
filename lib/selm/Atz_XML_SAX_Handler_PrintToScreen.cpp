/*

 Atz_XML_SAX_Handler_PrintToScreen.cpp
 
 Paul J. Atzberger
 http://atzberger.org/
 
*/

#include "Atz_XML_SAX_Handler_PrintToScreen.h"

Atz_XML_SAX_Handler_PrintToScreen::Atz_XML_SAX_Handler_PrintToScreen() {
  DataHandlerName = "PrintToScreen";
  DataHandlerType = "Atz_XML_SAX_Handler_PrintToScreen";
}

Atz_XML_SAX_Handler_PrintToScreen::~Atz_XML_SAX_Handler_PrintToScreen() {
  // TODO Auto-generated destructor stub
}


void Atz_XML_SAX_Handler_PrintToScreen::XML_startDocument(Atz_XML_SAX_DataHandler* sourceHandler) {

  cout << "XML Start Document" << endl;

}

void Atz_XML_SAX_Handler_PrintToScreen::XML_endDocument(Atz_XML_SAX_DataHandler* sourceHandler) {

  cout << "XML End Document" << endl;

}

void Atz_XML_SAX_Handler_PrintToScreen::XML_startElement(string qName, Atz_XML::AttributesType *attributes, Atz_XML_SAX_DataHandler* sourceHandler) {

  Atz_XML::AttributesType::iterator iter;

  cout << "XML Start Element : " << "TagName = \"" << qName << "\"" << endl;

  /* list the name, value pairs */
  for (iter = attributes->begin(); iter != attributes->end(); iter++) {
    cout << "                  ";
    cout << "AttrName = \"" << iter->first << "\"; ";
    cout << "AttrValue = \"" << iter->second << "\"" << endl;
  }

  xmlString_cur = "";

}

void Atz_XML_SAX_Handler_PrintToScreen::XML_characters(string xmlString_in, Atz_XML_SAX_DataHandler* sourceHandler) {

  xmlString_cur += xmlString_in;

  cout << "XML Characters : " << xmlString_cur << endl;

}

void Atz_XML_SAX_Handler_PrintToScreen::XML_endElement(string qName, Atz_XML_SAX_DataHandler* sourceHandler) {

  xmlString_last = xmlString_cur;
  xmlString_cur  = "";

  cout << "XML End Element : " << "TagName = \"" << qName << "\"" << endl;

}

void *Atz_XML_SAX_Handler_PrintToScreen::XML_getData() {

  string *xmlString_copy = new string(xmlString_last.c_str());

  cout << "XML Get Data Called" << endl;

  return xmlString_copy; /* this make the value safe for lists */
}
