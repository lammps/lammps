/*

 Atz_XML_SAX_Handler_PrintToScreen.h
 
 Paul J. Atzberger
 http://atzberger.org/
 
*/

#ifndef ATZ_XML_SAX_HANDLER_PRINTTOSCREEN_H_
#define ATZ_XML_SAX_HANDLER_PRINTTOSCREEN_H_

using namespace std; /* ensures standard template library in namespace */

#include "Atz_XML.h"
#include "Atz_XML_SAX_DataHandler.h"

class Atz_XML_SAX_Handler_PrintToScreen : public Atz_XML_SAX_DataHandler {
public:
  Atz_XML_SAX_Handler_PrintToScreen();
  virtual
  ~Atz_XML_SAX_Handler_PrintToScreen();

protected:
  string xmlString_cur;
  string xmlString_last;

public:
  void XML_startDocument(Atz_XML_SAX_DataHandler* sourceHandler);

  void XML_endDocument(Atz_XML_SAX_DataHandler* sourceHandler);

  void XML_startElement(string qName, Atz_XML::AttributesType *attributes, Atz_XML_SAX_DataHandler* sourceHandler);

  void XML_characters(string xmlString, Atz_XML_SAX_DataHandler* sourceHandler);

  void XML_endElement(string qName, Atz_XML_SAX_DataHandler* sourceHandler);

  void *XML_getData(); /* gets data from parsing the XML */

};

#endif /* ATZ_XML_SAX_HANDLER_PRINTTOSCREEN_H_ */
