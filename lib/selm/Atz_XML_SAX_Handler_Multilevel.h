/*
 * Atz_XML_SAX_Handler_Multilevel.h
 *
 * Paul J. Atzberger
 * http://atzberger.org/
 *
 */

#ifndef ATZ_XML_SAX_HANDLER_MULTILEVEL_H_
#define ATZ_XML_SAX_HANDLER_MULTILEVEL_H_

#include <stack>
#include <vector>
#include <string>

#include "Atz_XML_SAX_DataHandler.h"
#include "Atz_XML.h"

using namespace std; /* ensures standard template library in namespace */
using namespace Atz_XML;

class Atz_XML_SAX_Handler_Multilevel : public Atz_XML_SAX_DataHandler {

public:
  static const int PARSE_MODE_NULL          = 0;
  static const int PARSE_MODE_NEXT_TAG      = 1;
  static const int PARSE_MODE_CURRENT_SCOPE = 2;

protected:

  int flagVerbose;

  vector<Atz_XML_SAX_DataHandler*> dataHandlerList;
  stack<Atz_XML_SAX_DataHandler*>  dataHandlerStack;

  vector<int>                      scopeDepthList;
  stack<int>                       scopeDepthStack;

  vector<int>                      parseModeList;
  stack<int>                       parseModeStack;

  Atz_XML_SAX_DataHandler         *lastPoppedHandler;

  int                              scopeDepthCount;

public:
  Atz_XML_SAX_Handler_Multilevel();
  Atz_XML_SAX_Handler_Multilevel(Atz_XML_SAX_DataHandler* currentHandler_in);
  virtual ~Atz_XML_SAX_Handler_Multilevel();

public:
  void setupGeneric();
  void clearAllStacks();
  void setLevelOfVerbosity(int flagVerbose_in); /* level of output to print, if any */
  void setDataHandler(Atz_XML_SAX_DataHandler* currentHandler_in);
  void changeCurrentDataHandler(Atz_XML_SAX_DataHandler* newCurrentHandler_in);
  Atz_XML_SAX_DataHandler *getCurrentDataHandler();
  stack<Atz_XML_SAX_DataHandler*> *getDataHandlerStack();
  int getCurrentParseMode();
  const char *getCurrentParseModeStr();
  Atz_XML_SAX_DataHandler *getLastUsedDataHandler();
  Atz_XML_SAX_DataHandler *getLastPoppedDataHandler();
  void parseNextTagWithDataHandler(Atz_XML_SAX_DataHandler* currentHandler);
  void parseCurrentScopeWithDataHandler(Atz_XML_SAX_DataHandler* currentHandler);
  void pushDataHandler(Atz_XML_SAX_DataHandler* currentHandler);
  void pushDataHandler(Atz_XML_SAX_DataHandler* currentHandler, int parseMode);
  Atz_XML_SAX_DataHandler *popDataHandler();
  bool isEmptyDataHandlerStack();
  Atz_XML_SAX_DataHandler *peekDataHandler();
  void printCallInfo(const char *callName);


  /* ================================================================================= */
  /* ================================ SAX Event Handlers ============================= */
  void startDocument();
  void endDocument();
  void startElement(string qName, AttributesType *attributes);
  void characters(string xmlString);
  void endElement(string qName);

  /* ================================ SAX Event Handlers ============================= */
  /* ================================================================================= */


  /* ================================================================================= */
  /* ================== Atz_XML_SAX_DataHandler wrapper ===================== */
  virtual void  XML_startDocument(Atz_XML_SAX_DataHandler* sourceHandler);
  virtual void  XML_endDocument(Atz_XML_SAX_DataHandler* sourceHandler);
  virtual void  XML_startElement(string qName, AttributesType *attributes, Atz_XML_SAX_DataHandler* sourceHandler);
  virtual void  XML_endElement(string qName, Atz_XML_SAX_DataHandler* sourceHandler);
  virtual void  XML_characters(string xmlString, Atz_XML_SAX_DataHandler* sourceHandler);
  virtual void *XML_getData();

};

#endif /* ATZ_XML_SAX_HANDLER_MULTILEVEL_H_ */
