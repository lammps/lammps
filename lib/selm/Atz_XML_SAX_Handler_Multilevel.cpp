/*
 * Atz_XML_SAX_Handler_Multilevel.cpp
 *
 * Paul J. Atzberger
 * http://atzberger.org/
 *
 */

#include "Atz_XML_SAX_Handler_Multilevel.h"

Atz_XML_SAX_Handler_Multilevel::Atz_XML_SAX_Handler_Multilevel() {
  setupGeneric();
}


Atz_XML_SAX_Handler_Multilevel::Atz_XML_SAX_Handler_Multilevel(Atz_XML_SAX_DataHandler* currentHandler_in) {
  setupGeneric();
  setDataHandler(currentHandler_in);
  scopeDepthCount = 0;
}


Atz_XML_SAX_Handler_Multilevel::~Atz_XML_SAX_Handler_Multilevel() {
  // TODO Auto-generated destructor stub
}


void Atz_XML_SAX_Handler_Multilevel::setupGeneric() {

  //dataHandlerStack(dataHandlerList);
  //scopeDepthStack(scopeDepthList);
  //parseModeStack(parseModeList);
  DataHandlerName = "Multilevel";
  DataHandlerType = "Atz_XML_SAX_Handler_Multilevel";

  scopeDepthCount   = 0;
  flagVerbose       = 0; /* show output by default */
}


void Atz_XML_SAX_Handler_Multilevel::clearAllStacks() {

  /* pop off all elements until empty */
  while (!dataHandlerStack.empty()) {
    dataHandlerStack.pop();
  }

  while (!scopeDepthStack.empty()) {
    scopeDepthStack.pop();
  }

  while (!parseModeStack.empty()) {
    parseModeStack.pop();
  }

}

void Atz_XML_SAX_Handler_Multilevel::setLevelOfVerbosity(int flagVerbose_in) {
  flagVerbose = flagVerbose_in;
}

void Atz_XML_SAX_Handler_Multilevel::setDataHandler(Atz_XML_SAX_DataHandler* currentHandler_in) {
  clearAllStacks();
  pushDataHandler(currentHandler_in, PARSE_MODE_CURRENT_SCOPE); /* old way, no explicit parse mode */
  lastPoppedHandler = currentHandler_in; /* treat current as the last used (last popped) */
}


void Atz_XML_SAX_Handler_Multilevel::changeCurrentDataHandler(Atz_XML_SAX_DataHandler* newCurrentHandler_in) {
  popDataHandler();
  pushDataHandler(newCurrentHandler_in);
}


Atz_XML_SAX_DataHandler *Atz_XML_SAX_Handler_Multilevel::getCurrentDataHandler() {
  return (dataHandlerStack.top());
}


stack<Atz_XML_SAX_DataHandler*> *Atz_XML_SAX_Handler_Multilevel::getDataHandlerStack() {
  return &dataHandlerStack;
}


int Atz_XML_SAX_Handler_Multilevel::getCurrentParseMode() {
  return parseModeStack.top();
}

const char *Atz_XML_SAX_Handler_Multilevel::getCurrentParseModeStr() {

  switch (getCurrentParseMode()) {
  case PARSE_MODE_NEXT_TAG:
    return "PARSE_MODE_NEXT_TAG";
    break;
  case PARSE_MODE_CURRENT_SCOPE:
    return "PARSE_MODE_CURRENT_SCOPE";
    break;
  default:
    break;
  }

  return "UNKNOWN";
}



/* Returns the last data handler used to parse data */
Atz_XML_SAX_DataHandler *Atz_XML_SAX_Handler_Multilevel::getLastUsedDataHandler() {
  return getLastPoppedDataHandler();
}


Atz_XML_SAX_DataHandler *Atz_XML_SAX_Handler_Multilevel::getLastPoppedDataHandler() {
  return lastPoppedHandler;
}


/* Parses the next tag encountered at the current scope with the specified handler.
 */
void Atz_XML_SAX_Handler_Multilevel::parseNextTagWithDataHandler(Atz_XML_SAX_DataHandler* currentHandler) {
  pushDataHandler(currentHandler, PARSE_MODE_NEXT_TAG);
}


/* Parses all data at the current scope with the specified handler.
 */
void Atz_XML_SAX_Handler_Multilevel::parseCurrentScopeWithDataHandler(Atz_XML_SAX_DataHandler* currentHandler) {
  pushDataHandler(currentHandler, PARSE_MODE_CURRENT_SCOPE);
}


void Atz_XML_SAX_Handler_Multilevel::pushDataHandler(Atz_XML_SAX_DataHandler* currentHandler) {
  pushDataHandler(currentHandler, PARSE_MODE_NEXT_TAG); /* default parse mode used */
}


void Atz_XML_SAX_Handler_Multilevel::pushDataHandler(Atz_XML_SAX_DataHandler* currentHandler, int parseMode) {

  if (flagVerbose > 0) {
    cout << endl;
    printCallInfo("pushDataHandler():");
    cout << "Current scope = " << scopeDepthCount << endl;
    if (!dataHandlerStack.empty()) {
      cout << "Current data handler had name = " << getCurrentDataHandler()->DataHandlerName << endl;
      cout << "Current data handler had type = " << getCurrentDataHandler()->DataHandlerType << endl;
    } else {
      cout << "Stack is currently empty" << endl;
    }
    cout << "Pushing onto the stack" << endl;
  }

  dataHandlerStack.push(currentHandler);
  scopeDepthStack.push(scopeDepthCount);
  parseModeStack.push(parseMode);

  if (flagVerbose > 0) {
    cout << "New data handler is now name = " << getCurrentDataHandler()->DataHandlerName << endl;
    cout << "New data handler is now type = " << getCurrentDataHandler()->DataHandlerType << endl;
    cout << "New parser mode is now = " << getCurrentParseModeStr() << endl;
  }

}


Atz_XML_SAX_DataHandler *Atz_XML_SAX_Handler_Multilevel::popDataHandler() {

  int lastPoppedScope;
  int nowCurScope;

  if (flagVerbose > 0) {
    cout << endl;
    cout << "Atz_XML_SAX_Handler_Multilevel : popDataHandler():" << endl;
    cout << "this->DataHandlerName = " << this->DataHandlerName;
    cout << " : this->DataHandlerType = " << this->DataHandlerType;
    cout << "Last was scope = " << scopeDepthCount << endl;
    cout << "Last data handler had name = " << getCurrentDataHandler()->DataHandlerName << endl;
    cout << "Last data handler had type = " << getCurrentDataHandler()->DataHandlerType << endl;
    cout << "Popping the stack" << endl;
  }

  if (dataHandlerStack.empty() == false) {
    parseModeStack.pop();
    lastPoppedScope   = scopeDepthStack.top();
    scopeDepthStack.pop();
    nowCurScope       = scopeDepthStack.top();
    lastPoppedHandler = dataHandlerStack.top();
    dataHandlerStack.pop();
  } else {
    scopeDepthCount   = -1;
    lastPoppedHandler = NULL;
  }

  if (flagVerbose > 0) {
    cout << "Current data handler is now name = " << getCurrentDataHandler()->DataHandlerName << endl;
    cout << "Current data handler is now type = " << getCurrentDataHandler()->DataHandlerType << endl;
    cout << "Current scope should be = " << lastPoppedScope << endl;
    cout << "Current parser mode is now = " << getCurrentParseModeStr() << endl;
  }

  return lastPoppedHandler;
}

bool Atz_XML_SAX_Handler_Multilevel::isEmptyDataHandlerStack() {
  return dataHandlerStack.empty();
}


Atz_XML_SAX_DataHandler *Atz_XML_SAX_Handler_Multilevel::peekDataHandler() {
  return dataHandlerStack.top();
}


/* ================================================================================= */
/* ================================ SAX Event Handlers ============================= */

void Atz_XML_SAX_Handler_Multilevel::startDocument() {
  scopeDepthCount = 0;

  if (flagVerbose > 0) {
    cout << endl;
    printCallInfo("startDocument()");
    cout << "Current scope = " << scopeDepthCount << endl;
    cout << "Current data handler name = " << getCurrentDataHandler()->DataHandlerName << endl;
    cout << "Current data handler type = " << getCurrentDataHandler()->DataHandlerType << endl;
    cout << "Calling data handler XML_startDocument()" << endl;
  }

  getCurrentDataHandler()->XML_startDocument(this);
}

void Atz_XML_SAX_Handler_Multilevel::endDocument() {

  if (flagVerbose > 0) {
    cout << endl;
    printCallInfo("endDocument()");
    cout << "Current scope = " << scopeDepthCount << endl;
  }

  if (!isEmptyDataHandlerStack()) {

    if (flagVerbose > 0) {
      cout << "Current data handler name = " << getCurrentDataHandler()->DataHandlerName << endl;
      cout << "Current data handler type = " << getCurrentDataHandler()->DataHandlerType << endl;
      cout << "Calling data handler XML_endDocument()" << endl;
    }

    getCurrentDataHandler()->XML_endDocument(this);

  } else {
    if (flagVerbose > 0) {
      cout << "WARNING: Data handler stack empty so no calls to specific handler made." << endl;
    }
  }

  if (scopeDepthCount != 0) {
    if (flagVerbose > 0) {
      cout << endl;
      printCallInfo("endDocument()");
      cout << "Current scope = " << scopeDepthCount << endl;
      cout << "Current data handler type = " << getCurrentDataHandler()->DataHandlerType << endl;
      cout << "Calling data handler XML_endDocument()" << endl;
    } else {
      cout << "Atz_XML_SAX_Handler_Multilevel : endDocument():" << endl;
    }
    cout << "WARNING: Scope depth is not zero at end of file." << endl;
  }

}

//Event Handlers
void Atz_XML_SAX_Handler_Multilevel::startElement(string qName, AttributesType *attributes) {
  scopeDepthCount++;

  if (flagVerbose > 0) {
    cout << endl;
    printCallInfo("startElement()");
    cout << "Tag name = " << qName << endl;
    cout << "Parser mode = " << getCurrentParseModeStr() << endl;
    cout << "Current scope = " << scopeDepthCount << endl;
    cout << "Current data handler name = " << getCurrentDataHandler()->DataHandlerName << endl;
    cout << "Current data handler type = " << getCurrentDataHandler()->DataHandlerType << endl;
    cout << "Calling data handler XML_startElement()" << endl;
    /* @@@ print attributes here... */
  }

  getCurrentDataHandler()->XML_startElement(qName, attributes, this);
}

void Atz_XML_SAX_Handler_Multilevel::characters(string xmlString) {


  if (!isEmptyDataHandlerStack()) { /* pass characters along if data handler is specified */

    if (flagVerbose > 0) {
      cout << endl;
      printCallInfo("characters()");
      cout << "Current data handler name = " << getCurrentDataHandler()->DataHandlerName << endl;
      cout << "Current data handler type = " << getCurrentDataHandler()->DataHandlerType << endl;
      cout << "String of characters to process = " << xmlString << endl;
      cout << "Calling data handler XML_characters()" << endl;
    }

    getCurrentDataHandler()->XML_characters(xmlString, this);
  } else {
    /* ignore these XML characters possibly between tags */
    if (flagVerbose > 0) {
      cout << endl;
      printCallInfo("characters()");
      cout << "WARNING: Data handler stack empty so characters ignored." << endl;
    }
  }
}


void Atz_XML_SAX_Handler_Multilevel::endElement(string qName) {

  int lastScopeStacked = -1;

  /* process tag depending on mode */
  switch (getCurrentParseMode()) {

  case PARSE_MODE_NEXT_TAG:

    /* Parse until completion of the tag occuring just after the push.
     * This results in the pushed handler only operating on one data tag
     * and then returning control when the tag is completed.
     */

    if (scopeDepthStack.empty() == false) {
      lastScopeStacked = scopeDepthStack.top();
    }

    if (flagVerbose > 0) {
      cout << endl;
      printCallInfo("endElement()");
      cout << "Tag name = " << qName << endl;
      cout << "Parser mode = " << getCurrentParseModeStr() << endl;
      cout << "Current scope = " << scopeDepthCount << endl;
      cout << "Last scope stacked = " << lastScopeStacked << endl;
      cout << "Parser mode = PARSE_MODE_NEXT_TAG" << endl;
      cout << "Current data handler name = " << getCurrentDataHandler()->DataHandlerName << endl;
      cout << "Current data handler type = " << getCurrentDataHandler()->DataHandlerType << endl;
      cout << "Calling data handler XML_endElement()" << endl;
    }

    /* process the tag with the current data handler */
    getCurrentDataHandler()->XML_endElement(qName, this);

    /* check if we are about to return to the reference scope just before the push of dataHandler */
    if (scopeDepthCount <= lastScopeStacked + 1) { /* scope is about to return from child processing */
      popDataHandler(); /* pop the handler back to the original one before the push */
    }

    break; /* end PARSE_MODE_NEXT_TAG */

  case PARSE_MODE_CURRENT_SCOPE:

    /* Parse until an attempt is made to exit the current scope.
     * This allows for the pushed data handler to parse multiple tags
     * within the current scope.  When the reference scope is about to
     * be exited, control is returned to the previous data handler just
     * before the push.
     */

    if (scopeDepthStack.empty() == false) {
      lastScopeStacked = scopeDepthStack.top();
    }

    if (flagVerbose > 0) {
      printCallInfo("endElement()");
      cout << "Tag name = " << qName << endl;
      cout << "Parser mode = " << getCurrentParseModeStr() << endl;
      cout << "Current scope = " << scopeDepthCount << endl;
      cout << "Last scope stacked = " << lastScopeStacked << endl;
      cout << "Parser mode = PARSE_MODE_NEXT_TAG" << endl;
      cout << "Current data handler name = " << getCurrentDataHandler()->DataHandlerName << endl;
      cout << "Current data handler type = " << getCurrentDataHandler()->DataHandlerType << endl;
      cout << "Calling data handler XML_endElement()" << endl;
    }

    /* check if we are about to exit the scope just before the push of dataHandler */
    if (scopeDepthCount <= lastScopeStacked) { /* scope has returned from child processing */
      popDataHandler(); /* pop the handler back to the original one before the push */
    }

    /* use the data handler before the push to process the end tag */
    getCurrentDataHandler()->XML_endElement(qName, this);

    break; /* end PARSE_MODE_CURRENT_SCOPE */

  } /* end parseMode switch */

  scopeDepthCount--; /* just processed tag so decrease scope by one */

}


/* ================================ SAX Event Handlers ============================= */
/* ================================================================================= */


/* ================================================================================= */
/* ================== Atz_XML_SAX_DataHandler Interface        ===================== */
void Atz_XML_SAX_Handler_Multilevel::XML_startDocument(Atz_XML_SAX_DataHandler* sourceHandler) {
  startDocument();
}


void Atz_XML_SAX_Handler_Multilevel::XML_endDocument(Atz_XML_SAX_DataHandler* sourceHandler) {
  endDocument();
}


void Atz_XML_SAX_Handler_Multilevel::XML_startElement(string qName, AttributesType *attributes, Atz_XML_SAX_DataHandler* sourceHandler) {
  startElement(qName, attributes);
}


void Atz_XML_SAX_Handler_Multilevel::XML_endElement(string qName, Atz_XML_SAX_DataHandler* sourceHandler) {
  endElement(qName);
}


void Atz_XML_SAX_Handler_Multilevel::XML_characters(string xmlString, Atz_XML_SAX_DataHandler* sourceHandler) {
  characters(xmlString);
}

void *Atz_XML_SAX_Handler_Multilevel::XML_getData() {
  return this;
}

void Atz_XML_SAX_Handler_Multilevel::printCallInfo(const char *callName) {
  cout << this->DataHandlerType << " : " << callName << endl;
  cout << "this->DataHandlerName = " << this->DataHandlerName << endl;
}


/* ================== Atz_XML_SAX_DataHandlerInterface wrapper ===================== */
/* ================================================================================= */
