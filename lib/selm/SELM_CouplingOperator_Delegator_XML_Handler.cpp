/*
 SELM_CouplingOperators_Delegator_XML_Handler.cpp
 
 Paul J. Atzberger
 http://atzberger.org/
 
*/

#include "SELM_CouplingOperator_Delegator_XML_Handler.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

using namespace LAMMPS_NS;

SELM_CouplingOperator_Delegator_XML_Handler::SELM_CouplingOperator_Delegator_XML_Handler() {
  setup();
}

SELM_CouplingOperator_Delegator_XML_Handler::SELM_CouplingOperator_Delegator_XML_Handler(int               numLagrangianList_in,
                                                                                         SELM_Lagrangian **lagrangianList_in,
                                                                                         int               numEulerianList_in,
                                                                                         SELM_Eulerian   **eulerianList_in) {

  string numLagrangianList_str("numLagrangianList");
  string lagrangianList_str("lagrangianList");

  string numEulerianList_str("numEulerianList");
  string eulerianList_str("eulerianList");

  int    *numLagrangianList_ptr = (int *)malloc(sizeof(int));
  int    *numEulerianList_ptr   = (int *)malloc(sizeof(int));

  numLagrangianList_ptr[0] = numLagrangianList_in;
  numEulerianList_ptr[0]   = numEulerianList_in;

  setup();

  /* reference lists */
  extraData = new ExtraDataType();

  extraData->insert(ExtraDataPairType(numLagrangianList_str, (void *)numLagrangianList_ptr));
  extraData->insert(ExtraDataPairType(lagrangianList_str, (void *)lagrangianList_in));

  extraData->insert(ExtraDataPairType(numEulerianList_str, (void *)numEulerianList_ptr));
  extraData->insert(ExtraDataPairType(eulerianList_str, (void *)eulerianList_in));

}


void SELM_CouplingOperator_Delegator_XML_Handler::setup() {

  DataHandlerName                     = "Data handler for SELM_CouplingOperator_Delegator_XML_Handler";
  DataHandlerType                     = "SELM_CouplingOperator_Delegator_XML_Handler";

  xmlTagName_xml                      = "xml";
  xmlTagName_SELM_CouplingOperator    = "SELM_CouplingOperator";
  xmlTagName_CouplingOperatorName     = "CouplingOperatorName";
  xmlTagName_CouplingOperatorTypeStr  = "CouplingOperatorTypeStr";

  delegatee_dataHandler               = NULL;

  parseMode                           = PARSE_MODE_HANDLE_LOCALLY;

  extraData                           = NULL;

}

SELM_CouplingOperator_Delegator_XML_Handler::~SELM_CouplingOperator_Delegator_XML_Handler() {
  // TODO Auto-generated destructor stub
}

void SELM_CouplingOperator_Delegator_XML_Handler::XML_startDocument(Atz_XML_SAX_DataHandler* sourceHandler) {

}

void SELM_CouplingOperator_Delegator_XML_Handler::XML_endDocument(Atz_XML_SAX_DataHandler* sourceHandler) {

}

void SELM_CouplingOperator_Delegator_XML_Handler::XML_startElement(string qName, Atz_XML::AttributesType *attributes, Atz_XML_SAX_DataHandler* sourceHandler) {

  switch (parseMode) {

  case PARSE_MODE_HANDLE_LOCALLY:

    xmlAttributes = attributes;
    xmlString.clear();

    if (qName == xmlTagName_xml) {

    } else if (qName == xmlTagName_SELM_CouplingOperator) {

      /* setup the data structure */
      //lagrangian = new SELM_CouplingOperators_TABLE1();

    } else if (qName == xmlTagName_CouplingOperatorName) {

    } else if (qName == xmlTagName_CouplingOperatorTypeStr) {

    } else { /* unrecognized tags skip (avoid sub-tags triggering something here) */
      Atz_XML_SAX_Handler_Multilevel*    sourceHandler_Multilevel;
      Atz_XML_SAX_DataHandler*           dataHandler;

      sourceHandler_Multilevel = dynamic_cast<Atz_XML_SAX_Handler_Multilevel*>(sourceHandler);
      dataHandler              = new Atz_XML_Helper_Handler_SkipNextTag();
      sourceHandler_Multilevel->parseNextTagWithDataHandler(dataHandler);
    }

    break;

  case PARSE_MODE_DELEGATE:

    delegatee_dataHandler->XML_startElement(qName, attributes, sourceHandler);

    break;

  } /* end of switch */

}

void SELM_CouplingOperator_Delegator_XML_Handler::XML_characters(string xmlString_in, Atz_XML_SAX_DataHandler* sourceHandler) {

  switch (parseMode) {

  case PARSE_MODE_HANDLE_LOCALLY:
    xmlString += xmlString_in;
    break;

  case PARSE_MODE_DELEGATE:
    delegatee_dataHandler->XML_characters(xmlString_in,sourceHandler);
    break;

  } /* end switch */

}

void SELM_CouplingOperator_Delegator_XML_Handler::XML_endElement(string qName, Atz_XML_SAX_DataHandler* sourceHandler) {

  const char *error_str_code = "SELM_CouplingOperators_Delegator.cpp";
  const char *error_str_func = "XML_endElement()";

  switch (parseMode) {

  case PARSE_MODE_HANDLE_LOCALLY:

    if (qName == xmlTagName_CouplingOperatorName) {

      strcpy(SELM_CouplingOperatorName, Atz_XML_Helper_ParseData::getStringFromAttr(xmlAttributes)->c_str());

    } else if (qName == xmlTagName_CouplingOperatorTypeStr) {

      strcpy(SELM_CouplingOperatorTypeStr, Atz_XML_Helper_ParseData::getStringFromAttr(xmlAttributes)->c_str());

      /* Determine the type and delegate the rest of the parsing
       * to the appropriate XML data handler.
       */
      if (strcmp(SELM_CouplingOperatorTypeStr, SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::TYPE_STR) == 0) {
        delegatee_dataHandler = new SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1_XML_Handler(this,
                         *((int *)(SELM_Lagrangian **)extraData->find("numLagrangianList")->second),
                         (SELM_Lagrangian **)extraData->find("lagrangianList")->second,
                         *((int *)((SELM_Eulerian **)extraData->find("numEulerianList")->second)),
                         (SELM_Eulerian **)extraData->find("eulerianList")->second);
      } else {
        stringstream message;
        message << "CouplingOperators type was not recognized" << endl;
        message << "SELM_CouplingOperatorsTypeStr = " << SELM_CouplingOperatorTypeStr << endl;
        SELM_Package::packageError(error_str_code, error_str_func, message);
      }

      parseMode = PARSE_MODE_DELEGATE;

    } else { /* unrecognized tags skip (avoid sub-tags triggering something here) */

    }

    break;

  case PARSE_MODE_DELEGATE:

    delegatee_dataHandler->XML_endElement(qName, sourceHandler);

    if (qName == xmlTagName_SELM_CouplingOperator) {
      parseMode = PARSE_MODE_HANDLE_LOCALLY;
    }

    break;

  } /* end switch */

}

void *SELM_CouplingOperator_Delegator_XML_Handler::XML_getData() { /* gets data from parsing the XML */

  void *XML_Data_ptr = NULL;

  if (delegatee_dataHandler != NULL) {
    XML_Data_ptr = delegatee_dataHandler->XML_getData();
  }

  return XML_Data_ptr;

}
