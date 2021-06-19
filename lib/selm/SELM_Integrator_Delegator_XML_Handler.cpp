/*

  SELM_Integrator_Delegator_XML_Handler.cpp
 
  Paul J. Atzberger
  http://atzberger.org/
 
*/

#include "SELM_Integrator_Delegator_XML_Handler.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

using namespace LAMMPS_NS;

SELM_Integrator_Delegator_XML_Handler::SELM_Integrator_Delegator_XML_Handler() {

  DataHandlerName               = "Data handler for SELM_Integrator_Delegator_XML_Handler";
  DataHandlerType               = "SELM_Integrator_Delegator_XML_Handler";

  xmlTagName_xml                = "xml";
  xmlTagName_SELM_Integrator    = "SELM_Integrator";
  xmlTagName_IntegratorName     = "IntegratorName";
  xmlTagName_IntegratorTypeStr  = "IntegratorTypeStr";

  delegatee_dataHandler         = NULL;

  parseMode                     = PARSE_MODE_HANDLE_LOCALLY;

}

SELM_Integrator_Delegator_XML_Handler::~SELM_Integrator_Delegator_XML_Handler() {
  // TODO Auto-generated destructor stub
}

void SELM_Integrator_Delegator_XML_Handler::XML_startDocument(Atz_XML_SAX_DataHandler* sourceHandler) {

}

void SELM_Integrator_Delegator_XML_Handler::XML_endDocument(Atz_XML_SAX_DataHandler* sourceHandler) {

}

void SELM_Integrator_Delegator_XML_Handler::XML_startElement(string qName, Atz_XML::AttributesType *attributes, Atz_XML_SAX_DataHandler* sourceHandler) {

  switch (parseMode) {

  case PARSE_MODE_HANDLE_LOCALLY:

    xmlAttributes = attributes;
    xmlString.clear();

    if (qName == xmlTagName_xml) {

    } else if (qName == xmlTagName_SELM_Integrator) {

    } else if (qName == xmlTagName_IntegratorName) {

    } else if (qName == xmlTagName_IntegratorTypeStr) {

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

void SELM_Integrator_Delegator_XML_Handler::XML_characters(string xmlString_in, Atz_XML_SAX_DataHandler* sourceHandler) {

  switch (parseMode) {

  case PARSE_MODE_HANDLE_LOCALLY:
    xmlString += xmlString_in;
    break;

  case PARSE_MODE_DELEGATE:
    delegatee_dataHandler->XML_characters(xmlString_in,sourceHandler);
    break;

  } /* end switch */

}

void SELM_Integrator_Delegator_XML_Handler::XML_endElement(string qName, Atz_XML_SAX_DataHandler* sourceHandler) {

  const char *error_str_code = "SELM_Integrator_Delegator.cpp";
  const char *error_str_func = "XML_endElement()";

  switch (parseMode) {

  case PARSE_MODE_HANDLE_LOCALLY:

    if (qName == xmlTagName_IntegratorName) {

      strcpy(SELM_IntegratorName, Atz_XML_Helper_ParseData::getStringFromAttr(xmlAttributes)->c_str());

    } else if (qName == xmlTagName_IntegratorTypeStr) {

      strcpy(SELM_IntegratorTypeStr, Atz_XML_Helper_ParseData::getStringFromAttr(xmlAttributes)->c_str());

      /* Determine the type and delegate the rest of the parsing
       * to the appropriate XML data handler.
       */
      if (strcmp(SELM_IntegratorTypeStr, SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::TYPE_STR) == 0) {
        delegatee_dataHandler = new SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_XML_Handler(this);
      } else {
        stringstream message;
        message << "Integrator type was not recognized" << endl;
        message << "SELM_IntegratorTypeStr = " << SELM_IntegratorTypeStr << endl;
        SELM_Package::packageError(error_str_code, error_str_func, message);
      }

      parseMode = PARSE_MODE_DELEGATE;

    } else { /* unrecognized tags skip (avoid sub-tags triggering something here) */

    }

    break;

  case PARSE_MODE_DELEGATE:

    delegatee_dataHandler->XML_endElement(qName, sourceHandler);

    if (qName == xmlTagName_SELM_Integrator) {
      parseMode = PARSE_MODE_HANDLE_LOCALLY;
    }

    break;

  } /* end switch */

}

void *SELM_Integrator_Delegator_XML_Handler::XML_getData() { /* gets data from parsing the XML */

  void *XML_Data_ptr = NULL;

  if (delegatee_dataHandler != NULL) {
    XML_Data_ptr = delegatee_dataHandler->XML_getData();
  }

  return XML_Data_ptr;

}
