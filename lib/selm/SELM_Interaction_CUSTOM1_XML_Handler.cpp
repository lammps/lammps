/*
 * SELM_Interaction_CUSTOM1_XML_Handler.cpp
 *
 * Paul J. Atzberger
 * http://atzberger.org/
 *
 */

#include "SELM_Interaction_CUSTOM1_XML_Handler.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

using namespace LAMMPS_NS;

SELM_Interaction_CUSTOM1_XML_Handler::SELM_Interaction_CUSTOM1_XML_Handler() {

  setupDataHandler();

  interaction = NULL;

}

SELM_Interaction_CUSTOM1_XML_Handler::SELM_Interaction_CUSTOM1_XML_Handler(SELM_Interaction_Delegator_XML_Handler *delegatorHandler) {

  setupDataHandler();

  interaction = new SELM_Interaction_CUSTOM1();

  strcpy(interaction->nameStr, delegatorHandler->SELM_InteractionName);
  strcpy(interaction->typeStr, delegatorHandler->SELM_InteractionTypeStr);

}

void SELM_Interaction_CUSTOM1_XML_Handler::setupDataHandler() {

  DataHandlerName                    = "Data Handler for SELM_Interaction_CUSTOM1_XML_Handler";
  DataHandlerType                    = "SELM_Interaction_CUSTOM1_XML_Handler";

  xmlTagName_xml                     = "xml";
  xmlTagName_InteractionName         = "InteractionName";
  xmlTagName_InteractionTypeStr      = "InteractionTypeStr";
  xmlTagName_SELM_Interaction        = "SELM_Interaction";
  xmlTagName_numMembers              = "numMembers";
  xmlTagName_memberList_lagrangianI1 = "memberList_lagrangianI1";
  xmlTagName_memberList_ptI1         = "memberList_ptI1";
  xmlTagName_parameterDataList       = "parameterDataList";

}


SELM_Interaction_CUSTOM1_XML_Handler::~SELM_Interaction_CUSTOM1_XML_Handler() {
  // TODO Auto-generated destructor stub
}

void SELM_Interaction_CUSTOM1_XML_Handler::XML_startDocument(Atz_XML_SAX_DataHandler* sourceHandler) {

}

void SELM_Interaction_CUSTOM1_XML_Handler::XML_endDocument(Atz_XML_SAX_DataHandler* sourceHandler) {

}

void SELM_Interaction_CUSTOM1_XML_Handler::XML_startElement(string qName, Atz_XML::AttributesType *attributes, Atz_XML_SAX_DataHandler* sourceHandler) {

  xmlAttributes = attributes;

  xmlString.clear();

  if (qName == xmlTagName_xml) {

  } else if (qName == xmlTagName_SELM_Interaction) {

    /* setup the data structure */
    interaction = new SELM_Interaction_CUSTOM1();

  } else if (qName == xmlTagName_InteractionName) {

  } else if (qName == xmlTagName_InteractionTypeStr) {

  } else if (qName == xmlTagName_numMembers) {

  } else if (qName == xmlTagName_memberList_lagrangianI1) {

  } else if (qName == xmlTagName_memberList_ptI1) {

  } else if (qName == xmlTagName_parameterDataList) {

  } else { /* unrecognized tags skip (avoid sub-tags triggering something here) */
    Atz_XML_SAX_Handler_Multilevel* sourceHandler_Multilevel;
    Atz_XML_SAX_DataHandler* dataHandler;

    sourceHandler_Multilevel = dynamic_cast<Atz_XML_SAX_Handler_Multilevel*> (sourceHandler);
    dataHandler = new Atz_XML_Helper_Handler_SkipNextTag();
    sourceHandler_Multilevel->parseNextTagWithDataHandler(dataHandler);

  }

}

void SELM_Interaction_CUSTOM1_XML_Handler::XML_characters(string xmlString_in, Atz_XML_SAX_DataHandler* sourceHandler) {
  xmlString += xmlString_in;
}

void SELM_Interaction_CUSTOM1_XML_Handler::XML_endElement(string qName, Atz_XML_SAX_DataHandler* sourceHandler) {

  if (qName == xmlTagName_InteractionName) {

    strcpy(interaction->nameStr, Atz_XML_Helper_ParseData::getStringFromAttr(xmlAttributes)->c_str());

  } else if (qName == xmlTagName_InteractionTypeStr) {

  } else if (qName == xmlTagName_numMembers) {

    interaction->numMembers = Atz_XML_Helper_ParseData::getIntFromAttr(xmlAttributes);

  } else if (qName == xmlTagName_memberList_lagrangianI1) {

  } else if (qName == xmlTagName_memberList_ptI1) {

  } else if (qName == xmlTagName_parameterDataList) {

  } else { /* unrecognized tags skip (avoid sub-tags triggering something here) */

  }

}

void *SELM_Interaction_CUSTOM1_XML_Handler::XML_getData() { /* gets data from parsing the XML */
  return interaction;
}
