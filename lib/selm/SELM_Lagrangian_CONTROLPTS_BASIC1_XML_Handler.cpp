/*

 SELM_Lagrangian_CONTROLPTS_BASIC1_XML_Handler.cpp
 
 Paul J. Atzberger
 http://atzberger.org/
 
*/

#include "SELM_Lagrangian_CONTROLPTS_BASIC1_XML_Handler.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

using namespace LAMMPS_NS;

SELM_Lagrangian_CONTROLPTS_BASIC1_XML_Handler::SELM_Lagrangian_CONTROLPTS_BASIC1_XML_Handler() {

  setupDataHandler();

  lagrangian = NULL;

}

SELM_Lagrangian_CONTROLPTS_BASIC1_XML_Handler::SELM_Lagrangian_CONTROLPTS_BASIC1_XML_Handler(SELM_Lagrangian_Delegator_XML_Handler *delegatorHandler) {

  setupDataHandler();

  lagrangian = new SELM_Lagrangian_CONTROLPTS_BASIC1();

  strcpy(lagrangian->nameStr, delegatorHandler->SELM_LagrangianName);
  strcpy(lagrangian->typeStr, delegatorHandler->SELM_LagrangianTypeStr);

}

void SELM_Lagrangian_CONTROLPTS_BASIC1_XML_Handler::setupDataHandler() {

  DataHandlerName            = "Data Handler for SELM_Lagrangian_CONTROLPTS_BASIC1";
  DataHandlerType            = "SELM_Lagrangian_CONTROLPTS_BASIC1_XML_Handler";

  xmlTagName_xml             = "xml";
  xmlTagName_LagrangianName  = "LagrangianName";
  xmlTagName_SELM_Lagrangian = "SELM_Lagrangian";
  xmlTagName_num_dim         = "num_dim";
  xmlTagName_numControlPts   = "numControlPts";
  xmlTagName_pt_X            = "pt_X";
  xmlTagName_pt_Vel          = "pt_Vel";
  xmlTagName_pt_Energy       = "pt_Energy";
  xmlTagName_pt_Force        = "pt_Force";
  xmlTagName_pt_type         = "pt_type";
  xmlTagName_pt_type_extras  = "pt_type_extras";

}


SELM_Lagrangian_CONTROLPTS_BASIC1_XML_Handler::~SELM_Lagrangian_CONTROLPTS_BASIC1_XML_Handler() {
  // TODO Auto-generated destructor stub
}

void SELM_Lagrangian_CONTROLPTS_BASIC1_XML_Handler::XML_startDocument(Atz_XML_SAX_DataHandler* sourceHandler) {

}

void SELM_Lagrangian_CONTROLPTS_BASIC1_XML_Handler::XML_endDocument(Atz_XML_SAX_DataHandler* sourceHandler) {

}

void SELM_Lagrangian_CONTROLPTS_BASIC1_XML_Handler::XML_startElement(string qName, Atz_XML::AttributesType *attributes, Atz_XML_SAX_DataHandler* sourceHandler) {

  xmlAttributes = attributes;

  xmlString.clear();

  if (qName == xmlTagName_xml) {

  } else if (qName == xmlTagName_SELM_Lagrangian) {

    /* setup the data structure */
    lagrangian = new SELM_Lagrangian_CONTROLPTS_BASIC1();

  } else if (qName == xmlTagName_LagrangianName) {

  } else if (qName == xmlTagName_num_dim) {

  } else if (qName == xmlTagName_numControlPts) {

  } else if (qName == xmlTagName_pt_X) {

  } else if (qName == xmlTagName_pt_Vel) {

  } else if (qName == xmlTagName_pt_type) {

  } else if (qName == xmlTagName_pt_type_extras) {

  } else { /* unrecognized tags skip (avoid sub-tags triggering something here) */
    Atz_XML_SAX_Handler_Multilevel*    sourceHandler_Multilevel;
    Atz_XML_SAX_DataHandler*           dataHandler;

    sourceHandler_Multilevel = dynamic_cast<Atz_XML_SAX_Handler_Multilevel*>(sourceHandler);
    dataHandler              = new Atz_XML_Helper_Handler_SkipNextTag();
    sourceHandler_Multilevel->parseNextTagWithDataHandler(dataHandler);
  }

}

void SELM_Lagrangian_CONTROLPTS_BASIC1_XML_Handler::XML_characters(string xmlString_in, Atz_XML_SAX_DataHandler* sourceHandler) {
  xmlString += xmlString_in;
}

void SELM_Lagrangian_CONTROLPTS_BASIC1_XML_Handler::XML_endElement(string qName, Atz_XML_SAX_DataHandler* sourceHandler) {

  if (qName == xmlTagName_LagrangianName) {

    strcpy(lagrangian->nameStr, Atz_XML_Helper_ParseData::getStringFromAttr(xmlAttributes)->c_str());

  } else if (qName == xmlTagName_num_dim) {

    lagrangian->num_dim = Atz_XML_Helper_ParseData::getIntFromAttr(xmlAttributes);

  } else if (qName == xmlTagName_numControlPts) {

    lagrangian->numControlPts       = Atz_XML_Helper_ParseData::getIntFromAttr(xmlAttributes);
    //lagrangian->numControlPts_alloc = lagrangian->numControlPts;

  } else if (qName == xmlTagName_pt_X) {

    int N;
    lagrangian->pt_X = NULL; /* indicates parse routine should allocate array */
    Atz_XML_Helper_ParseData::parseDoubleArrayFromString(xmlString, &lagrangian->pt_X, &N);

  } else if (qName == xmlTagName_pt_Vel) {

    int N;
    lagrangian->pt_Vel = NULL; /* indicates parse routine should allocate array */
    Atz_XML_Helper_ParseData::parseDoubleArrayFromString(xmlString, &lagrangian->pt_Vel, &N);

  } else if (qName == xmlTagName_pt_type) {

    /* parse list from string of characters */

  } else if (qName == xmlTagName_pt_type_extras) {

    /* type specific data... */

  } else { /* unrecognized tags skip (avoid sub-tags triggering something here) */

  }

}

void *SELM_Lagrangian_CONTROLPTS_BASIC1_XML_Handler::XML_getData() { /* gets data from parsing the XML */
  return lagrangian;
}
