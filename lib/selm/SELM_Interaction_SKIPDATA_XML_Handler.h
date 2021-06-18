/* ----------------------------------------------------------------------
 Custom interaction.

 Paul J. Atzberger
 http://atzberger.org/
 
------------------------------------------------------------------------- */

#ifndef SELM_INTERACTION_SKIPDATA_XML_HANDLER_H
#define SELM_INTERACTION_SKIPDATA_XML_HANDLER_H

#include "Atz_XML_Package.h"
#include "SELM_Interaction_Delegator_XML_Handler.h"
#include "SELM_Interaction_SKIPDATA.h"

namespace LAMMPS_NS {

class SELM_Interaction_Delegator_XML_Handler; /* declare forward reference of class, since this
                                                also refers to the current class below */

class SELM_Interaction_SKIPDATA_XML_Handler : public Atz_XML_SAX_DataHandler {

 public:

  /* ======================== Function prototypes ======================= */
  SELM_Interaction_SKIPDATA_XML_Handler();
  SELM_Interaction_SKIPDATA_XML_Handler(SELM_Interaction_Delegator_XML_Handler *delegatorHandler);

  virtual ~SELM_Interaction_SKIPDATA_XML_Handler();

  string         xmlTagName_xml;
  string         xmlTagName_SELM_Interaction;
  string         xmlTagName_InteractionName;
  string         xmlTagName_InteractionTypeStr;
  string         xmlTagName_numMembers;
  string         xmlTagName_memberList_lagrangianI1;
  string         xmlTagName_memberList_ptI1;
  string         xmlTagName_parameterDataList;

  friend class SELM_Interaction_SKIPDATA;

  SELM_Interaction_SKIPDATA *interaction; /* object to construct */

  Atz_XML::AttributesType *xmlAttributes;
  string xmlString;

 public:

   void setupDataHandler();

   void XML_startDocument(Atz_XML_SAX_DataHandler* sourceHandler);

   void XML_endDocument(Atz_XML_SAX_DataHandler* sourceHandler);

   void XML_startElement(string qName, Atz_XML::AttributesType *attributes, Atz_XML_SAX_DataHandler* sourceHandler);

   void XML_characters(string xmlString, Atz_XML_SAX_DataHandler* sourceHandler);

   void XML_endElement(string qName, Atz_XML_SAX_DataHandler* sourceHandler);

   void *XML_getData(); /* gets data from parsing the XML */

};

}

#endif
