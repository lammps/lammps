/*
 * This class delegates the parsing of XML files to the appropriate
 * SELM_Interaction class associated with the given data type.  This
 * allows for objects to be instantiated from XML files for the
 * specified types.
 *
 * Paul J. Atzberger
 * http://atzberger.org/
 *
 */

#ifndef SELM_INTERACTION_DELEGATOR_XML_HANDLER_H
#define SELM_INTERACTION_DELEGATOR_XML_HANDLER_H

#include "Atz_XML_Package.h"
#include "SELM_Interaction.h"
#include "SELM_Interaction_CUSTOM1_XML_Handler.h"
#include "SELM_Interaction_SKIPDATA_XML_Handler.h"

namespace LAMMPS_NS {

class SELM_Interaction_Delegator_XML_Handler : public Atz_XML_SAX_DataHandler {

 public:

  /* ======================== Constants ======================= */
  static const int PARSE_MODE_NULL           = 0;
  static const int PARSE_MODE_HANDLE_LOCALLY = 1;
  static const int PARSE_MODE_DELEGATE       = 2;

  /* ======================== Function prototypes ======================= */
  SELM_Interaction_Delegator_XML_Handler();
  virtual ~SELM_Interaction_Delegator_XML_Handler();

  string         xmlTagName_xml;
  string         xmlTagName_SELM_Interaction;
  string         xmlTagName_InteractionName;
  string         xmlTagName_InteractionTypeStr;

  char           SELM_InteractionName[1000];
  char           SELM_InteractionTypeStr[1000];

  Atz_XML::AttributesType *xmlAttributes;
  string xmlString;

  int    parseMode;

  Atz_XML_SAX_DataHandler*  delegatee_dataHandler;

 public:
   void XML_startDocument(Atz_XML_SAX_DataHandler* sourceHandler);

   void XML_endDocument(Atz_XML_SAX_DataHandler* sourceHandler);

   void XML_startElement(string qName, Atz_XML::AttributesType *attributes, Atz_XML_SAX_DataHandler* sourceHandler);

   void XML_characters(string xmlString, Atz_XML_SAX_DataHandler* sourceHandler);

   void XML_endElement(string qName, Atz_XML_SAX_DataHandler* sourceHandler);

   void *XML_getData(); /* gets data from parsing the XML */

};

}

#endif
