/* ----------------------------------------------------------------------
 XML Handling.

 Paul J. Atzberger
 http://atzberger.org/
------------------------------------------------------------------------- */

#ifndef SELM_LAGRANGIAN_LAMMPS_ATOM_ANGLE_STYLE_XML_HANDLER_H
#define SELM_LAGRANGIAN_LAMMPS_ATOM_ANGLE_STYLE_XML_HANDLER_H

#include "Atz_XML_Package.h"
#include "SELM_Lagrangian_Delegator_XML_Handler.h"
#include "SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE.h"

namespace LAMMPS_NS {


class SELM_Lagrangian_Delegator_XML_Handler; /* declare forward reference of class, since this
                                                also refers to the current class below */

class SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE_XML_Handler : public Atz_XML_SAX_DataHandler {

 public:

  /* ======================== Function prototypes ======================= */
  SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE_XML_Handler();
  SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE_XML_Handler(SELM_Lagrangian_Delegator_XML_Handler *delegatorHandler);

  virtual ~SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE_XML_Handler();

  static const char   *error_str_code;

  string         xmlTagName_xml;
  string         xmlTagName_SELM_Lagrangian;
  string         xmlTagName_LagrangianName;
  string         xmlTagName_num_dim;
  string         xmlTagName_numControlPts;
  string         xmlTagName_ptsX;
  string         xmlTagName_atomID;
  string         xmlTagName_moleculeID;
  string         xmlTagName_typeID;
  string         xmlTagName_atomMass;
  string         xmlTagName_pt_Vel;
  string         xmlTagName_pt_Energy;
  string         xmlTagName_pt_Force;
  string         xmlTagName_pt_type;
  string         xmlTagName_pt_type_extras;

  string         xmlTagName_flagWriteVTK;
  string         xmlTagName_flagWriteSimulationData;
  string         xmlTagName_saveSkipSimulationData;
  string         xmlTagName_outputSimulationData;

  friend class SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE;

  SELM_Lagrangian_LAMMPS_ATOM_ANGLE_STYLE *lagrangian; /* object to construct */

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
