/* ----------------------------------------------------------------------
 * Coupling operators.
 * 
 * Paul J. Atzberger
 * http://atzberger.org/
 *
------------------------------------------------------------------------- */

#ifndef SELM_COUPLINGOPERATOR_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1_XML_HANDLER_H
#define SELM_COUPLINGOPERATOR_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1_XML_HANDLER_H

#include "Atz_XML_Package.h"
#include "SELM_CouplingOperator_Delegator_XML_Handler.h"
#include "SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1_XML_Handler.h"
#include "SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1.h"
#include "SELM_Eulerian.h"
#include "SELM_Lagrangian.h"

namespace LAMMPS_NS {

class SELM_CouplingOperator_Delegator_XML_Handler; /* declare forward reference of class, since this
                                                       also refers to the current class below */

class SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1_XML_Handler : public Atz_XML_SAX_DataHandler {

 public:

  /* ======================== Function prototypes ======================= */
  SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1_XML_Handler();
  SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1_XML_Handler(SELM_CouplingOperator_Delegator_XML_Handler *delegatorHandler);
  SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1_XML_Handler(SELM_CouplingOperator_Delegator_XML_Handler *delegatorHandler,
                                           int               numLagrangianList,
                                           SELM_Lagrangian **lagrangianList_in,
                                           int               numEulerianList,
                                           SELM_Eulerian   **eulerianList_in);

  virtual ~SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1_XML_Handler();

  string         xmlTagName_xml;
  string         xmlTagName_SELM_CouplingOperator;
  string         xmlTagName_CouplingOperatorName;
  string         xmlTagName_CouplingOperatorTypeStr;
  string         xmlTagName_operatorData;

  string         xmlTagName_numCoupleList;

  string         xmlTagName_lagrangianList;

  string         xmlTagName_SELM_Lagrangian_Ref;

  string         xmlTagName_LagrangianName;
  string         xmlTagName_LagrangianTypeStr;

  string         xmlTagName_eulerianList;

  string         xmlTagName_SELM_Eulerian_Ref;

  string         xmlTagName_EulerianName;
  string         xmlTagName_EulerianTypeStr;

  string         xmlTagName_T_KERNEL_1_weightTableFilename;
  string         xmlTagName_flagWriteSimulationData;
  string         xmlTagName_saveSkipSimulationData;

  friend class SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1;

  SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1 *couplingOp; /* object to construct */
  int                           flagDeterminedCouplingOpType;

  Atz_XML::AttributesType *xmlAttributes;
  string xmlString;

  int               numEulerianList;
  SELM_Eulerian   **eulerianList;

  int               numLagrangianList;
  SELM_Lagrangian **lagrangianList;

  char            **lagrangianNameList;
  int               lagrangianNameList_I;

  char            **lagrangianTypeStrList;
  int               lagrangianTypeStrList_I;

  char            **eulerianNameList;
  int               eulerianNameList_I;

  char            **eulerianTypeStrList;
  int               eulerianTypeStrList_I;

  char basePath[10000];

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
