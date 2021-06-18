/* ----------------------------------------------------------------------
 Handle XML data.

 Paul J. Atzberger
 http://atzberger.org/
 
------------------------------------------------------------------------- */

#ifndef FIX_SELM_XML_HANDLER_H
#define FIX_SELM_XML_HANDLER_H

#include "driver_selm.h"
#include "Atz_XML_Package.h"
#include "SELM_Lagrangian_Delegator_XML_Handler.h"
#include "SELM_Eulerian_Delegator_XML_Handler.h"
#include "SELM_CouplingOperator_Delegator_XML_Handler.h"
#include "SELM_Integrator_Delegator_XML_Handler.h"
#include "SELM_Interaction_Delegator_XML_Handler.h"

namespace LAMMPS_NS {

class Driver_SELM_XML_Handler : public Atz_XML_SAX_DataHandler {

 public:

  /* ======================== Function prototypes ======================= */
  Driver_SELM_XML_Handler(DriverSELM *Driver_SELM_Data_ptr);

  virtual ~Driver_SELM_XML_Handler();

  int              parseMode;

  static const int TAG_NOT_RECOGNIZED               = 0;
  static const int TAG_RECOGNIZED                   = 1;

  static const int PARSE_MODE_NULL                        = 0;
  static const int PARSE_MODE_FIX_SELM                    = 1;
  static const int PARSE_MODE_SELM_Lagrangian_List        = 2;
  static const int PARSE_MODE_SELM_Eulerian_List          = 3;
  static const int PARSE_MODE_SELM_CouplingOperator_List  = 4;
  static const int PARSE_MODE_SELM_Integrator             = 5;
  static const int PARSE_MODE_SELM_Interaction_List       = 6;


  typedef struct SELM_Lagrangian_ParamsType {
    char SELM_LagrangianName[1000];
    char SELM_LagrangianTypeStr[1000];
  } SELM_Lagrangian_ParamsType;

  typedef struct SELM_Eulerian_ParamsType {
    char SELM_EulerianName[1000];
    char SELM_EulerianTypeStr[1000];
  } SELM_Eulerian_ParamsType;

  typedef struct SELM_CouplingOperator_ParamsType {
    char SELM_CouplingOperatorName[1000];
    char SELM_CouplingOperatorTypeStr[1000];
  } SELM_CouplingOperator_ParamsType;

  typedef struct SELM_Interaction_ParamsType {
    char SELM_InteractionName[1000];
    char SELM_InteractionTypeStr[1000];
  } SELM_Interaction_ParamsType;

  typedef struct SELM_Integrator_ParamsType {
    char SELM_IntegratorName[1000];
    char SELM_IntegratorTypeStr[1000];
  } SELM_Integrator_ParamsType;

  string         xmlTagName_xml;
  string         xmlTagName_FixSELM;

  string         xmlTagName_SELM_Version;
  string         xmlTagName_SELM_Run_Description;
  string         xmlTagName_SELM_BasePath;
  string         xmlTagName_SELM_BaseFilename;
  string         xmlTagName_SELM_Seed;

  string         xmlTagName_SELM_Lagrangian_List;
  string         xmlTagName_SELM_Lagrangian;
  string         xmlTagName_SELM_LagrangianName;
  string         xmlTagName_SELM_LagrangianTypeStr;

  string         xmlTagName_SELM_Eulerian_List;
  string         xmlTagName_SELM_Eulerian;
  string         xmlTagName_SELM_EulerianName;
  string         xmlTagName_SELM_EulerianTypeStr;

  string         xmlTagName_SELM_CouplingOperator_List;
  string         xmlTagName_SELM_CouplingOperator;
  string         xmlTagName_SELM_CouplingOperatorName;
  string         xmlTagName_SELM_CouplingOperatorTypeStr;

  string         xmlTagName_SELM_Interaction_List;
  string         xmlTagName_SELM_Interaction;
  string         xmlTagName_SELM_InteractionName;
  string         xmlTagName_SELM_InteractionTypeStr;

  string         xmlTagName_SELM_Integrator;
  string         xmlTagName_SELM_IntegratorName;
  string         xmlTagName_SELM_IntegratorTypeStr;

  string         xmlTagName_flagWriteSimulationData;
  string         xmlTagName_saveSkipSimulationData;

  friend class                       DriverSELM;
  DriverSELM                        *driver_SELM_Data; /* object to construct */

  int                                SELM_Lagrangian_List_Params_N;
  int                                SELM_Lagrangian_List_Params_I;
  SELM_Lagrangian_ParamsType       **SELM_Lagrangian_List_Params;

  int                                SELM_Eulerian_List_Params_N;
  int                                SELM_Eulerian_List_Params_I;
  SELM_Eulerian_ParamsType         **SELM_Eulerian_List_Params;

  int                                SELM_CouplingOperator_List_Params_N;
  int                                SELM_CouplingOperator_List_Params_I;
  SELM_CouplingOperator_ParamsType **SELM_CouplingOperator_List_Params;

  int                                SELM_Interaction_List_Params_N;
  int                                SELM_Interaction_List_Params_I;
  SELM_Interaction_ParamsType      **SELM_Interaction_List_Params;

  SELM_Integrator_ParamsType        *SELM_Integrator_Params;

  Atz_XML::AttributesType *xmlAttributes;
  string                   xmlString;


 public:
   void XML_startDocument(Atz_XML_SAX_DataHandler* sourceHandler);

   void XML_endDocument(Atz_XML_SAX_DataHandler* sourceHandler);

   void XML_startElement(string qName, Atz_XML::AttributesType *attributes, Atz_XML_SAX_DataHandler* sourceHandler);

   void XML_characters(string xmlString, Atz_XML_SAX_DataHandler* sourceHandler);

   void XML_endElement(string qName, Atz_XML_SAX_DataHandler* sourceHandler);

   void *XML_getData(); /* gets data from parsing the XML */

   /* setup the DriverSELM from the parameters data */
   void setupDriverSELM_From_Params();

};

}

#endif
