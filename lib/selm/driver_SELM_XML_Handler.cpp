/*
 Driver_SELM_XML_Handler.cpp
 
 Paul J. Atzberger
 http://atzberger.org/
 
*/

#include "driver_SELM_XML_Handler.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <sys/stat.h>

using namespace LAMMPS_NS;

Driver_SELM_XML_Handler::Driver_SELM_XML_Handler(DriverSELM *driver_SELM_Data_ptr) {

  DataHandlerName                         = "Default Data Handler";
  DataHandlerType                         = "Driver_SELM_XML_Handler";

  xmlTagName_xml                          = "xml";
  xmlTagName_FixSELM                      = "FixSELM";

  xmlTagName_SELM_Version                 = "SELM_Version";
  xmlTagName_SELM_Run_Description         = "SELM_Run_Description";

  xmlTagName_SELM_BasePath                = "SELM_BasePath";
  xmlTagName_SELM_BaseFilename            = "SELM_BaseFilename";

  xmlTagName_SELM_Seed                    = "SELM_Seed";

  xmlTagName_SELM_Lagrangian_List         = "SELM_Lagrangian_List";
  xmlTagName_SELM_Lagrangian              = "SELM_Lagrangian";
  xmlTagName_SELM_LagrangianName          = "SELM_LagrangianName";
  xmlTagName_SELM_LagrangianTypeStr       = "SELM_LagrangianTypeStr";

  xmlTagName_SELM_Eulerian_List           = "SELM_Eulerian_List";
  xmlTagName_SELM_Eulerian                = "SELM_Eulerian";
  xmlTagName_SELM_EulerianName            = "SELM_EulerianName";
  xmlTagName_SELM_EulerianTypeStr         = "SELM_EulerianTypeStr";

  xmlTagName_SELM_CouplingOperator_List   = "SELM_CouplingOperator_List";
  xmlTagName_SELM_CouplingOperator        = "SELM_CouplingOperator";
  xmlTagName_SELM_CouplingOperatorName    = "SELM_CouplingOperatorName";
  xmlTagName_SELM_CouplingOperatorTypeStr = "SELM_CouplingOperatorTypeStr";

  xmlTagName_SELM_Interaction_List        = "SELM_Interaction_List";
  xmlTagName_SELM_Interaction             = "SELM_Interaction";
  xmlTagName_SELM_InteractionName         = "SELM_InteractionName";
  xmlTagName_SELM_InteractionTypeStr      = "SELM_InteractionTypeStr";

  xmlTagName_SELM_Integrator              = "SELM_Integrator";
  xmlTagName_SELM_IntegratorName          = "SELM_IntegratorName";
  xmlTagName_SELM_IntegratorTypeStr       = "SELM_IntegratorTypeStr";

  xmlTagName_flagWriteSimulationData      = "flagWriteSimulationData";
  xmlTagName_saveSkipSimulationData       = "saveSkipSimulationData";

  driver_SELM_Data                        = driver_SELM_Data_ptr;
  parseMode                               = PARSE_MODE_FIX_SELM;

}

Driver_SELM_XML_Handler::~Driver_SELM_XML_Handler() {
  // TODO Auto-generated destructor stub
}

void Driver_SELM_XML_Handler::XML_startDocument(Atz_XML_SAX_DataHandler* sourceHandler) {
  parseMode = PARSE_MODE_FIX_SELM;
}

void Driver_SELM_XML_Handler::XML_endDocument(Atz_XML_SAX_DataHandler* sourceHandler) {

}

void Driver_SELM_XML_Handler::XML_startElement(string qName, Atz_XML::AttributesType *attributes, Atz_XML_SAX_DataHandler* sourceHandler) {

  const char *error_str_code = "Driver_SELM_XML_Handler.cpp";
  const char *error_str_func = "XML_startElement()";

  int flagRecognizedTag = TAG_NOT_RECOGNIZED;

  /* setup the parser */
  xmlAttributes = attributes;
  xmlString.clear();

  /* check tags to determine initial parser mode */
  if ((qName == xmlTagName_xml) || (qName == xmlTagName_FixSELM)) {
    flagRecognizedTag = TAG_RECOGNIZED;

    /* setup the data structure */
    //driver_SELM_Data = new DriverSELM();
    parseMode = PARSE_MODE_FIX_SELM;
  }

  /* perform parsing based on mode */
  switch (parseMode) {

  case PARSE_MODE_FIX_SELM:

    if (qName == xmlTagName_xml) {
      flagRecognizedTag = TAG_RECOGNIZED;
    } else if (qName == xmlTagName_FixSELM) {
      flagRecognizedTag = TAG_RECOGNIZED;
    } else if (qName == xmlTagName_SELM_Version) {
      flagRecognizedTag = TAG_RECOGNIZED;
    } else if (qName == xmlTagName_SELM_Run_Description) {
      flagRecognizedTag = TAG_RECOGNIZED;
    } else if (qName == xmlTagName_SELM_BasePath) {
      flagRecognizedTag = TAG_RECOGNIZED;
    } else if (qName == xmlTagName_SELM_BaseFilename) {
      flagRecognizedTag = TAG_RECOGNIZED;
    } else if (qName == xmlTagName_SELM_Seed) {
      flagRecognizedTag = TAG_RECOGNIZED;
    } else if (qName == xmlTagName_SELM_Lagrangian_List) {

      /* indicate tag recognized */
      flagRecognizedTag = TAG_RECOGNIZED;

      /* setup the list data structure for parameter data */
      int N = Atz_XML_Helper_ParseData::getIntFromAttr("numEntries", xmlAttributes);
      SELM_Lagrangian_List_Params_N = N;
      SELM_Lagrangian_List_Params_I = 0;
      SELM_Lagrangian_List_Params
        = (Driver_SELM_XML_Handler::SELM_Lagrangian_ParamsType **)
          malloc(sizeof(Driver_SELM_XML_Handler::SELM_Lagrangian_ParamsType *)*N);

      /* allocate the entries in the list */
      for (int I = 0; I < N; I++) {
        SELM_Lagrangian_List_Params[I]
          = (Driver_SELM_XML_Handler::SELM_Lagrangian_ParamsType *)
            malloc(sizeof(Driver_SELM_XML_Handler::SELM_Lagrangian_ParamsType));
      } /* end of I loop */

      /* parse tags containing list of names + types */
      parseMode = PARSE_MODE_SELM_Lagrangian_List;

    } else if (qName == xmlTagName_SELM_Eulerian_List) {

      /* indicate tag recognized */
      flagRecognizedTag = TAG_RECOGNIZED;

      /* setup the list data structure for parameter data */
      int N = Atz_XML_Helper_ParseData::getIntFromAttr("numEntries", xmlAttributes);
      SELM_Eulerian_List_Params_N = N;
      SELM_Eulerian_List_Params_I = 0;
      SELM_Eulerian_List_Params
        = (Driver_SELM_XML_Handler::SELM_Eulerian_ParamsType **)
          malloc(sizeof(Driver_SELM_XML_Handler::SELM_Eulerian_ParamsType *)*N);

      /* allocate the entries in the list */
      for (int I = 0; I < N; I++) {
        SELM_Eulerian_List_Params[I]
          = (Driver_SELM_XML_Handler::SELM_Eulerian_ParamsType *)
            malloc(sizeof(Driver_SELM_XML_Handler::SELM_Eulerian_ParamsType));
      } /* end of I loop */

      /* parse tags containing list of names + types */
      parseMode = PARSE_MODE_SELM_Eulerian_List;

    } else if (qName == xmlTagName_SELM_CouplingOperator_List) {

      /* indicate tag recognized */
      flagRecognizedTag = TAG_RECOGNIZED;

      /* setup the list data structure for parameter data */
      int N = Atz_XML_Helper_ParseData::getIntFromAttr("numEntries", xmlAttributes);
      SELM_CouplingOperator_List_Params_N = N;
      SELM_CouplingOperator_List_Params_I = 0;
      SELM_CouplingOperator_List_Params
        = (Driver_SELM_XML_Handler::SELM_CouplingOperator_ParamsType **)
          malloc(sizeof(Driver_SELM_XML_Handler::SELM_CouplingOperator_ParamsType *)*N);

      /* allocate the entries in the list */
      for (int I = 0; I < N; I++) {
        SELM_CouplingOperator_List_Params[I]
          = (Driver_SELM_XML_Handler::SELM_CouplingOperator_ParamsType *)
            malloc(sizeof(Driver_SELM_XML_Handler::SELM_CouplingOperator_ParamsType));
      } /* end of I loop */

      /* parse tags containing list of names + types */
      parseMode = PARSE_MODE_SELM_CouplingOperator_List;

    } else if (qName == xmlTagName_SELM_Interaction_List) {

      /* indicate tag recognized */
      flagRecognizedTag = TAG_RECOGNIZED;

      /* setup the list data structure for parameter data */
      int N = Atz_XML_Helper_ParseData::getIntFromAttr("numEntries",
                                                       xmlAttributes);
      SELM_Interaction_List_Params_N = N;
      SELM_Interaction_List_Params_I = 0;
      SELM_Interaction_List_Params
        = (Driver_SELM_XML_Handler::SELM_Interaction_ParamsType **)
          malloc(sizeof(Driver_SELM_XML_Handler::SELM_Interaction_ParamsType *)* N);

      /* allocate the entries in the list */
      for (int I = 0; I < N; I++) {
        SELM_Interaction_List_Params[I]
          = (Driver_SELM_XML_Handler::SELM_Interaction_ParamsType *)
            malloc(sizeof(Driver_SELM_XML_Handler::SELM_Interaction_ParamsType));
      } /* end of I loop */

      /* parse tags containing list of names + types */
      parseMode = PARSE_MODE_SELM_Interaction_List;

    } else if (qName == xmlTagName_SELM_Integrator) {

      /* indicate tag recognized */
      flagRecognizedTag = TAG_RECOGNIZED;

      /* allocate data structure for parameter data */
      SELM_Integrator_Params
        = (Driver_SELM_XML_Handler::SELM_Integrator_ParamsType *)
          malloc(sizeof(Driver_SELM_XML_Handler::SELM_Integrator_ParamsType));

      /* parse tags containing list of names + types */
      parseMode = PARSE_MODE_SELM_Integrator;

    } else if (qName == xmlTagName_flagWriteSimulationData) {
      flagRecognizedTag = TAG_RECOGNIZED;
    } else if (qName == xmlTagName_saveSkipSimulationData) {
      flagRecognizedTag = TAG_RECOGNIZED;
    } else { /* unrecognized tags skip (avoid sub-tags triggering something here) */
      flagRecognizedTag = TAG_NOT_RECOGNIZED;
    }

    break;

  case PARSE_MODE_SELM_Lagrangian_List:

    if (qName == xmlTagName_SELM_Lagrangian) {
      flagRecognizedTag = TAG_RECOGNIZED;
    } else if (qName == xmlTagName_SELM_LagrangianName) {
      flagRecognizedTag = TAG_RECOGNIZED;
    } else if (qName == xmlTagName_SELM_LagrangianTypeStr) {
      flagRecognizedTag = TAG_RECOGNIZED;
    }

    break;

  case PARSE_MODE_SELM_Eulerian_List:

    if (qName == xmlTagName_SELM_Eulerian) {
      flagRecognizedTag = TAG_RECOGNIZED;
    } else if (qName == xmlTagName_SELM_EulerianName) {
      flagRecognizedTag = TAG_RECOGNIZED;
    } else if (qName == xmlTagName_SELM_EulerianTypeStr) {
      flagRecognizedTag = TAG_RECOGNIZED;
    }

    break;

  case PARSE_MODE_SELM_Interaction_List:

    if (qName == xmlTagName_SELM_Interaction) {
      flagRecognizedTag = TAG_RECOGNIZED;
    } else if (qName == xmlTagName_SELM_InteractionName) {
      flagRecognizedTag = TAG_RECOGNIZED;
    } else if (qName == xmlTagName_SELM_InteractionTypeStr) {
      flagRecognizedTag = TAG_RECOGNIZED;
    }

    break;

  case PARSE_MODE_SELM_CouplingOperator_List:

    if (qName == xmlTagName_SELM_CouplingOperator) {
      flagRecognizedTag = TAG_RECOGNIZED;
    } else if (qName == xmlTagName_SELM_CouplingOperatorName) {
      flagRecognizedTag = TAG_RECOGNIZED;
    } else if (qName == xmlTagName_SELM_CouplingOperatorTypeStr) {
      flagRecognizedTag = TAG_RECOGNIZED;
    }

    break;

  case PARSE_MODE_SELM_Integrator:

    if (qName == xmlTagName_SELM_Integrator) {
      flagRecognizedTag = TAG_RECOGNIZED;
    } else if (qName == xmlTagName_SELM_IntegratorName) {
      flagRecognizedTag = TAG_RECOGNIZED;
    } else if (qName == xmlTagName_SELM_IntegratorTypeStr) {
      flagRecognizedTag = TAG_RECOGNIZED;
    }

    break;

  } /* end of switch */

  if (flagRecognizedTag == TAG_NOT_RECOGNIZED) {
    Atz_XML_SAX_Handler_Multilevel*    sourceHandler_Multilevel;
    Atz_XML_SAX_DataHandler*           dataHandler;

    sourceHandler_Multilevel = dynamic_cast<Atz_XML_SAX_Handler_Multilevel*>(sourceHandler);
    dataHandler              = new Atz_XML_Helper_Handler_SkipNextTag();
    sourceHandler_Multilevel->parseNextTagWithDataHandler(dataHandler);
  } /* end tag not recognized */

}

void Driver_SELM_XML_Handler::XML_characters(string xmlString_in, Atz_XML_SAX_DataHandler* sourceHandler) {
  xmlString += xmlString_in;
}

void Driver_SELM_XML_Handler::XML_endElement(string qName, Atz_XML_SAX_DataHandler* sourceHandler) {

  const char *error_str_code = "Driver_SELM_XML_Handler.cpp";
  const char *error_str_func = "XML_endElement()";

  if (driver_SELM_Data != NULL) {
    /* */
  }

  switch (parseMode) {

  case PARSE_MODE_FIX_SELM:

    if (qName == xmlTagName_SELM_Version) {
      driver_SELM_Data->SELM_Version = Atz_XML_Helper_ParseData::getIntFromAttr(xmlAttributes);
    } else if (qName == xmlTagName_SELM_Run_Description) {
      driver_SELM_Data->SELM_Run_Description = *Atz_XML_Helper_ParseData::getStringFromAttr(xmlAttributes);
    } else if (qName == xmlTagName_SELM_BasePath) {
    
      int rv;    
      string *cpp_str = Atz_XML_Helper_ParseData::getStringFromAttr(xmlAttributes);
      const char *c_str = cpp_str->c_str();
      int len         = strlen(c_str);
            
      driver_SELM_Data->SELM_BasePath = (char *)malloc(sizeof(char)*(len + 1));
      strcpy(driver_SELM_Data->SELM_BasePath, c_str);
      
      // make the directory
      printf("Making directory: %s \n",driver_SELM_Data->SELM_BasePath);
      rv = mkdir(driver_SELM_Data->SELM_BasePath,0755);
      
      if ((rv == -1) && (errno != EEXIST)) {  
        stringstream message;
        message << "Failed making directory path = " << driver_SELM_Data->SELM_BasePath << endl;
        SELM_Package::packageError(error_str_code, error_str_func, message);
      }

      stringstream dir_output;
      dir_output << driver_SELM_Data->SELM_BasePath << "/sim_data";
      const char *c_str2 = dir_output.str().c_str();
      int len2         = strlen(c_str2);
      driver_SELM_Data->SELM_dir_sim_data = (char *)malloc(sizeof(char)*(len2 + 1));
      strcpy(driver_SELM_Data->SELM_dir_sim_data, c_str2);
      
      // make the directory
      printf("Making directory: %s \n",driver_SELM_Data->SELM_dir_sim_data);
      rv = mkdir(driver_SELM_Data->SELM_dir_sim_data,0755);
      
      if ((rv == -1) && (errno != EEXIST)) {  
        stringstream message;
        message << "Failed making directory path = " << driver_SELM_Data->SELM_dir_sim_data << endl;
        SELM_Package::packageError(error_str_code, error_str_func, message);
      }

    } else if (qName == xmlTagName_SELM_BaseFilename) {
      string *cpp_str = Atz_XML_Helper_ParseData::getStringFromAttr(xmlAttributes);
      const char *c_str = cpp_str->c_str();
      int len         = strlen(c_str);
      driver_SELM_Data->SELM_BaseFilename = (char *)malloc(sizeof(char)*(len + 1));
      strcpy(driver_SELM_Data->SELM_BaseFilename, c_str);
    } else if (qName == xmlTagName_SELM_Seed) {
      driver_SELM_Data->SELM_Seed = Atz_XML_Helper_ParseData::getIntFromAttr(xmlAttributes);
    } else if (qName == xmlTagName_FixSELM) {

      /* Setup the fix_SELM data structures that include:
       * Eulerian, Lagrangian, CouplingOp, Integrator.
       *
       * This may require parsing additional files and
       * cross-referencing named objects, etc...
       */
      setupDriverSELM_From_Params();

    }

    break;

  case PARSE_MODE_SELM_Lagrangian_List:

    if (qName == xmlTagName_SELM_LagrangianName) {
      int I = SELM_Lagrangian_List_Params_I;
      strcpy((char *)SELM_Lagrangian_List_Params[I]->SELM_LagrangianName,
             Atz_XML_Helper_ParseData::getStringFromAttr(xmlAttributes)->c_str());
    } else if (qName == xmlTagName_SELM_LagrangianTypeStr) {
      int I = SELM_Lagrangian_List_Params_I;
      strcpy((char *)SELM_Lagrangian_List_Params[I]->SELM_LagrangianTypeStr,
             Atz_XML_Helper_ParseData::getStringFromAttr(xmlAttributes)->c_str());
    } else if (qName == xmlTagName_SELM_Lagrangian) {

      /* done parsing this data entry, so increment the counter */
      SELM_Lagrangian_List_Params_I++;

    } else if (qName == xmlTagName_SELM_Lagrangian_List) {

      /* list is done, so return parser back to default mode */
      parseMode = PARSE_MODE_FIX_SELM;

    }

    break;

  case PARSE_MODE_SELM_Eulerian_List:

    if (qName == xmlTagName_SELM_EulerianName) {
      int I = SELM_Eulerian_List_Params_I;
      strcpy((char *)SELM_Eulerian_List_Params[I]->SELM_EulerianName,
             Atz_XML_Helper_ParseData::getStringFromAttr(xmlAttributes)->c_str());
    } else if (qName == xmlTagName_SELM_EulerianTypeStr) {
      int I = SELM_Eulerian_List_Params_I;
      strcpy((char *)SELM_Eulerian_List_Params[I]->SELM_EulerianTypeStr,
             Atz_XML_Helper_ParseData::getStringFromAttr(xmlAttributes)->c_str());
    } else if (qName == xmlTagName_SELM_Eulerian) {

      /* done parsing this data entry, so increment the counter */
      SELM_Eulerian_List_Params_I++;

    } else if (qName == xmlTagName_SELM_Eulerian_List) {

      /* list is done, so return parser back to default mode */
      parseMode = PARSE_MODE_FIX_SELM;

    }

    break;

  case PARSE_MODE_SELM_CouplingOperator_List:

    if (qName == xmlTagName_SELM_CouplingOperatorName) {
      int I = SELM_CouplingOperator_List_Params_I;
      strcpy((char *)SELM_CouplingOperator_List_Params[I]->SELM_CouplingOperatorName,
             Atz_XML_Helper_ParseData::getStringFromAttr(xmlAttributes)->c_str());
    } else if (qName == xmlTagName_SELM_CouplingOperatorTypeStr) {
      int I = SELM_CouplingOperator_List_Params_I;
      strcpy((char *)SELM_CouplingOperator_List_Params[I]->SELM_CouplingOperatorTypeStr,
             Atz_XML_Helper_ParseData::getStringFromAttr(xmlAttributes)->c_str());
    } else if (qName == xmlTagName_SELM_CouplingOperator) {

      /* done parsing this data entry, so increment the counter */
      SELM_CouplingOperator_List_Params_I++;

    } else if (qName == xmlTagName_SELM_CouplingOperator_List) {

      /* list is done, so return parser back to default mode */
      parseMode = PARSE_MODE_FIX_SELM;

    }

    break;

  case PARSE_MODE_SELM_Interaction_List:

    if (qName == xmlTagName_SELM_InteractionName) {
      int I = SELM_Interaction_List_Params_I;
      strcpy((char *) SELM_Interaction_List_Params[I]->SELM_InteractionName,
             Atz_XML_Helper_ParseData::getStringFromAttr(xmlAttributes)->c_str());
    } else if (qName == xmlTagName_SELM_InteractionTypeStr) {
      int I = SELM_Interaction_List_Params_I;
      strcpy((char *) SELM_Interaction_List_Params[I]->SELM_InteractionTypeStr,
             Atz_XML_Helper_ParseData::getStringFromAttr(xmlAttributes)->c_str());
    } else if (qName == xmlTagName_SELM_Interaction) {

      /* done parsing this data entry, so increment the counter */
      SELM_Interaction_List_Params_I++;

    } else if (qName == xmlTagName_SELM_Interaction_List) {

      /* list is done, so return parser back to default mode */
      parseMode = PARSE_MODE_FIX_SELM;

    }

    break;

  case PARSE_MODE_SELM_Integrator:

    if (qName == xmlTagName_SELM_IntegratorName) {
      strcpy((char *)SELM_Integrator_Params->SELM_IntegratorName,
             Atz_XML_Helper_ParseData::getStringFromAttr(xmlAttributes)->c_str());
    } else if (qName == xmlTagName_SELM_IntegratorTypeStr) {
      strcpy((char *)SELM_Integrator_Params->SELM_IntegratorTypeStr,
             Atz_XML_Helper_ParseData::getStringFromAttr(xmlAttributes)->c_str());
    } else if (qName == xmlTagName_SELM_Integrator) {

      /* list is done, so return parser back to default mode */
      parseMode = PARSE_MODE_FIX_SELM;

    }

    break;

  } /* end switch */

}

void *Driver_SELM_XML_Handler::XML_getData() { /* gets data from parsing the XML */
  return driver_SELM_Data;
}


void Driver_SELM_XML_Handler::setupDriverSELM_From_Params() {

  const char *error_str_code = "Driver_SELM_XML_Handler.cpp";
  const char *error_str_func = "setupDriverSELM_From_Params()";

  int  N;
  char filename[10000];

  /* == setup the Lagrangian data structures */
  driver_SELM_Data->SELM_Lagrangian_List_N = SELM_Lagrangian_List_Params_N;
  N = driver_SELM_Data->SELM_Lagrangian_List_N;
  driver_SELM_Data->SELM_Lagrangian_List
    = (SELM_Lagrangian **) malloc(sizeof(SELM_Lagrangian *)*N);

  /* -- loop over SELM_Lagrangian specifications and parse the
   * required data from XML files for SELM_Lagrangian objects.
   */
  N = SELM_Lagrangian_List_Params_N;
  for (int I = 0; I < N; I++) {
    sprintf(filename, "%s%s.SELM_Lagrangian",
            driver_SELM_Data->SELM_BasePath,
            SELM_Lagrangian_List_Params[I]->SELM_LagrangianName);

    SELM_Lagrangian_Delegator_XML_Handler *SELM_Lagrangian_DataHandler
      = new SELM_Lagrangian_Delegator_XML_Handler();

    Atz_XML_SAX_Handler_Multilevel *dataHandler
      = new Atz_XML_SAX_Handler_Multilevel(SELM_Lagrangian_DataHandler);

    dataHandler->parseCurrentScopeWithDataHandler(SELM_Lagrangian_DataHandler); /* ensure current scope parsed with handler */

    Atz_XML_Parser::parse(filename, dataHandler);

    SELM_Lagrangian *lagrangian
      = (SELM_Lagrangian *) SELM_Lagrangian_DataHandler->XML_getData();

    if (lagrangian == NULL) {

      /* indicates some problem */
      stringstream message;
      message << "Lagrangian data returned is NULL." << endl;
      message << "(lagrangian == NULL)" << endl;
      message << "This could indicate the type is not recognized." << endl;
      message << "filename = " << filename << endl;
      SELM_Package::packageError(error_str_code, error_str_func, message);

    } else {
      lagrangian->setGlobalRefs(driver_SELM_Data->lammps, driver_SELM_Data);

      //delete Atz_XML_SAX_Handler_Multilevel;
      //delete SELM_Lagrangian_DataHandler;
    }

    driver_SELM_Data->SELM_Lagrangian_List[I] = lagrangian;

  }

  /* == setup the Eulerian data structures */
  driver_SELM_Data->SELM_Eulerian_List_N = SELM_Eulerian_List_Params_N;
  N = driver_SELM_Data->SELM_Eulerian_List_N;
  driver_SELM_Data->SELM_Eulerian_List
  = (SELM_Eulerian **) malloc(sizeof(SELM_Eulerian *)*N);

  /* -- loop over SELM_Eulerian specifications and parse the
   * required data from XML files for SELM_Eulerian objects.
   */
  N = SELM_Eulerian_List_Params_N;
  for (int I = 0; I < N; I++) {
    sprintf(filename, "%s%s.SELM_Eulerian",
            driver_SELM_Data->SELM_BasePath,
            SELM_Eulerian_List_Params[I]->SELM_EulerianName);

    SELM_Eulerian_Delegator_XML_Handler *SELM_Eulerian_DataHandler
    = new SELM_Eulerian_Delegator_XML_Handler();

    Atz_XML_SAX_Handler_Multilevel *dataHandler
    = new Atz_XML_SAX_Handler_Multilevel(SELM_Eulerian_DataHandler);

    dataHandler->parseCurrentScopeWithDataHandler(SELM_Eulerian_DataHandler); /* ensure current scope parsed with handler */

    Atz_XML_Parser::parse(filename, dataHandler);

    SELM_Eulerian *eulerian
    = (SELM_Eulerian *) SELM_Eulerian_DataHandler->XML_getData();

    eulerian->setGlobalRefs(driver_SELM_Data->lammps, driver_SELM_Data);

    //delete Atz_XML_SAX_Handler_Multilevel;
    //delete SELM_Eulerian_DataHandler;

    driver_SELM_Data->SELM_Eulerian_List[I] = eulerian;
  }

  /* == setup the CouplingOperators data structures */
  driver_SELM_Data->SELM_CouplingOperator_List_N = SELM_CouplingOperator_List_Params_N;
  N = driver_SELM_Data->SELM_CouplingOperator_List_N;
  driver_SELM_Data->SELM_CouplingOperator_List
  = (SELM_CouplingOperator **) malloc(sizeof(SELM_CouplingOperator *)*N);

  /* -- loop over SELM_CouplingOperators specifications and parse the
   * required data from XML files for SELM_CouplingOperators objects.
   */
  N = SELM_CouplingOperator_List_Params_N;
  for (int I = 0; I < N; I++) {
    sprintf(filename, "%s%s.SELM_CouplingOperator",
            driver_SELM_Data->SELM_BasePath,
            SELM_CouplingOperator_List_Params[I]->SELM_CouplingOperatorName);

    SELM_CouplingOperator_Delegator_XML_Handler *SELM_CouplingOperators_DataHandler
    = new SELM_CouplingOperator_Delegator_XML_Handler(driver_SELM_Data->SELM_Lagrangian_List_N,
                                                      driver_SELM_Data->SELM_Lagrangian_List,
                                                      driver_SELM_Data->SELM_Eulerian_List_N,
                                                      driver_SELM_Data->SELM_Eulerian_List);

    /* specify the base path to use for internal file references during parsing */
    strcpy(SELM_CouplingOperators_DataHandler->basePath, driver_SELM_Data->SELM_BasePath);

    Atz_XML_SAX_Handler_Multilevel *dataHandler
    = new Atz_XML_SAX_Handler_Multilevel(SELM_CouplingOperators_DataHandler);

    dataHandler->parseCurrentScopeWithDataHandler(SELM_CouplingOperators_DataHandler); /* ensure current scope parsed with handler */

    Atz_XML_Parser::parse(filename, dataHandler);

    SELM_CouplingOperator *couplingOp
    = (SELM_CouplingOperator *) SELM_CouplingOperators_DataHandler->XML_getData();

    couplingOp->setGlobalRefs(driver_SELM_Data->lammps, driver_SELM_Data);

    //delete Atz_XML_SAX_Handler_Multilevel;
    //delete SELM_CouplingOperators_DataHandler;

    driver_SELM_Data->SELM_CouplingOperator_List[I] = couplingOp;
  }

  /* == setup the Interaction data structures */
  driver_SELM_Data->SELM_Interaction_List_N = SELM_Interaction_List_Params_N;
  N = driver_SELM_Data->SELM_Interaction_List_N;
  driver_SELM_Data->SELM_Interaction_List
      = (SELM_Interaction **) malloc(sizeof(SELM_Interaction *) * N);

  /* -- loop over SELM_Interaction specifications and parse the
   * required data from XML files for SELM_Interaction objects.
   */
  N = SELM_Interaction_List_Params_N;
  for (int I = 0; I < N; I++) {
    sprintf(filename, "%s%s.SELM_Interaction", driver_SELM_Data->SELM_BasePath,
            SELM_Interaction_List_Params[I]->SELM_InteractionName);

    SELM_Interaction_Delegator_XML_Handler *SELM_Interaction_DataHandler =
        new SELM_Interaction_Delegator_XML_Handler();

    Atz_XML_SAX_Handler_Multilevel *dataHandler =
        new Atz_XML_SAX_Handler_Multilevel(SELM_Interaction_DataHandler);

    dataHandler->parseCurrentScopeWithDataHandler(SELM_Interaction_DataHandler); /* ensure current scope parsed with handler */

    Atz_XML_Parser::parse(filename, dataHandler);

    SELM_Interaction *interaction =
        (SELM_Interaction *) SELM_Interaction_DataHandler->XML_getData();

    interaction->setGlobalRefs(driver_SELM_Data->lammps, driver_SELM_Data);

    //delete Atz_XML_SAX_Handler_Multilevel;
    //delete SELM_Interaction_DataHandler;

    driver_SELM_Data->SELM_Interaction_List[I] = interaction;
  }

  /* == setup the Integrator data structures */
  sprintf(filename, "%s%s.SELM_Integrator",
          driver_SELM_Data->SELM_BasePath,
          SELM_Integrator_Params->SELM_IntegratorName);

  SELM_Integrator_Delegator_XML_Handler *SELM_Integrator_DataHandler
  = new SELM_Integrator_Delegator_XML_Handler();

  Atz_XML_SAX_Handler_Multilevel *dataHandler
  = new Atz_XML_SAX_Handler_Multilevel(SELM_Integrator_DataHandler);

  dataHandler->parseCurrentScopeWithDataHandler(SELM_Integrator_DataHandler); /* ensure current scope parsed with handler */

  Atz_XML_Parser::parse(filename, dataHandler);

  SELM_Integrator *integrator
  = (SELM_Integrator *) SELM_Integrator_DataHandler->XML_getData();

  integrator->setGlobalRefs(driver_SELM_Data->lammps, driver_SELM_Data);

  //delete Atz_XML_SAX_Handler_Multilevel;
  //delete SELM_Integrator_DataHandler;

  driver_SELM_Data->SELM_IntegratorData = integrator;

}
