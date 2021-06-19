/*
 * SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_XML_Handler.cpp
 *
 * Paul J. Atzberger
 * http://atzberger.org/
 *
 */

#include "SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_XML_Handler.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

using namespace LAMMPS_NS;

const char* SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_XML_Handler::error_str_code = "SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_XML_Handler";

SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_XML_Handler::SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_XML_Handler() {

  setupDataHandler();

}

SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_XML_Handler::SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_XML_Handler(SELM_Integrator_Delegator_XML_Handler *delegatorHandler) {

  typedef  SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_ParamsType SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_ParamsType;

  setupDataHandler();

  integrator = new SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3();

   // setup parameters for XML parsing
  integrator->SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Params
    = (SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_ParamsType *)
    malloc(sizeof(SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_ParamsType));
  
  // initialize to zero  
  //memset(integrator->SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Params,
  //       0,sizeof(SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_ParamsType));
      
  strcpy(integrator->nameStr, delegatorHandler->SELM_IntegratorName);
  strcpy(integrator->typeStr, delegatorHandler->SELM_IntegratorTypeStr);

}

void SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_XML_Handler::setupDataHandler() {

  DataHandlerName                         = "Default Data Handler";
  DataHandlerType                         = "SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_XML_Handler";

  xmlTagName_xml                          = "xml";
  xmlTagName_SELM_Integrator              = "SELM_Integrator";
  xmlTagName_IntegratorName               = "IntegratorName";
  xmlTagName_maxTimeStepIndex             = "maxTimeStepIndex";
  xmlTagName_deltaT                       = "deltaT";
  xmlTagName_mu                           = "mu";
  xmlTagName_rho                          = "rho";
  xmlTagName_KB                           = "KB";
  xmlTagName_T                            = "T";
  //xmlTagName_flagShearModeStr             = "flagShearModeStr";

  xmlTagName_shearData                    = "shearData";

  //xmlTagName_RM_SHEAR1                    = "RM_SHEAR1";
  xmlTagName_shearRate                    = "shearRate";
  xmlTagName_shearDir                     = "shearDir";
  xmlTagName_shearVelDir                  = "shearVelDir";
  xmlTagName_shearDist                    = "shearDist";

  //xmlTagName_RM_OSC1                      = "RM_OSC1";
  xmlTagName_shearOmega                   = "shearOmega";
  xmlTagName_shearRateAmplitude           = "shearRateAmplitude";

  xmlTagName_flagStochasticDriving        = "flagStochasticDriving";
  xmlTagName_flagIncompressibleFluid      = "flagIncompressibleFluid";
  xmlTagName_flagWriteSimulationData      = "flagWriteSimulationData";
  xmlTagName_saveSkipSimulationData       = "saveSkipSimulationData";

  integrator = NULL;

  parseMode  = PARSE_MODE_DEFAULT;

}

SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_XML_Handler::~SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_XML_Handler() {
  // TODO Auto-generated destructor stub
}

void SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_XML_Handler::XML_startDocument(Atz_XML_SAX_DataHandler* sourceHandler) {

}

void SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_XML_Handler::XML_endDocument(Atz_XML_SAX_DataHandler* sourceHandler) {

}

void SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_XML_Handler::XML_startElement(string qName, Atz_XML::AttributesType *attributes, Atz_XML_SAX_DataHandler* sourceHandler) {

  const char *error_str_func = "XML_startElement()";

  typedef SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::ShearData_RM_SHEAR1_Type ShearData_RM_SHEAR1_Type;
  typedef SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::ShearData_RM_OSC1_Type   ShearData_RM_OSC1_Type;
  typedef SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_ParamsType SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_ParamsType;

  ShearData_RM_SHEAR1_Type *shearData_RM_SHEAR1;
  ShearData_RM_OSC1_Type   *shearData_RM_OSC1;

  SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_ParamsType *integratorData;

  xmlAttributes = attributes;

  xmlString.clear();

  if (integrator != NULL) {
    integratorData = integrator->SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Params;
  }

  if (qName == xmlTagName_xml) {

  } else if (qName == xmlTagName_SELM_Integrator) {

    /* setup the data structure */
    integrator = new SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3();

  } else if (qName == xmlTagName_IntegratorName) {

  } else if (qName == xmlTagName_maxTimeStepIndex) {

  } else if (qName == xmlTagName_deltaT) {

  } else if (qName == xmlTagName_mu) {

  } else if (qName == xmlTagName_rho) {

  } else if (qName == xmlTagName_KB) {

  } else if (qName == xmlTagName_T) {

  //} else if (qName == xmlTagName_flagShearModeStr) {

  } else if (qName == xmlTagName_shearData) {

    strcpy(integratorData->flagShearModeStr, Atz_XML_Helper_ParseData::getStringFromAttr(xmlAttributes)->c_str());

    integratorData->flagShearMode = getFlagShearModeFromStr(integratorData->flagShearModeStr);

    switch (integratorData->flagShearMode) {

      case SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::SHEAR_MODE_TYPE_RM_SHEAR1:
        parseMode = PARSE_MODE_RM_SHEAR1;
        integratorData->shearData
          = (ShearData_RM_SHEAR1_Type *) malloc(sizeof(ShearData_RM_SHEAR1_Type));

        shearData_RM_SHEAR1
          = (ShearData_RM_SHEAR1_Type *)integratorData->shearData;

        shearData_RM_SHEAR1->shearRate          = 0;
        shearData_RM_SHEAR1->shearDir           = -1;
        shearData_RM_SHEAR1->shearVelDir        = -1;
        shearData_RM_SHEAR1->shearDist          = 0;
        shearData_RM_SHEAR1->shearDist_last     = 0;

      break;

      case SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::SHEAR_MODE_TYPE_RM_OSC1:
        parseMode = PARSE_MODE_RM_OSC1;
        integratorData->shearData
          = (ShearData_RM_OSC1_Type *) malloc(sizeof(ShearData_RM_OSC1_Type));

        shearData_RM_OSC1
          = (ShearData_RM_OSC1_Type *)integratorData->shearData;

        shearData_RM_OSC1->shearRate          = 0;
        shearData_RM_OSC1->shearDir           = -1;
        shearData_RM_OSC1->shearVelDir        = -1;
        shearData_RM_OSC1->shearDist          = 0;
        shearData_RM_OSC1->shearDist_last     = 0;
        shearData_RM_OSC1->shearRateAmplitude = 0;
        shearData_RM_OSC1->shearOmega         = 0;

      break;

      default:
        stringstream message;
        message << "The shear mode specified is not recognized or supported yet" << endl;
        message << "flagShearModeStr = " << integratorData->flagShearModeStr << endl;
        SELM_Package::packageError(error_str_code, error_str_func, message);
        break;

    } /* end switch */

  } else if (qName == xmlTagName_shearRate) {
    int a = 1;

  } else if (qName == xmlTagName_shearDir) {
    int a = 1;

  } else if (qName == xmlTagName_shearVelDir) {
    int a = 1;

  } else if (qName == xmlTagName_shearDist) {
    int a = 1;

  //} else if (qName == xmlTagName_RM_OSC1) {

  } else if (qName == xmlTagName_shearOmega) {

  } else if (qName == xmlTagName_shearRateAmplitude) {

  } else if (qName == xmlTagName_flagStochasticDriving) {

  } else if (qName == xmlTagName_flagIncompressibleFluid) {

  } else if (qName == xmlTagName_flagWriteSimulationData) {

  } else if (qName == xmlTagName_saveSkipSimulationData) {

  } else { /* unrecognized tags skip (avoid sub-tags triggering something here) */
    Atz_XML_SAX_Handler_Multilevel*    sourceHandler_Multilevel;
    Atz_XML_SAX_DataHandler*           dataHandler;

    sourceHandler_Multilevel = dynamic_cast<Atz_XML_SAX_Handler_Multilevel*>(sourceHandler);
    dataHandler              = new Atz_XML_Helper_Handler_SkipNextTag();
    sourceHandler_Multilevel->parseNextTagWithDataHandler(dataHandler);
  }

}

void SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_XML_Handler::XML_characters(string xmlString_in, Atz_XML_SAX_DataHandler* sourceHandler) {
  xmlString += xmlString_in;
}

void SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_XML_Handler::XML_endElement(string qName, Atz_XML_SAX_DataHandler* sourceHandler) {

  const char *error_str_func = "XML_endElement()";

  typedef SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::ShearData_RM_SHEAR1_Type ShearData_RM_SHEAR1_Type;
  typedef SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::ShearData_RM_OSC1_Type   ShearData_RM_OSC1_Type;
  typedef SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_ParamsType SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_ParamsType;

  ShearData_RM_SHEAR1_Type *shearData_RM_SHEAR1;
  ShearData_RM_OSC1_Type   *shearData_RM_OSC1;

  SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_ParamsType *integratorData = NULL;

  if (integrator != NULL) {
    integratorData = integrator->SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_Params;
  } else {
    stringstream message;
    message << "The integrator object was not created yet." << endl;
    SELM_Package::packageError(error_str_code, error_str_func, message);
  }

  switch (parseMode) {

  case PARSE_MODE_DEFAULT:

    if (qName == xmlTagName_IntegratorName) {
      //strcpy(integrator->nameStr, Atz_XML_Helper_ParseData::getStringFromAttr(xmlAttributes)->c_str());
    } else if (qName == xmlTagName_maxTimeStepIndex) {
      integratorData->maxTimeStepIndex = Atz_XML_Helper_ParseData::getIntFromAttr(xmlAttributes);
    } else if (qName == xmlTagName_deltaT) {
      integratorData->deltaT = Atz_XML_Helper_ParseData::getDoubleFromAttr(xmlAttributes);
    } else if (qName == xmlTagName_mu) {
      integratorData->mu = Atz_XML_Helper_ParseData::getDoubleFromAttr(xmlAttributes);
    } else if (qName == xmlTagName_rho) {
      integratorData->rho = Atz_XML_Helper_ParseData::getDoubleFromAttr(xmlAttributes);
    } else if (qName == xmlTagName_KB) {
      integratorData->KB = Atz_XML_Helper_ParseData::getDoubleFromAttr(xmlAttributes);
    } else if (qName == xmlTagName_T) {
      integratorData->T = Atz_XML_Helper_ParseData::getDoubleFromAttr(xmlAttributes);
    //} else if (qName == xmlTagName_flagShearModeStr) {

    } else if (qName == xmlTagName_flagStochasticDriving) {
      integratorData->flagStochasticDriving = Atz_XML_Helper_ParseData::getIntFromAttr(xmlAttributes);
    } else if (qName == xmlTagName_flagIncompressibleFluid) {
      integratorData->flagIncompressibleFluid = Atz_XML_Helper_ParseData::getIntFromAttr(xmlAttributes);
    } else if (qName == xmlTagName_flagWriteSimulationData) {
      integratorData->flagWriteSimulationData = Atz_XML_Helper_ParseData::getIntFromAttr(xmlAttributes);
    } else if (qName == xmlTagName_saveSkipSimulationData) {
      integratorData->saveSkipSimulationData = Atz_XML_Helper_ParseData::getIntFromAttr(xmlAttributes);
    } else { /* unrecognized tags skip (avoid sub-tags triggering something here) */

    }

    break;

  case PARSE_MODE_RM_SHEAR1:

    if (integratorData->shearData != NULL) {
      shearData_RM_SHEAR1
        = (ShearData_RM_SHEAR1_Type *)integratorData->shearData;
    }

    if (qName == xmlTagName_shearRate) {
      shearData_RM_SHEAR1->shearRate   = Atz_XML_Helper_ParseData::getDoubleFromAttr(xmlAttributes);
    } else if (qName == xmlTagName_shearDir) {
      shearData_RM_SHEAR1->shearDir    = Atz_XML_Helper_ParseData::getIntFromAttr(xmlAttributes);
    } else if (qName == xmlTagName_shearVelDir) {
      shearData_RM_SHEAR1->shearVelDir = Atz_XML_Helper_ParseData::getIntFromAttr(xmlAttributes);
    } else if (qName == xmlTagName_shearDist) {
      shearData_RM_SHEAR1->shearDist   = Atz_XML_Helper_ParseData::getDoubleFromAttr(xmlAttributes);
    } else if (qName == xmlTagName_shearData) {
      parseMode = PARSE_MODE_DEFAULT; /* return to default parsing mode */
    }

    break;

  case PARSE_MODE_RM_OSC1:

    if (integratorData->shearData != NULL) {
      shearData_RM_OSC1
      = (ShearData_RM_OSC1_Type *)integratorData->shearData;
    }
    
    if (qName == xmlTagName_shearDir) {
      shearData_RM_OSC1->shearDir = Atz_XML_Helper_ParseData::getIntFromAttr(xmlAttributes);
    } else if (qName == xmlTagName_shearVelDir) {
      shearData_RM_OSC1->shearVelDir = Atz_XML_Helper_ParseData::getIntFromAttr(xmlAttributes);
    } else if (qName == xmlTagName_shearOmega) {
      shearData_RM_OSC1->shearOmega = Atz_XML_Helper_ParseData::getDoubleFromAttr(xmlAttributes);
    } else if (qName == xmlTagName_shearRateAmplitude) {
      shearData_RM_OSC1->shearRateAmplitude =  Atz_XML_Helper_ParseData::getDoubleFromAttr(xmlAttributes);
    } else if (qName == xmlTagName_shearDist) {
      shearData_RM_OSC1->shearDist = Atz_XML_Helper_ParseData::getDoubleFromAttr(xmlAttributes);  
    } else if (qName == xmlTagName_shearData) {
      parseMode = PARSE_MODE_DEFAULT; /* return to default parsing mode */
    }

    break;

  } /* end switch parseMode */


}



int SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_XML_Handler::getFlagShearModeFromStr(const char *flagShearModeStr) {

  if (strcmp(flagShearModeStr, SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::SHEAR_MODE_TYPE_STR_RM_SHEAR1) == 0) {
    return SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::SELM_INTEGRATOR_TYPE_OVERDAMPED1_RM_SHEAR1;
  } else if (strcmp(flagShearModeStr, SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::SHEAR_MODE_TYPE_STR_RM_OSC1) == 0) {
    return SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::SHEAR_MODE_TYPE_RM_OSC1;
  } else {
    return SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3::SHEAR_MODE_TYPE_NULL;
  }

}

void *SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_XML_Handler::XML_getData() { /* gets data from parsing the XML */
  return integrator;
}
