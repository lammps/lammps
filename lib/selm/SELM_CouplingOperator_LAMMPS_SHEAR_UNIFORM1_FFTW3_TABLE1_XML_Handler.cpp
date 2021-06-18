/*
 * SELM_CouplingOperators_TABLE1_XML_Handler.cpp
 *
 * Paul J. Atzberger
 * http://atzberger.org/
 *
 */

#include "SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1_XML_Handler.h"
#include "SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <malloc.h>

using namespace std;
using namespace LAMMPS_NS;

SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1_XML_Handler::SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1_XML_Handler() {

  setupDataHandler();

}


SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1_XML_Handler::SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1_XML_Handler(SELM_CouplingOperator_Delegator_XML_Handler *delegatorHandler) {

  setupDataHandler();

  couplingOp = new SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1();

  strcpy(couplingOp->nameStr, delegatorHandler->SELM_CouplingOperatorName);
  strcpy(couplingOp->typeStr, delegatorHandler->SELM_CouplingOperatorTypeStr);
  strcpy(basePath,            delegatorHandler->basePath);

}

SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1_XML_Handler::SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1_XML_Handler(SELM_CouplingOperator_Delegator_XML_Handler *delegatorHandler,
                                                                                   int numLagrangianList_in,
                                                                                   SELM_Lagrangian **lagrangianList_in,
                                                                                   int numEulerianList_in,
                                                                                   SELM_Eulerian **eulerianList_in) {
  setupDataHandler();

  couplingOp = new SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1();

  strcpy(couplingOp->nameStr, delegatorHandler->SELM_CouplingOperatorName);
  strcpy(couplingOp->typeStr, delegatorHandler->SELM_CouplingOperatorTypeStr);
  strcpy(basePath,            delegatorHandler->basePath);

  numLagrangianList = numLagrangianList_in;
  lagrangianList    = lagrangianList_in;

  numEulerianList   = numEulerianList_in;
  eulerianList      = eulerianList_in;

}

void SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1_XML_Handler::setupDataHandler() {

  DataHandlerName                           = "Default Data Handler";
  DataHandlerType                           = "SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1_XML_Handler";

  xmlTagName_xml                            = "xml";
  xmlTagName_SELM_CouplingOperator          = "SELM_CouplingOperator";
  xmlTagName_CouplingOperatorName           = "CouplingOperatorName";
  xmlTagName_CouplingOperatorTypeStr        = "CouplingOperatorTypeStr";
  xmlTagName_operatorData                   = "operatorData";

  xmlTagName_numCoupleList                  = "numCoupleList";

  xmlTagName_lagrangianList                 = "lagrangianList";

  xmlTagName_SELM_Lagrangian_Ref            = "SELM_Lagrangian_Ref";

  xmlTagName_LagrangianName                 = "LagrangianName";
  xmlTagName_LagrangianTypeStr              = "LagrangianTypeStr";

  xmlTagName_eulerianList                   = "eulerianList";

  xmlTagName_SELM_Eulerian_Ref              = "SELM_Eulerian_Ref";

  xmlTagName_EulerianName                   = "EulerianName";
  xmlTagName_EulerianTypeStr                = "EulerianTypeStr";

  xmlTagName_T_KERNEL_1_weightTableFilename = "weightTableFilename";
  xmlTagName_flagWriteSimulationData        = "flagWriteSimulationData";
  xmlTagName_saveSkipSimulationData         = "saveSkipSimulationData";

  couplingOp                   = NULL;
  flagDeterminedCouplingOpType = 0;

  strcpy(basePath,"");

}


SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1_XML_Handler::~SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1_XML_Handler() {
  // TODO Auto-generated destructor stub
}

void SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1_XML_Handler::XML_startDocument(Atz_XML_SAX_DataHandler* sourceHandler) {

}

void SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1_XML_Handler::XML_endDocument(Atz_XML_SAX_DataHandler* sourceHandler) {

}

void SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1_XML_Handler::XML_startElement(string qName, Atz_XML::AttributesType *attributes, Atz_XML_SAX_DataHandler* sourceHandler) {

  xmlAttributes = attributes;

  xmlString.clear();

  if (qName == xmlTagName_xml) {

  } else if (qName == xmlTagName_SELM_CouplingOperator) {

    /* setup the data structure */
    couplingOp = new SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1();
    flagDeterminedCouplingOpType = 0;

  } else if (qName == xmlTagName_CouplingOperatorName) {

  } else if (qName == xmlTagName_CouplingOperatorTypeStr) {

  } else if (qName == xmlTagName_operatorData) {
    strcpy(couplingOp->operatorTypeStr, Atz_XML_Helper_ParseData::getStringFromAttr(xmlAttributes)->c_str());
    couplingOp->operatorType     = couplingOp->getOperatorTypeFromStr(couplingOp->operatorTypeStr);
    flagDeterminedCouplingOpType = 1;
  } else if (qName == xmlTagName_T_KERNEL_1_weightTableFilename) {

  } else if (qName == xmlTagName_flagWriteSimulationData) {

  } else if (qName == xmlTagName_saveSkipSimulationData) {

  } else if (qName == xmlTagName_numCoupleList) {

  } else if (qName == xmlTagName_lagrangianList) {

  } else if (qName == xmlTagName_SELM_Lagrangian_Ref) {

  } else if (qName == xmlTagName_LagrangianName) {

  } else if (qName == xmlTagName_LagrangianTypeStr) {

  } else if (qName == xmlTagName_eulerianList) {

  } else if (qName == xmlTagName_SELM_Eulerian_Ref) {

  } else if (qName == xmlTagName_EulerianName) {

  } else if (qName == xmlTagName_EulerianTypeStr) {

  } else { /* unrecognized tags skip (avoid sub-tags triggering something here) */
    Atz_XML_SAX_Handler_Multilevel*    sourceHandler_Multilevel;
    Atz_XML_SAX_DataHandler*           dataHandler;

    sourceHandler_Multilevel = dynamic_cast<Atz_XML_SAX_Handler_Multilevel*>(sourceHandler);
    dataHandler              = new Atz_XML_Helper_Handler_SkipNextTag();
    sourceHandler_Multilevel->parseNextTagWithDataHandler(dataHandler);
  }

}

void SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1_XML_Handler::XML_characters(string xmlString_in, Atz_XML_SAX_DataHandler* sourceHandler) {
  xmlString += xmlString_in;
}

void SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1_XML_Handler::XML_endElement(string qName, Atz_XML_SAX_DataHandler* sourceHandler) {

  const char *error_str_code = "SELM_CouplingOperators_TABLE1_XML_Handler.cpp";
  const char *error_str_func = "XML_endElement()";

  if (couplingOp != NULL) {
    /* */
  }

  if (qName == xmlTagName_CouplingOperatorName) {
    strcpy(couplingOp->nameStr, Atz_XML_Helper_ParseData::getStringFromAttr(xmlAttributes)->c_str());
  } else if (qName == xmlTagName_CouplingOperatorTypeStr) {
    strcpy(couplingOp->typeStr, Atz_XML_Helper_ParseData::getStringFromAttr(xmlAttributes)->c_str());
  } else if (qName == xmlTagName_numCoupleList) {

    couplingOp->numCoupleList = Atz_XML_Helper_ParseData::getIntFromAttr(xmlAttributes);

    int N = couplingOp->numCoupleList;

    lagrangianNameList      = (char **) malloc(sizeof(char *)*N);
    lagrangianNameList_I    = 0;
    lagrangianTypeStrList   = (char **) malloc(sizeof(char *)*N);
    lagrangianTypeStrList_I = 0;

    eulerianNameList        = (char **) malloc(sizeof(char *)*N);
    eulerianNameList_I      = 0;
    eulerianTypeStrList     = (char **) malloc(sizeof(char *)*N);
    eulerianTypeStrList_I   = 0;

  } else if (qName == xmlTagName_LagrangianName) {
    int I = lagrangianNameList_I;
    const char *ss = Atz_XML_Helper_ParseData::getStringFromAttr(xmlAttributes)->c_str();
    int nn = strlen(ss);
    lagrangianNameList[I] = (char *)malloc(sizeof(char)*(nn + 1));
    strcpy(lagrangianNameList[I],ss);
    lagrangianNameList_I++;
  } else if (qName == xmlTagName_LagrangianTypeStr) {
    int I = lagrangianTypeStrList_I;
    const char *ss = Atz_XML_Helper_ParseData::getStringFromAttr(xmlAttributes)->c_str();
    int nn = strlen(ss);
    lagrangianTypeStrList[I] = (char *)malloc(sizeof(char)*(nn + 1));
    //strcpy(lagrangianTypeStrList[I], Atz_XML_Helper_ParseData::getStringFromAttr(xmlAttributes)->c_str());
    strcpy(lagrangianTypeStrList[I],ss);
    lagrangianTypeStrList_I++;
  } else if (qName == xmlTagName_lagrangianList) {

    /* resolve the name references to the master list of lagrangian DOF */
    int N = couplingOp->numCoupleList;
    int M = numLagrangianList;

    couplingOp->lagrangianList = (SELM_Lagrangian **) malloc(sizeof(SELM_Lagrangian *)*N);

    for (int k = 0; k < N; k++) {

      char *nameRef = lagrangianNameList[k];

      couplingOp->lagrangianList[k] = NULL;

      for (int j = 0; j < M; j++) {

        SELM_Lagrangian *lagrangian = lagrangianList[j];
        char            *nameMaster = lagrangian->nameStr;

        if (strcmp(nameRef, nameMaster) == 0) {
          couplingOp->lagrangianList[k] = lagrangian;
        }

      } /* end j loop */

      free(lagrangianNameList[k]);
      free(lagrangianTypeStrList[k]);

    } /* end k loop */

    free(lagrangianNameList);
    free(lagrangianTypeStrList);

  } else if (qName == xmlTagName_EulerianName) {
    int I = eulerianNameList_I;
    const char *ss = Atz_XML_Helper_ParseData::getStringFromAttr(xmlAttributes)->c_str();
    int nn = strlen(ss);
    eulerianNameList[I] = (char *)malloc(sizeof(char)*(nn + 1));
    strcpy(eulerianNameList[I], ss);
    eulerianNameList_I++;
  } else if (qName == xmlTagName_EulerianTypeStr) {
    int I = eulerianTypeStrList_I;
    const char *ss = Atz_XML_Helper_ParseData::getStringFromAttr(xmlAttributes)->c_str();
    int nn = strlen(ss);
    eulerianTypeStrList[I] = (char *)malloc(sizeof(char)*(nn + 1));
    strcpy(eulerianTypeStrList[I], ss);
    eulerianTypeStrList_I++;
  } else if (qName == xmlTagName_eulerianList) {

    /* resolve the name references to the master list of eulerian DOF */
    int N = couplingOp->numCoupleList;
    int M = numEulerianList;

    couplingOp->eulerianList = (SELM_Eulerian **) malloc(sizeof(SELM_Eulerian *)*N);

    for (int k = 0; k < N; k++) {

      char *nameRef = eulerianNameList[k];

      couplingOp->eulerianList[k] = NULL;

      for (int j = 0; j < M; j++) {

        SELM_Eulerian *eulerian = eulerianList[j];
        char *nameMaster        = eulerian->nameStr;

        if (strcmp(nameRef,nameMaster) == 0) {
          couplingOp->eulerianList[k] = eulerian;
        }

      } /* end j loop */

      free(eulerianNameList[k]);
      free(eulerianTypeStrList[k]);

    } /* end k loop */

    free(eulerianNameList);
    free(eulerianTypeStrList);

  } else if (qName == xmlTagName_operatorData) {

  } else if (flagDeterminedCouplingOpType) {

    /* switch based on the type of operator if determined */
    switch (couplingOp->operatorType) {

    case SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::OPERATOR_TYPE_T_KERNEL_1:

      SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::operatorDataType_T_KERNEL_1 *opData;

      /* check is operator data allocated yet */
      if (couplingOp->operatorData == NULL) {
        couplingOp->operatorData
          = (SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::operatorDataType_T_KERNEL_1 *)
            malloc(sizeof(SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::operatorDataType_T_KERNEL_1));
      }

      opData = (SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1::operatorDataType_T_KERNEL_1 *) couplingOp->operatorData;

      /* read the weight filename */
      if (qName == xmlTagName_T_KERNEL_1_weightTableFilename) {
        sprintf(opData->weightTableFilename, "%s%s",
                basePath, Atz_XML_Helper_ParseData::getStringFromAttr(xmlAttributes)->c_str());
        opData->weightTable = NULL;
        couplingOp->readWeightTable(opData->weightTableFilename,
                                    &opData->weightTable);
      } else {

      }

      break;

    default:
      stringstream message;
      message << "Invalid operator type was specified." << endl;
      message << "operatorTypeStr = " << couplingOp->operatorTypeStr <<  endl;
      message << "operatorType = " << couplingOp->operatorType << endl;
      message << "Case may not be implemented yet." << endl;
      SELM_Package::packageError(error_str_code, error_str_func, message);
    } /* end of switch */

  } else { /* unrecognized tags skip (avoid sub-tags triggering something here) */

  }

}

void *SELM_CouplingOperator_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1_XML_Handler::XML_getData() { /* gets data from parsing the XML */
  return couplingOp;
}
