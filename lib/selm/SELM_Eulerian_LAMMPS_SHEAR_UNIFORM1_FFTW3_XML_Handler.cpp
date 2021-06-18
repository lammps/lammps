/*
 * SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_XML_Handler.cpp
 *
 * Paul J. Atzberger
 * http://atzberger.org/
 *
 */

#include "SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_XML_Handler.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

using namespace LAMMPS_NS;

SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_XML_Handler::SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_XML_Handler() {

  setupDataHandler();

}

SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_XML_Handler::SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_XML_Handler(SELM_Eulerian_Delegator_XML_Handler *delegatorHandler) {

  typedef  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3::SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_ParamsType SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_ParamsType;

  setupDataHandler();

  /* allocate the eulerian object and the parameters data structure */
  eulerian = new SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3();

  eulerian->SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Params
    = (SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_ParamsType *)
    malloc(sizeof(SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_ParamsType));
  
  memset(eulerian->SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Params,
         0,sizeof(SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_ParamsType));

  strcpy(eulerian->nameStr, delegatorHandler->SELM_EulerianName);
  strcpy(eulerian->typeStr, delegatorHandler->SELM_EulerianTypeStr);

}

SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_XML_Handler::~SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_XML_Handler() {
  // TODO Auto-generated destructor stub
}

void SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_XML_Handler::setupDataHandler() {

  DataHandlerName                         = "Default Data Handler";
  DataHandlerType                         = "SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_XML_Handler";

  xmlTagName_xml                          = "xml";
  xmlTagName_SELM_Eulerian                = "SELM_Eulerian";
  xmlTagName_EulerianName                 = "EulerianName";
  xmlTagName_num_dim                      = "num_dim";
  xmlTagName_numMeshPtsPerDir             = "numMeshPtsPerDir";
  xmlTagName_meshDeltaX                   = "meshDeltaX";
  xmlTagName_meshCenterX0                 = "meshCenterX0";
  xmlTagName_shearDir                     = "shearDir";
  xmlTagName_shearVelDir                  = "shearVelDir";
  xmlTagName_shearRate                    = "shearRate";
  xmlTagName_shearDist                    = "shearDist";
  xmlTagName_flagWriteSimulationData      = "flagWriteSimulationData";
  xmlTagName_saveSkipSimulationData       = "saveSkipSimulationData";

  xmlTagName_flagWriteFluidVel_VTK        = "flagWriteFluidVel_VTK";
  xmlTagName_flagWriteFluidForce_VTK      = "flagWriteFluidForce_VTK";
  xmlTagName_flagWriteFluidPressure_VTK   = "flagWriteFluidPressure_VTK";

  eulerian = NULL;

}

void SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_XML_Handler::XML_startDocument(Atz_XML_SAX_DataHandler* sourceHandler) {

}

void SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_XML_Handler::XML_endDocument(Atz_XML_SAX_DataHandler* sourceHandler) {

}

void SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_XML_Handler::XML_startElement(string qName, Atz_XML::AttributesType *attributes, Atz_XML_SAX_DataHandler* sourceHandler) {

  typedef SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3::SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_ParamsType SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_ParamsType;

  xmlAttributes = attributes;

  xmlString.clear();

  if (qName == xmlTagName_xml) {

  } else if (qName == xmlTagName_SELM_Eulerian) {

    /* setup the data structure */
    /* allocate the eulerian object and the parameters data structure */
    eulerian = new SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3();

    eulerian->SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Params
      = (SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_ParamsType *)
      malloc(sizeof(SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_ParamsType));
      
    memset(eulerian->SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Params,
           0,sizeof(SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_ParamsType));

  } else if (qName == xmlTagName_EulerianName) {

  } else if (qName == xmlTagName_num_dim) {

  } else if (qName == xmlTagName_numMeshPtsPerDir) {

  } else if (qName == xmlTagName_meshDeltaX) {

  } else if (qName == xmlTagName_meshCenterX0) {

  } else if (qName == xmlTagName_shearDir) {

  } else if (qName == xmlTagName_shearVelDir) {

  } else if (qName == xmlTagName_shearRate) {

  } else if (qName == xmlTagName_shearDist) {

  } else if (qName == xmlTagName_flagWriteSimulationData) {

  } else if (qName == xmlTagName_saveSkipSimulationData) {

  } else if (qName == xmlTagName_flagWriteFluidVel_VTK) {

  } else if (qName == xmlTagName_flagWriteFluidForce_VTK) {

  } else if (qName == xmlTagName_flagWriteFluidPressure_VTK) {


  } else { /* unrecognized tags skip (avoid sub-tags triggering something here) */
    Atz_XML_SAX_Handler_Multilevel*    sourceHandler_Multilevel;
    Atz_XML_SAX_DataHandler*           dataHandler;

    sourceHandler_Multilevel = dynamic_cast<Atz_XML_SAX_Handler_Multilevel*>(sourceHandler);
    dataHandler              = new Atz_XML_Helper_Handler_SkipNextTag();
    sourceHandler_Multilevel->parseNextTagWithDataHandler(dataHandler);
  }

}

void SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_XML_Handler::XML_characters(string xmlString_in, Atz_XML_SAX_DataHandler* sourceHandler) {
  xmlString += xmlString_in;
}

void SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_XML_Handler::XML_endElement(string qName, Atz_XML_SAX_DataHandler* sourceHandler) {

  SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3::SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_ParamsType *eulerianData;

  if (eulerian != NULL) {
    eulerianData = eulerian->SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_Params;
  }

  if (qName == xmlTagName_EulerianName) {
    //strcpy(eulerian->nameStr, Atz_XML_Helper_ParseData::getStringFromAttr(xmlAttributes)->c_str());
  } else if (qName == xmlTagName_num_dim) {
    eulerianData->num_dim = Atz_XML_Helper_ParseData::getIntFromAttr(xmlAttributes);
  } else if (qName == xmlTagName_numMeshPtsPerDir) {
    int N             = eulerianData->num_dim;                  /* indicates size of pre-alloc array */
    int *intArray_ptr = (int *)&eulerianData->numMeshPtsPerDir; /* indicates use pre-alloc array */
    Atz_XML_Helper_ParseData::parseIntArrayFromString(xmlString, &intArray_ptr, &N);
  } else if (qName == xmlTagName_meshDeltaX) {
    eulerianData->meshDeltaX = Atz_XML_Helper_ParseData::getDoubleFromAttr(xmlAttributes);
  } else if (qName == xmlTagName_meshCenterX0) {
    int     N               = eulerianData->num_dim;                 /* indicates size of pre-alloc array */
    double *doubleArray_ptr = (double *)&eulerianData->meshCenterX0; /* indicates use pre-alloc array */
    Atz_XML_Helper_ParseData::parseDoubleArrayFromString(xmlString, &doubleArray_ptr , &N);
  } else if (qName == xmlTagName_shearDir) {
    eulerianData->shearDir = Atz_XML_Helper_ParseData::getIntFromAttr(xmlAttributes);
  } else if (qName == xmlTagName_shearVelDir) {
    eulerianData->shearVelDir = Atz_XML_Helper_ParseData::getIntFromAttr(xmlAttributes);
  } else if (qName == xmlTagName_shearRate) {
    eulerianData->shearRate = Atz_XML_Helper_ParseData::getDoubleFromAttr(xmlAttributes);
  } else if (qName == xmlTagName_shearDist) {
    eulerianData->shearDist = Atz_XML_Helper_ParseData::getDoubleFromAttr(xmlAttributes);
  } else if (qName == xmlTagName_flagWriteSimulationData) {
    eulerian->flagWriteSimulationData = Atz_XML_Helper_ParseData::getIntFromAttr(xmlAttributes);
  } else if (qName == xmlTagName_saveSkipSimulationData) {
    eulerian->saveSkipSimulationData = Atz_XML_Helper_ParseData::getIntFromAttr(xmlAttributes);
  } else if (qName == xmlTagName_flagWriteFluidVel_VTK) {
    eulerian->flagWriteFluidVel_VTK = Atz_XML_Helper_ParseData::getIntFromAttr(xmlAttributes);
  } else if (qName == xmlTagName_flagWriteFluidForce_VTK) {
    eulerian->flagWriteFluidForce_VTK = Atz_XML_Helper_ParseData::getIntFromAttr(xmlAttributes);
  } else if (qName == xmlTagName_flagWriteFluidPressure_VTK) {
    eulerian->flagWriteFluidPressure_VTK = Atz_XML_Helper_ParseData::getIntFromAttr(xmlAttributes);
  } else if (qName == xmlTagName_SELM_Eulerian) {
    //eulerian->setup();
  } else { /* unrecognized tags skip (avoid sub-tags triggering something here) */

  }

}

void *SELM_Eulerian_LAMMPS_SHEAR_UNIFORM1_FFTW3_XML_Handler::XML_getData() { /* gets data from parsing the XML */
  return eulerian;
}
