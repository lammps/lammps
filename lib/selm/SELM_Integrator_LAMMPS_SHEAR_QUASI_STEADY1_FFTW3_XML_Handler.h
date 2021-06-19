/* ----------------------------------------------------------------------
 XML Handling.

 Paul J. Atzberger
 http://atzberger.org/
------------------------------------------------------------------------- */

#ifndef SELM_INTEGRATOR_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_XML_HANDLER_H
#define SELM_INTEGRATOR_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_XML_HANDLER_H

#include "Atz_XML_Package.h"
#include "SELM_Integrator_Delegator_XML_Handler.h"
#include "SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3.h"

namespace LAMMPS_NS {

class SELM_Integrator_Delegator_XML_Handler; /* declare forward reference of class, since this
                                                also refers to the current class below */

class SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_XML_Handler : public Atz_XML_SAX_DataHandler {

 public:
   static const int PARSE_MODE_DEFAULT   = 0;
   static const int PARSE_MODE_RM_SHEAR1 = 1;
   static const int PARSE_MODE_RM_OSC1   = 2;

   static const char* error_str_code;

 public:

  /* ======================== Function prototypes ======================= */
  SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_XML_Handler();
  SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_XML_Handler(SELM_Integrator_Delegator_XML_Handler *delegatorHandler);

  virtual ~SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3_XML_Handler();

  int            parseMode;

  string         xmlTagName_xml;
  string         xmlTagName_SELM_Integrator;
  string         xmlTagName_IntegratorName;
  string         xmlTagName_maxTimeStepIndex;
  string         xmlTagName_deltaT;
  string         xmlTagName_mu;
  string         xmlTagName_rho;
  string         xmlTagName_KB;
  string         xmlTagName_T;
  //string         xmlTagName_flagShearModeStr;

  string         xmlTagName_shearData;

  string         xmlTagName_shearRate;
  string         xmlTagName_shearDir;
  string         xmlTagName_shearVelDir;
  string         xmlTagName_shearDist;

  //string         xmlTagName_RM_OSC1;
  string         xmlTagName_shearOmega;
  string         xmlTagName_shearRateAmplitude;

  string         xmlTagName_flagStochasticDriving;
  string         xmlTagName_flagIncompressibleFluid;
  string         xmlTagName_flagWriteSimulationData;
  string         xmlTagName_saveSkipSimulationData;

  friend class SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3;

  SELM_Integrator_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3 *integrator; /* object to construct */

  Atz_XML::AttributesType *xmlAttributes;
  string xmlString;

 public:

   void setupDataHandler();

   /* XML interface */
   void XML_startDocument(Atz_XML_SAX_DataHandler* sourceHandler);

   void XML_endDocument(Atz_XML_SAX_DataHandler* sourceHandler);

   void XML_startElement(string qName, Atz_XML::AttributesType *attributes, Atz_XML_SAX_DataHandler* sourceHandler);

   void XML_characters(string xmlString, Atz_XML_SAX_DataHandler* sourceHandler);

   void XML_endElement(string qName, Atz_XML_SAX_DataHandler* sourceHandler);

   void *XML_getData(); /* gets data from parsing the XML */

   /* routines special to this handler */
   int getFlagShearModeFromStr(const char *flagShearModeStr);

};

}

#endif
