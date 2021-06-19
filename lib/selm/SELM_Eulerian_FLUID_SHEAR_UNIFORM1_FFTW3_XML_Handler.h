/* ----------------------------------------------------------------------

  XML Handling.

  Paul J. Atzberger
  http://atzberger.org/
 
------------------------------------------------------------------------- */

#ifndef SELM_EULERIAN_FLUID_SHEAR_UNIFORM1_FFTW3_XML_HANDLER_H
#define SELM_EULERIAN_FLUID_SHEAR_UNIFORM1_FFTW3_XML_HANDLER_H

#include "Atz_XML_Package.h"
#include "SELM_Eulerian_Delegator_XML_Handler.h"
#include "SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3.h"

namespace LAMMPS_NS {

class SELM_Eulerian_Delegator_XML_Handler; /* declare forward reference of class, since this
                                              also refers to the current class below */

class SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_XML_Handler : public Atz_XML_SAX_DataHandler {

 public:

  /* ======================== Function prototypes ======================= */
  SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_XML_Handler();
  SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_XML_Handler(SELM_Eulerian_Delegator_XML_Handler *delegatorHandler);

  virtual ~SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3_XML_Handler();

  string         xmlTagName_xml;
  string         xmlTagName_SELM_Eulerian;
  string         xmlTagName_EulerianName;
  string         xmlTagName_num_dim;
  string         xmlTagName_numMeshPtsPerDir;
  string         xmlTagName_meshDeltaX;
  string         xmlTagName_meshCenterX0;
  string         xmlTagName_shearDir;
  string         xmlTagName_shearVelDir;
  string         xmlTagName_shearRate;
  string         xmlTagName_shearDist;
  string         xmlTagName_flagWriteSimulationData;
  string         xmlTagName_saveSkipSimulationData;

  friend class SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3;

  SELM_Eulerian_FLUID_SHEAR_UNIFORM1_FFTW3 *eulerian; /* object to construct */

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
