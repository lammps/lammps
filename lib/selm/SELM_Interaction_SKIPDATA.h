/* ----------------------------------------------------------------------
 Custom interaction.

 Paul J. Atzberger
 http://atzberger.org/

------------------------------------------------------------------------- */

#ifndef SELM_INTERACTION_SKIPDATA_H
#define SELM_INTERACTION_SKIPDATA_H

#ifndef FFT_FFTW3 /* determine if FFTW package specified */
  #warning "No FFTW package specified for SELM codes.  The SELM functionality will be disabled."
#else
  #define USE_PACKAGE_FFTW3
#endif

#ifndef DYNAMIC_LINK_PLUG_IN /* determine if plug-in features enabled */
  #warning "Dynamic link plug-ins for the SELM package are no enabled.  Must specify compiler flag DYNAMIC_LINK_PLUG_IN to use these features. "
#else
  //#warning "Dynamic link plug-ins for the SELM package are disabled.  Must specify compiler flag DYNAMIC_LINK_PLUG_IN to use these features. "
#endif


#include "SELM_Interaction.h"
#include "SELM_Interaction_Types.h"

#include "SELM_Lagrangian.h"


namespace LAMMPS_NS {

class SELM_Interaction_SKIPDATA : public SELM_Interaction {

 public:

  /* ======================== Constants ======================= */
  static const int   TYPE;
  static const char* TYPE_STR;

  static const char *error_str_code;

  /* ======================== Function prototypes ======================= */
  SELM_Interaction_SKIPDATA();
  SELM_Interaction_SKIPDATA(int, char **);
  SELM_Interaction_SKIPDATA(class LAMMPS *lmps, class DriverSELM *fx);

  // ~SELM_Interaction_SKIPDATA();

  void init();

  void init_library();

  void parse_ParameterFile(const char *filename);

  void setup();

  void computeForceAndEnergy();

  void writeSimulationDataToDisk(const char *baseFilename, int timeIndex);



  /* ======================== Data structure type definitions ======================== */
  typedef map<string, void *>   paramDataType;

  typedef pair<string, void *>   paramDataPairType;

  /* =============================== Variables =========================== */
  string            libName;

  paramDataType    *paramData;

};

}

#endif
