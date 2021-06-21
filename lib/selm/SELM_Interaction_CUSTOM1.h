/* ----------------------------------------------------------------------

 Custom interactions.

 Paul J. Atzberger
 http://atzberger.org/
 
------------------------------------------------------------------------- */

#ifndef SELM_INTERACTION_CUSTOM1_H
#define SELM_INTERACTION_CUSTOM1_H

#ifndef FFT_FFTW3 /* determine if FFTW package specified */
  #warning "No FFTW package specified for SELM codes.  The SELM functionality will be disabled."
#else
  #define USE_PACKAGE_FFTW3
#endif

#include "SELM_Interaction.h"
#include "SELM_Interaction_Types.h"

#include "SELM_Lagrangian.h"


namespace LAMMPS_NS {

class SELM_Interaction_CUSTOM1 : public SELM_Interaction {

 public:

  /* ======================== Constants ======================= */
  static const int   TYPE;
  static const char* TYPE_STR;

  static const char *error_str_code;

  /* ======================== Function prototypes ======================= */
  SELM_Interaction_CUSTOM1();
  SELM_Interaction_CUSTOM1(int, char **);
  SELM_Interaction_CUSTOM1(class LAMMPS *lmps, class DriverSELM *fx);

  // ~SELM_Interaction_CUSTOM1();

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

  int               numMembers;
  SELM_Lagrangian **membraneList_lagrangianI1;
  int              *memberList_ptI1;

  paramDataType    *paramData;

};

}

#endif
