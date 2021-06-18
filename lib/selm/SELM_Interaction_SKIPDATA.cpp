/* ----------------------------------------------------------------------

  Stochastic Eulerian Lagrangian Methods

  Paul J. Atzberger
  http://atzberger.org/
 
------------------------------------------------------------------------- */

#include <cstdlib>
#include "stdio.h"
#include "string.h"
#include "wrapper_selm.h"

#include "SELM_Parser1.h"

#include "SELM_Interaction.h"
#include "SELM_Interaction_SKIPDATA.h"

#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "lammps.h"
#include "domain.h"
#include "wrapper_selm.h"

#include <malloc.h>

#ifdef DYNAMIC_LINK_PLUG_IN /* determine if plug-in features enabled */
  #include <iostream>
  #include <dlfcn.h>
#endif


using namespace LAMMPS_NS;
using std::cout;
using std::cerr;

/* ---------------------------------------------------------------------- */
const int   SELM_Interaction_SKIPDATA::TYPE     = SELM_Interaction_Types::TYPE_SKIPDATA;
const char* SELM_Interaction_SKIPDATA::TYPE_STR = SELM_Interaction_Types::TYPE_STR_SKIPDATA;

const char *SELM_Interaction_SKIPDATA::error_str_code = "SELM_Interaction_SKIPDATA.cpp";

SELM_Interaction_SKIPDATA::SELM_Interaction_SKIPDATA() {

  /* WARNING: may need to call parent constructor */

  /* generic initialization */
  init();

}



SELM_Interaction_SKIPDATA::SELM_Interaction_SKIPDATA(int narg, char **arg) : SELM_Interaction(narg, arg) {

  const char *error_str_func = "SELM_Interaction_SKIPDATA()";

  /* generic initialization */
  init();

  /*
  printf("Creating SELM_Interaction_SKIPDATA : not yet implemented. \n");
  exit(1);
  */

}


SELM_Interaction_SKIPDATA::SELM_Interaction_SKIPDATA(class LAMMPS *lmps, class DriverSELM *fx) {

  /* generic initialization */
  init();

  setGlobalRefs(lammps, driverSELM);

}

void SELM_Interaction_SKIPDATA::init() {

  type = SELM_Interaction_SKIPDATA::TYPE;
  strcpy(typeStr, SELM_Interaction_SKIPDATA::TYPE_STR);

  strcpy(nameStr, "No Name");

  /*
  num_dim             = 3;

  numControlPts       = 0;
  numControlPts_alloc = 0;

  ptsX                = NULL;
  pt_Vel              = NULL;

  atomID              = NULL;
  moleculeID          = NULL;
  typeID              = NULL;
  atomMass            = NULL;

  pt_Energy           = 0;
  pt_Force            = NULL;

  pt_type             = NULL;
  pt_type_extras      = NULL;

  numEntriesOpGammaVel = 0;
  opGammaVel           = NULL;
  */

  setGlobalRefs(NULL,NULL);

}

void SELM_Interaction_SKIPDATA::init_library() {

  // do nothing

}

void SELM_Interaction_SKIPDATA::parse_ParameterFile(const char *filename) {

}

void SELM_Interaction_SKIPDATA::setup() {

  const char *error_str_func = "setup()";

  /* load and initialize the dynamic link library */
  init_library();

}

void SELM_Interaction_SKIPDATA::computeForceAndEnergy() {

}

void SELM_Interaction_SKIPDATA::writeSimulationDataToDisk(const char *baseFilename, int timeIndex) {

}
