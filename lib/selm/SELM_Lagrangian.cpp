/* ----------------------------------------------------------------------

  Stochastic Eulerian Lagrangian Methods

  Paul J. Atzberger
  http://atzberger.org/

------------------------------------------------------------------------- */

#include "driver_selm.h"
#include "SELM_Lagrangian.h"

#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "lammps.h"

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <malloc.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

SELM_Lagrangian::SELM_Lagrangian() {

  const char *error_str_code = "SELM_Lagrangian.cpp";
  const char *error_str_func = "SELM_Lagrangian()";

  flagWriteSimulationData = 0;
  saveSkipSimulationData  = 0;

}

SELM_Lagrangian::SELM_Lagrangian(int narg, char **arg) {

  const char *error_str_code = "SELM_Lagrangian.cpp";
  const char *error_str_func = "SELM_Lagrangian(narg,arg)";

  flagWriteSimulationData = 0;
  saveSkipSimulationData  = 0;

}

SELM_Lagrangian::SELM_Lagrangian(class LAMMPS *lmps, class DriverSELM *fx) {

  const char *error_str_code = "SELM_Lagrangian.cpp";
  const char *error_str_func = "SELM_Lagrangian(narg,arg)";

  setGlobalRefs(lmps,fx);

  flagWriteSimulationData = 0;
  saveSkipSimulationData  = 0;
}

SELM_Lagrangian::~SELM_Lagrangian() {

}


void SELM_Lagrangian::setGlobalRefs(LAMMPS *lmps, DriverSELM *fx) {

  /* setup the LAMMPS related information */
  lammps  = lmps;
  driverSELM = fx;

}



