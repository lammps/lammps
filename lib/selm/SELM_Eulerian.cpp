/* ----------------------------------------------------------------------

  Stochastic Eulerian Lagrangian Methods

  Paul J. Atzberger
  http://atzberger.org/
 

------------------------------------------------------------------------- */

#include <cstdlib>
#include "stdio.h"
#include "string.h"
#include "SELM_Eulerian.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"

#include <malloc.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

SELM_Eulerian::SELM_Eulerian() {

  const char *error_str_code = "SELM_Eulerian.cpp";
  const char *error_str_func = "SELM_Eulerian()";

}

SELM_Eulerian::SELM_Eulerian(int narg, char **arg) {

  const char *error_str_code = "SELM_Eulerian.cpp";
  const char *error_str_func = "SELM_Eulerian(narg,arg)";

}

SELM_Eulerian::~SELM_Eulerian() {

}

void SELM_Eulerian::setGlobalRefs(LAMMPS *lmps, DriverSELM *fx) {

  lammps  = lmps;
  driverSELM = fx;

}


