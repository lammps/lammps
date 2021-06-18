/* ----------------------------------------------------------------------

  Stochastic Eulerian Lagrangian Methods

  Paul J. Atzberger
  http://atzberger.org/
 
------------------------------------------------------------------------- */

#include "driver_selm.h"
#include "lammps.h"

#include "SELM_CouplingOperator.h"

#include "atom.h"
#include "error.h"
#include "force.h"
#include "update.h"
#include "respa.h"

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <malloc.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

SELM_CouplingOperator::SELM_CouplingOperator() {

  const char *error_str_code = "SELM_CouplingOperator.cpp";
  const char *error_str_func = "SELM_CouplingOperator()";

}

SELM_CouplingOperator::SELM_CouplingOperator(int narg, char **arg) {

  const char *error_str_code = "SELM_CouplingOperator.cpp";
  const char *error_str_func = "SELM_CouplingOperator(narg,arg)";

}

SELM_CouplingOperator::~SELM_CouplingOperator() {

}

void SELM_CouplingOperator::setGlobalRefs(LAMMPS *lmps, DriverSELM *fx) {

  lammps  = lmps;
  driverSELM = fx;

}

