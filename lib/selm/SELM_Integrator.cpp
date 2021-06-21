/* ----------------------------------------------------------------------

  Stochastic Eulerian Lagrangian Methods

  Paul J. Atzberger
  http://atzberger.org/
 
------------------------------------------------------------------------- */

#include "driver_selm.h"
#include "SELM_Integrator.h"

#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"

#include <malloc.h>
#include <cstdlib>
#include "stdio.h"
#include "string.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

SELM_Integrator::SELM_Integrator() {

  const char *error_str_code = "SELM_Integrator.cpp";
  const char *error_str_func = "SELM_Integrator()";

}

SELM_Integrator::SELM_Integrator(int narg, char **arg) {

  const char *error_str_code = "SELM_Integrator.cpp";
  const char *error_str_func = "SELM_Integrator(narg,arg)";

}

SELM_Integrator::~SELM_Integrator() {

}

void SELM_Integrator::setGlobalRefs(LAMMPS *lmps, DriverSELM *fx) {

  /* setup the LAMMPS related information */
  lammps  = lmps;
  driver_selm = fx;

}

void SELM_Integrator::post_force(int vflag) {
  // perform operations after computing forces
}

void SELM_Integrator::setup_LAMMPS(int vflag) {
  // perform operations after computing forces
}


