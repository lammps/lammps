/* ----------------------------------------------------------------------

  Stochastic Eulerian Lagrangian Methods

  Paul J. Atzberger
  http://atzberger.org/
 

------------------------------------------------------------------------- */
#include "SELM_Interaction.h"
#include "driver_selm.h"

#include "lammps.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "respa.h"
#include "error.h"

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <malloc.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

SELM_Interaction::SELM_Interaction() {

  const char *error_str_code = "SELM_Interaction.cpp";
  const char *error_str_func = "SELM_Interaction()";

}

SELM_Interaction::SELM_Interaction(int narg, char **arg) {

  const char *error_str_code = "SELM_Interaction.cpp";
  const char *error_str_func = "SELM_Interaction(narg,arg)";

}

SELM_Interaction::~SELM_Interaction() {

}

void SELM_Interaction::setGlobalRefs(LAMMPS *lmps, DriverSELM *fx) {

  /* setup the LAMMPS related information */
  lammps  = lmps;
  driverSELM = fx;

}

