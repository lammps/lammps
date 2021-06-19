/* ----------------------------------------------------------------------

  Stochastic Eulerian Lagrangian Methods

  Paul J. Atzberger
  http://atzberger.org/
 

------------------------------------------------------------------------- */

#include "SELM_Package.h"

#include "comm.h"

#include <cstdio>
#include <cstdlib>

using namespace std; /* ensures standard template library in namespace */
using namespace LAMMPS_NS;


/* ---------------------------------------------------------------------- */
/* Constant valued strings                                                */
/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

/* initialize static variables */
LAMMPS *SELM_Package::lammps = NULL;

SELM_Package::SELM_Package() {

  const char *error_str_code = "SELM_Package.cpp";
  const char *error_str_func = "SELM_Package()";

  lammps = NULL;
}

void SELM_Package::setLAMMPS(LAMMPS *lammps_in) {
  lammps = lammps_in;
}

void SELM_Package::packageError(const char *error_str_code, const char *error_str_func, stringstream &message) {
  packageError(error_str_code, error_str_func, message.str());
}

void SELM_Package::packageError(const char *error_str_code, const char *error_str_func, const string &message) {
  int code = 1;

  lammps->error->one(error_str_code,code,"{}: {}",error_str_func,message);
}

void SELM_Package::packageWarning(const char *error_str_code, const char *error_str_func, stringstream &message) {
  packageWarning(error_str_code, error_str_func, message.str());
}

void SELM_Package::packageWarning(const char *error_str_code, const char *error_str_func, const string &message) {
  int code = 1;
  if (lammps->comm->me == 0)
    lammps->error->warning(error_str_code,code,"{}: {}",error_str_func,message);
}
