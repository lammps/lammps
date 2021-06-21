/* ----------------------------------------------------------------------
 SELM Package.
 
 Paul J. Atzberger
 http://atzberger.org/
 
------------------------------------------------------------------------- */

#ifndef SELM_PACKAGE_H
#define SELM_PACKAGE_H

#include "mpi.h"
#include "lammps.h"
#include "error.h"
#include "mpi.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <stdio.h>
#include <stdlib.h>

using namespace std; /* ensures standard template library names used */

namespace LAMMPS_NS {

class SELM_Package {

 public:
  /* ================ Variables ================= */
  static LAMMPS *lammps;

  /* ================ Function prototypes ================= */
  SELM_Package();
  virtual ~SELM_Package(){};

  static void    setLAMMPS(LAMMPS *lammps_in);

  static void    packageError(const char *error_str_code, const char *error_str_func, stringstream &message);
  static void    packageError(const char *error_str_code, const char *error_str_func, const string &message);
  static void    packageWarning(const char *error_str_code, const char *error_str_func, stringstream &message);
  static void    packageWarning(const char *error_str_code, const char *error_str_func, const string &message);
};
}

#endif
