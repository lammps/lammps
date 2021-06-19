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

  /* ================ Constants ================= */

  /* ================ Data structure type definitions ================= */

  /* ================ Variables ================= */
  static LAMMPS *lammps;

  /* ================ Function prototypes ================= */
  SELM_Package();
  virtual ~SELM_Package();

  static void    setLAMMPS(LAMMPS *lammps_in);
  static LAMMPS *getLAMMPS();

  static void    packageError(const char *error_str_code, const char *error_str_func, stringstream &message);
  static void    packageError(LAMMPS *lammps_in, const char *error_str_code, const char *error_str_func, stringstream &message);

  static void    packageError(const char *error_str_code, const char *error_str_func, string &message);
  static void    packageError(LAMMPS *lammps_in, const char *error_str_code, const char *error_str_func, string &message);

  static void    packageError(const char *error_str_code, const char *error_str_func, const char *message);
  static void    packageError(LAMMPS *lammps_in, const char *error_str_code, const char *error_str_func, const char *message);

  static void    packageWarning(const char *error_str_code, const char *error_str_func, stringstream &message);
  static void    packageWarning(const char *error_str_code, const char *error_str_func, string &message);
  static void    packageWarning(const char *error_str_code, const char *error_str_func, const char *message);
  static void    packageWarning(LAMMPS *lammps_in, const char *error_str_code, const char *error_str_func, const char *message);

};

}

#endif
