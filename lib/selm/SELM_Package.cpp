/* ----------------------------------------------------------------------

  Stochastic Eulerian Lagrangian Methods

  Paul J. Atzberger
  http://atzberger.org/
 

------------------------------------------------------------------------- */

#include "SELM_Package.h"
#include <stdio.h>
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

SELM_Package::~SELM_Package() {

}


LAMMPS *SELM_Package::getLAMMPS() {
  return lammps;
}


void SELM_Package::setLAMMPS(LAMMPS *lammps_in) {
  lammps = lammps_in;
}

void SELM_Package::packageError(const char *error_str_code, const char *error_str_func, stringstream &message) {
  packageError(lammps, error_str_code, error_str_func, message.str().c_str());
}

void SELM_Package::packageError(LAMMPS *lammps_in, const char *error_str_code, const char *error_str_func, stringstream &message) {
  packageError(lammps_in, error_str_code, error_str_func, message.str().c_str());
}

void SELM_Package::packageError(const char *error_str_code, const char *error_str_func, string &message) {
  packageError(lammps, error_str_code, error_str_func, message.c_str());
}

void SELM_Package::packageError(LAMMPS *lammps_in, const char *error_str_code, const char *error_str_func, string &message) {
  packageError(lammps_in, error_str_code, error_str_func, message.c_str());
}

void SELM_Package::packageError(const char *error_str_code, const char *error_str_func, const char *message) {
  packageError(lammps, error_str_code, error_str_func, message);
}

void SELM_Package::packageError(LAMMPS *lammps_in, const char *error_str_code, const char *error_str_func, const char *message) {

  int code = 1;

  stringstream outputStream;
  string       outputStr;

  /* generate the error message string */
  outputStream << error_str_code << " : " << error_str_func << endl;
  outputStream << message << endl;

  outputStr = outputStream.str();

  /* output the error message */
  /*
  if (lammps_in != NULL) {
    lammps_in->error->all(outputStr.c_str());
  } else {
    cerr << "ERROR: " << outputStr.c_str() << endl;
  }
  */
  cerr << "ERROR: " << outputStr.c_str() << endl;

  /* finalize any processes and exit the program */
  MPI_Finalize();  /* be sure to close all MPI processes */

  std::exit(code);

}

void SELM_Package::packageWarning(const char *error_str_code, const char *error_str_func, stringstream &message) {
  packageWarning(lammps, error_str_code, error_str_func, message.str().c_str());
}

void SELM_Package::packageWarning(const char *error_str_code, const char *error_str_func, string &message) {
  packageWarning(lammps, error_str_code, error_str_func, message.c_str());
}

void SELM_Package::packageWarning(const char *error_str_code, const char *error_str_func, const char *message) {
  packageError(lammps, error_str_code, error_str_func, message);
}

void SELM_Package::packageWarning(LAMMPS *lammps_in, const char *error_str_code, const char *error_str_func, const char *message) {
  int code = 1;
  int logFlag = 1;

  stringstream outputStream;
  string       outputStr;

  /* generate the error message string */
  outputStream << error_str_code << " : " << error_str_func << endl;
  outputStream << message << endl;

  outputStr = outputStream.str();

  /* output the error message */
  /*
  if (lammps_in != NULL) {
    logFlag = 1;
    lammps_in->error->warning(outputStr.c_str(), logFlag);
  } else {
    cerr << "WARNING: " << outputStr.c_str() << endl;
  }
  */
  cerr << "WARNING: " << outputStr.c_str() << endl;

  /* finalize any processes and exit the program */
  //MPI_Finalize();  /* be sure to close all MPI processes */
  //exit(code);

}


