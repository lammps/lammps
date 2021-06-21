/* ----------------------------------------------------------------------

  Stochastic Eulerian Lagrangian Methods : Integrator Types

  Paul J. Atzberger
  http://atzberger.org/
 

------------------------------------------------------------------------- */
#include "SELM_Integrator_Types.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */
const char* SELM_Integrator_Types::TYPE_STR_NULL                              = "NULL";
const char* SELM_Integrator_Types::TYPE_STR_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3  = "LAMMPS_SHEAR_QUASI_STEADY1_FFTW3";

SELM_Integrator_Types::SELM_Integrator_Types() {

}

SELM_Integrator_Types::~SELM_Integrator_Types() {

}
