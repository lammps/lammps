/* ----------------------------------------------------------------------

  Stochastic Eulerian Lagrangian Methods

  Paul J. Atzberger
  http://atzberger.org/
 

------------------------------------------------------------------------- */
#include "SELM_CouplingOperator_Types.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */
const int   SELM_CouplingOperator_Types::TYPE_NULL                                   = 0;
const char* SELM_CouplingOperator_Types::TYPE_STR_NULL                               = "NULL";

const int   SELM_CouplingOperator_Types::TYPE_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1     = 1;
const char* SELM_CouplingOperator_Types::TYPE_STR_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1 = "LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1";

SELM_CouplingOperator_Types::SELM_CouplingOperator_Types() {

}

SELM_CouplingOperator_Types::~SELM_CouplingOperator_Types() {

}
