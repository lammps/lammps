/* ----------------------------------------------------------------------

  Stochastic Eulerian Lagrangian Methods

  Paul J. Atzberger
  http://atzberger.org/
 

------------------------------------------------------------------------- */
#include "SELM_Eulerian_Types.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */
const char* SELM_Eulerian_Types::TYPE_STR_NULL                        = "NULL";
const char* SELM_Eulerian_Types::TYPE_STR_FLUID_SHEAR_UNIFORM1_FFTW3  = "FLUID_SHEAR_UNIFORM1_FFTW3";
const char* SELM_Eulerian_Types::TYPE_STR_LAMMPS_SHEAR_UNIFORM1_FFTW3 = "LAMMPS_SHEAR_UNIFORM1_FFTW3";
const char* SELM_Eulerian_Types::TYPE_STR_StaggeredGrid1              = "StaggeredGrid1";
const char* SELM_Eulerian_Types::TYPE_STR_Uniform1_Periodic           = "Uniform1_Periodic";

SELM_Eulerian_Types::SELM_Eulerian_Types() {

}

SELM_Eulerian_Types::~SELM_Eulerian_Types() {

}
