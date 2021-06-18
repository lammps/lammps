/* ----------------------------------------------------------------------

  Stochastic Eulerian Lagrangian Methods

  Paul J. Atzberger
  http://atzberger.org/


------------------------------------------------------------------------- */
#include "SELM_Lagrangian_Types.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */
const char* SELM_Lagrangian_Types::TYPE_STR_NULL                              = "NULL";
const char* SELM_Lagrangian_Types::TYPE_STR_CONTROLPTS_BASIC1                 = "CONTROLPTS_BASIC1";
const char* SELM_Lagrangian_Types::TYPE_STR_LAMMPS_ATOM_ANGLE_STYLE           = "LAMMPS_ATOM_ANGLE_STYLE";
const char* SELM_Lagrangian_Types::TYPE_STR_LAMMPS_ATOM_STYLE_ELLIPSOID       = "LAMMPS_ATOM_STYLE_ELLIPSOID";
const char* SELM_Lagrangian_Types::TYPE_STR_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE  = "LAMMPS_HYBRID_CHARGE_ANGLE_STYLE";

SELM_Lagrangian_Types::SELM_Lagrangian_Types() {

}

SELM_Lagrangian_Types::~SELM_Lagrangian_Types() {

}
