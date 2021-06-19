/* ----------------------------------------------------------------------
   Table of types for SELM_CouplingOperators derived classes.
   
   Paul J. Atzberger
   http://atzberger.org/
 
------------------------------------------------------------------------- */

#ifndef SELM_COUPLINGOPERATOR_TYPES_H
#define SELM_COUPLINGOPERATOR_TYPES_H

namespace LAMMPS_NS {

class SELM_CouplingOperator_Types {

 public:
  SELM_CouplingOperator_Types();
  virtual ~SELM_CouplingOperator_Types();

  static const int   TYPE_NULL;
  static const char *TYPE_STR_NULL;

  static const int   TYPE_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1;
  static const char *TYPE_STR_LAMMPS_SHEAR_UNIFORM1_FFTW3_TABLE1;

};

}

#endif
