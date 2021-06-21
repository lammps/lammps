/* ----------------------------------------------------------------------
   Table of types for SELM_Integrator derived classes.
   
   Paul J. Atzberger
   http://atzberger.org/
 
------------------------------------------------------------------------- */

#ifndef SELM_INTEGRATOR_TYPES_H
#define SELM_INTEGRATOR_TYPES_H

namespace LAMMPS_NS {

class SELM_Integrator_Types {

 public:
  SELM_Integrator_Types();
  virtual ~SELM_Integrator_Types();

  static const int   TYPE_NULL = 0;
  static const char       *TYPE_STR_NULL;

  static const int   TYPE_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3 = 2;
  static const char       *TYPE_STR_LAMMPS_SHEAR_QUASI_STEADY1_FFTW3;

};

}

#endif
