/* ----------------------------------------------------------------------
   Table of types for SELM_Eulerian derived classes.
   
   Paul J. Atzberger
   http://atzberger.org/
    
------------------------------------------------------------------------- */

#ifndef SELM_EULERIAN_TYPES_H
#define SELM_EULERIAN_TYPES_H

namespace LAMMPS_NS {

class SELM_Eulerian_Types {

 public:
  SELM_Eulerian_Types();
  virtual ~SELM_Eulerian_Types();

  static const int        TYPE_NULL = 0;
  static const char      *TYPE_STR_NULL;

  static const int        TYPE_FLUID_SHEAR_UNIFORM1_FFTW3 = 1;
  static const char      *TYPE_STR_FLUID_SHEAR_UNIFORM1_FFTW3;

  static const int        TYPE_LAMMPS_SHEAR_UNIFORM1_FFTW3 = 2;
  static const char      *TYPE_STR_LAMMPS_SHEAR_UNIFORM1_FFTW3;

  static const int        TYPE_StaggeredGrid1 = 3;
  static const char      *TYPE_STR_StaggeredGrid1;

  static const int        TYPE_Uniform1_Periodic = 4;
  static const char      *TYPE_STR_Uniform1_Periodic;

};

}

#endif
