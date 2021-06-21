/* ----------------------------------------------------------------------

 Table of types for SELM_Lagrangian derived classes.
   
 Paul J. Atzberger
 http://atzberger.org/
 
------------------------------------------------------------------------- */

#ifndef SELM_LAGRANGIAN_TYPES_H
#define SELM_LAGRANGIAN_TYPES_H

namespace LAMMPS_NS {

class SELM_Lagrangian_Types {

 public:
  SELM_Lagrangian_Types();
  virtual ~SELM_Lagrangian_Types();

  static const int   TYPE_NULL = 0;
  static const char *TYPE_STR_NULL;

  static const int   TYPE_CONTROLPTS_BASIC1 = 1;
  static const char *TYPE_STR_CONTROLPTS_BASIC1;

  static const int   TYPE_LAMMPS_ATOM_ANGLE_STYLE = 2;
  static const char *TYPE_STR_LAMMPS_ATOM_ANGLE_STYLE;

  static const int   TYPE_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE = 3;
  static const char *TYPE_STR_LAMMPS_HYBRID_CHARGE_ANGLE_STYLE;
  
  static const int   TYPE_LAMMPS_ATOM_STYLE_ELLIPSOID = 4;
  static const char *TYPE_STR_LAMMPS_ATOM_STYLE_ELLIPSOID;

};

}

#endif
