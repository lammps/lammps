/* ----------------------------------------------------------------------
 Table of types for SELM_Interaction derived classes.
   
 Paul J. Atzberger
 http://atzberger.org/
    
------------------------------------------------------------------------- */

#ifndef SELM_INTERACTION_TYPES_H
#define SELM_INTERACTION_TYPES_H

namespace LAMMPS_NS {

class SELM_Interaction_Types {

 public:
  SELM_Interaction_Types();
  virtual ~SELM_Interaction_Types();

  static const int   TYPE_NULL = 0;
  static const char *TYPE_STR_NULL;

  static const int   TYPE_CUSTOM1 = 1;
  static const char *TYPE_STR_CUSTOM1;

  static const int   TYPE_SKIPDATA = 2;
  static const char *TYPE_STR_SKIPDATA;

};

}

#endif
