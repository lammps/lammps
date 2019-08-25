/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(controller,FixController)

#else

#ifndef LMP_FIX_CONTROLLER_H
#define LMP_FIX_CONTROLLER_H

#include "fix.h"

namespace LAMMPS_NS {

class FixController : public Fix {
 public:
  FixController(class LAMMPS *, int, char **);
  ~FixController();
  int setmask();
  void init();
  void end_of_step();
  void reset_dt();
  double compute_vector(int);

 private:
  double kp,ki,kd,alpha,tau;
  double setpoint;
  int pvwhich,pvindex;
  char *pvID,*cvID;
  int firsttime;

  double control,err,olderr,deltaerr,sumerr;

  class Compute *pcompute;
  class Fix *pfix;
  int pvar,cvar;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

UNDOCUMENTED

E: Compute ID for fix controller does not exist

UNDOCUMENTED

E: Fix controller compute does not calculate a global scalar or vector

UNDOCUMENTED

E: Fix controller compute vector is accessed out-of-range

UNDOCUMENTED

E: Fix ID for fix controller does not exist

UNDOCUMENTED

E: Fix controller fix does not calculate a global scalar or vector

UNDOCUMENTED

E: Fix controller fix vector is accessed out-of-range

UNDOCUMENTED

E: Variable name for fix controller does not exist

UNDOCUMENTED

E: Fix controller variable is not equal-style variable

UNDOCUMENTED

E: Fix controller variable is not internal-style variable

UNDOCUMENTED

*/
