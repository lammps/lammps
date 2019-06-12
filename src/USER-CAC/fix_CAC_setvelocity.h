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

FixStyle(CAC/setvelocity,FixCAC_Set_Velocity)

#else

#ifndef LMP_FIX_CAC_SETVELOCITY_H
#define LMP_FIX_CAC_SETVELOCITY_H

#include "fix.h"

namespace LAMMPS_NS {

class FixCAC_Set_Velocity : public Fix {
 public:
	 FixCAC_Set_Velocity(class LAMMPS *, int, char **);
  virtual ~FixCAC_Set_Velocity();
  int setmask();
  virtual void init();
  void setup(int);
  void min_setup(int);
  virtual void initial_integrate(int vflag);
  void min_post_force(int);
  double compute_vector(int);

  double memory_usage();

 protected:
  double xvalue,yvalue,zvalue;
  int varflag,iregion;
  char *xstr,*ystr,*zstr;
  char *idregion;
  int xvar,yvar,zvar,xstyle,ystyle,zstyle;
  double voriginal[3],voriginal_all[3];
  int force_flag;

  int maxatom;
  double **svelocity;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Region ID for fix setforce does not exist

Self-explanatory.

E: Variable name for fix setforce does not exist

Self-explanatory.

E: Variable for fix setforce is invalid style

Only equal-style variables can be used.

E: Cannot use non-zero forces in an energy minimization

Fix setforce cannot be used in this manner.  Use fix addforce
instead.

*/
