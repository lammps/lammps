/* ----------------------------------------------------------------------
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

FixStyle(efield,FixEfield)

#else

#ifndef LMP_FIX_EFIELD_H
#define LMP_FIX_EFIELD_H

#include "fix.h"

namespace LAMMPS_NS {

class FixEfield : public Fix {
 public:
  FixEfield(class LAMMPS *, int, char **);
  ~FixEfield();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double memory_usage();
  double compute_scalar();
  double compute_vector(int);

 private:
  double ex,ey,ez;
  int varflag;
  char *xstr,*ystr,*zstr,*estr;
  int xvar,yvar,zvar,evar,xstyle,ystyle,zstyle,estyle;
  int nlevels_respa;
  double qe2f;
  double fdotx;
  int qflag,muflag;

  int maxatom;
  double **efield;

  int force_flag;
  double fsum[4],fsum_all[4];
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix efield requires atom attribute q or mu

Self-explanatory.

E: Variable name for fix efield does not exist

Self-explanatory.

E: Variable for fix efield is invalid style

Only equal-style or atom-style variables can be used.

E: Fix efield with dipoles cannot use atom-style variables

This feature is not yet supported.

W: The minimizer does not re-orient dipoles when using fix efield

Self-explanatory.

E: Cannot use variable energy with constant efield in fix efield

Self-explanatory.

E: Must use variable energy with fix efield

One or more variables are defined for fix efield, which require
variable energy when using the minimizer.

*/
