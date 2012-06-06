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

FixStyle(aveforce,FixAveForce)

#else

#ifndef LMP_FIX_AVEFORCE_H
#define LMP_FIX_AVEFORCE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixAveForce : public Fix {
 public:
  FixAveForce(class LAMMPS *, int, char **);
  ~FixAveForce();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double compute_vector(int);

 private:
  double xvalue,yvalue,zvalue;
  int varflag;
  char *xstr,*ystr,*zstr;
  char *idregion;
  int xvar,yvar,zvar,xstyle,ystyle,zstyle;
  int iregion;
  double foriginal_all[4];
  int nlevels_respa;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Region ID for fix aveforce does not exist

Self-explanatory.

E: Variable name for fix aveforce does not exist

Self-explanatory.

E: Variable for fix aveforce is invalid style

Only equal-style variables can be used.

*/
