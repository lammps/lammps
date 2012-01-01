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
  void post_force(int);
  void post_force_respa(int, int, int);
  double memory_usage();

 private:
  double ex,ey,ez;
  int varflag;
  char *xstr,*ystr,*zstr;
  int xvar,yvar,zvar,xstyle,ystyle,zstyle;
  int nlevels_respa;
  double qe2f;

  int maxatom;
  double **efield;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix efield requires atom attribute q

Self-explanatory.

E: Variable name for fix efield does not exist

Self-explanatory.

E: Variable for fix efield is invalid style

Only equal-style variables can be used.

*/
