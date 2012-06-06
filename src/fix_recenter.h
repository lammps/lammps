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

FixStyle(recenter,FixRecenter)

#else

#ifndef LMP_FIX_RECENTER_H
#define LMP_FIX_RECENTER_H

#include "fix.h"

namespace LAMMPS_NS {

class FixRecenter : public Fix {
 public:
  FixRecenter(class LAMMPS *, int, char **);
  int setmask();
  void init();
  void initial_integrate(int);

 private:
  int group2bit,scaleflag;
  int xflag,yflag,zflag;
  int xinitflag,yinitflag,zinitflag;
  double xcom,ycom,zcom,xinit,yinit,zinit,masstotal;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Could not find fix recenter group ID

A group ID used in the fix recenter command does not exist.

E: Use of fix recenter with undefined lattice

Must use lattice command with fix recenter command if units option is
set to lattice.

E: Fix recenter group has no atoms

Self-explanatory.

W: Fix recenter should come after all other integration fixes

Other fixes may change the position of the center-of-mass, so
fix recenter should come last.

*/
