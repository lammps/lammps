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

FixStyle(bfield,FixBfield)

#else

#ifndef LMP_FIX_BFIELD_H
#define LMP_FIX_BFIELD_H

#include "fix.h"
#include "region.h"

namespace LAMMPS_NS {

class FixBfield : public Fix {
 public:
  FixBfield(class LAMMPS *, int, char **);
  ~FixBfield();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void initial_integrate(int);
  void post_integrate();
  void post_force(int);
  double memory_usage();
  double compute_scalar();

 private:
  Region *region;
  char *idregion;
  int varflag,iregion;
  char *xstr,*ystr,*zstr,*estr;
  int xvar,yvar,zvar,evar,xstyle,ystyle,zstyle,estyle;
  int maxatom;
  double dtf;
  double **v0, **fb;
  double B[3];
  double omega[3];
  //double qBm2f;
  int bmuflag;
  int qflag;
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

E: Fix bfield requires atom attribute q 

The atom style defined does not have this attribute.

*/
