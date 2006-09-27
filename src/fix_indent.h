/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef FIX_INDENT_H
#define FIX_INDENT_H

#include "fix.h"

class FixIndent : public Fix {
 public:
  FixIndent(int, char **);
  ~FixIndent() {}
  int setmask();
  void init();
  void setup();
  void min_setup();
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  int thermo_fields(int, int *, char **);
  int thermo_compute(double *);

 private:
  int ntimestep_initial,istyle,scaleflag,radflag,thermo_flag,eflag_on;
  double k,k3,eng,etotal;
  double x0,y0,z0,r0_stop,r0_start;
  int cdim;
  double c1,c2;
  double vx,vy,vz;
  int nlevels_respa;

  void options(int, char **);
};

#endif
