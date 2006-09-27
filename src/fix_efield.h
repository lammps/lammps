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

#ifndef FIX_EFIELD_H
#define FIX_EFIELD_H

#include "fix.h"

class FixEfield : public Fix {
 public:
  FixEfield(int, char **);
  ~FixEfield() {}
  int setmask();
  void init();
  void setup();
  void post_force(int);
  void post_force_respa(int, int, int);

 private:
  double ex,ey,ez;
  int nlevels_respa;
};

#endif
