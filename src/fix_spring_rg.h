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

#ifndef FIX_SPRING_RG_H
#define FIX_SPRING_RG_H

#include "fix.h"

class FixSpringRG : public Fix {
 public:
  FixSpringRG(int, char **);
  ~FixSpringRG() {}
  int setmask();
  void init();
  void setup();
  void post_force(int);
  void post_force_respa(int, int, int);

 private:
  int nlevels_respa,rg0_flag;
  double rg0,k,masstotal;
};

#endif
