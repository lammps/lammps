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

#ifndef FIX_SPRING_H
#define FIX_SPRING_H

#include "fix.h"

class FixSpring : public Fix {
 public:
  FixSpring(int, char **);
  ~FixSpring() {}
  int setmask();
  void init();
  void setup();
  void post_force(int);
  void post_force_respa(int, int, int);

 private:
  double xc,yc,zc,r0;
  double k_spring;
  int xflag,yflag,zflag;
  int styleflag;
  int igroup2,group2bit;
  double masstotal,masstotal2;
  int nlevels_respa;

  void spring_tether();
  void spring_couple();
};

#endif
