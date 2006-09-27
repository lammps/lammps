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

#ifndef FIX_RECENTER_H
#define FIX_RECENTER_H

#include "fix.h"

class FixRecenter : public Fix {
 public:
  FixRecenter(int, char **);
  ~FixRecenter() {}
  int setmask();
  void init();
  void initial_integrate();

 private:
  int group2bit,scaleflag;
  int xflag,yflag,zflag;
  int xinitflag,yinitflag,zinitflag;
  double xcom,ycom,zcom,xinit,yinit,zinit,masstotal;
};

#endif
