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

#ifndef FIX_GRAVITY_H
#define FIX_GRAVITY_H

#include "fix.h"

class FixGravity : public Fix {
  friend class FixPour;

 public:
  FixGravity(int, char **);
  ~FixGravity() {}
  int setmask();
  void init();
  void setup();
  void post_force(int);

 private:
  double phi,theta,phigrad,thetagrad;
  int dynamic,time_initial;
  int granular;                      // 0 if non-granular, 1 if granular
  double magnitude,xdir,ydir,zdir;
  double dt;
  double xgrav,ygrav,zgrav;
  double degree2rad;
};

#endif
