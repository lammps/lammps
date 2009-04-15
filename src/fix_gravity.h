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

#ifndef FIX_GRAVITY_H
#define FIX_GRAVITY_H

#include "fix.h"

namespace LAMMPS_NS {

class FixGravity : public Fix {
  friend class FixPour;

 public:
  FixGravity(class LAMMPS *, int, char **);
  int setmask();
  void init();
  void setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);

 private:
  int style,time_initial;
  double magnitude,dt;
  double phi,theta,phigrad,thetagrad;
  double xdir,ydir,zdir;
  double xgrav,ygrav,zgrav,xacc,yacc,zacc;
  double degree2rad;
  int nlevels_respa;
};

}

#endif
