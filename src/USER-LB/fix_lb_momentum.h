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

FixStyle(lb/momentum,FixLbMomentum)

#else

#ifndef LMP_FIX_LB_MOMENTUM_H
#define LMP_FIX_LB_MOMENTUM_H

#include "fix.h"

namespace LAMMPS_NS {

class FixLbMomentum : public Fix {
 public:
  FixLbMomentum(class LAMMPS *, int, char **);
  int setmask();
  void init();
  void end_of_step();

 private:
  int linear;
  int xflag,yflag,zflag;
  double masstotal;

  class FixLbFluid *fix_lb_fluid;
};

}

#endif
#endif
