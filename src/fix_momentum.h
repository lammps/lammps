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

#ifndef FIX_MOMENTUM_H
#define FIX_MOMENTUM_H

#include "fix.h"

namespace LAMMPS_NS {

class FixMomentum : public Fix {
 public:
  FixMomentum(class LAMMPS *, int, char **);
  int setmask();
  void init();
  void end_of_step();

 private:
  int linear,angular;
  int xflag,yflag,zflag;
  double masstotal;
};

}

#endif
