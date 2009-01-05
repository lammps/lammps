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

#ifndef FIX_WALL_REFLECT_H
#define FIX_WALL_REFLECT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixWallReflect : public Fix {
 public:
  FixWallReflect(class LAMMPS *, int, char **);
  int setmask();
  void post_integrate();

 private:
  int xloflag,xhiflag,yloflag,yhiflag,zloflag,zhiflag;
};

}

#endif
