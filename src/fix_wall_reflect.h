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

FixStyle(wall/reflect,FixWallReflect)

#else

#ifndef LMP_FIX_WALL_REFLECT_H
#define LMP_FIX_WALL_REFLECT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixWallReflect : public Fix {
 public:
  FixWallReflect(class LAMMPS *, int, char **);
  ~FixWallReflect();
  int setmask();
  void init();
  void post_integrate();

 private:
  int nwall;
  int wallwhich[6],wallstyle[6];
  double coord0[6];
  char *varstr[6];
  int varindex[6];
  int varflag;
};

}

#endif
#endif
