/* -*- c++ -*- ----------------------------------------------------------
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

FixStyle(ave/spatial,FixAveSpatial)

#else

#ifndef LMP_FIX_AVE_SPATIAL_H
#define LMP_FIX_AVE_SPATIAL_H

#include <stdio.h>
#include "fix.h"

namespace LAMMPS_NS {

class FixAveSpatial : public Fix {
 public:
  FixAveSpatial(class LAMMPS *, int, char **);
  ~FixAveSpatial() {}
  int setmask() {}
  void init() {}
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: The fix ave/spatial command has been removed from LAMMPS

It has been replaced by the more flexible fix ave/chunk and compute
chunk/atom commands.  All the fix ave/spatial keywords and options are
available in those two newer commands.

*/
