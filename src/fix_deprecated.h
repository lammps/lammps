/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
// list all deprecated and removed fix styles here
FixStyle(DEPRECATED,FixDeprecated);
FixStyle(ave/spatial,FixDeprecated);
FixStyle(ave/spatial/sphere,FixDeprecated);
// clang-format on
#else

#ifndef LMP_FIX_DEPRECATED_H
#define LMP_FIX_DEPRECATED_H

#include "fix.h"

namespace LAMMPS_NS {

class FixDeprecated : public Fix {
 public:
  FixDeprecated(class LAMMPS *, int, char **);
  ~FixDeprecated() {}
  int setmask() { return 0; }
  void init() {}
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: This fix command has been removed from LAMMPS

UNDOCUMENTED

U: The fix ave/spatial command has been removed from LAMMPS

It has been replaced by the more flexible fix ave/chunk and compute
chunk/atom commands.  All the fix ave/spatial keywords and options are
available in those two newer commands.

*/
