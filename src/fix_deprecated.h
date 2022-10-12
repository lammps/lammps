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
FixStyle(lb/pc,FixDeprecated);
FixStyle(lb/rigid/pc/sphere,FixDeprecated);
FixStyle(client/md,FixDeprecated);
// clang-format on
#else

#ifndef LMP_FIX_DEPRECATED_H
#define LMP_FIX_DEPRECATED_H

#include "fix.h"

namespace LAMMPS_NS {

class FixDeprecated : public Fix {
 public:
  FixDeprecated(class LAMMPS *, int, char **);

  int setmask() override { return 0; }
  void init() override {}
};

}    // namespace LAMMPS_NS

#endif
#endif
