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

#ifdef KSPACE_CLASS
// clang-format off
KSpaceStyle(DEPRECATED,KSpaceDeprecated);
// clang-format on
#else

#ifndef LMP_KSPACE_DEPRECATED_H
#define LMP_KSPACE_DEPRECATED_H

#include "kspace.h"

namespace LAMMPS_NS {

class KSpaceDeprecated : public KSpace {
 public:
  KSpaceDeprecated(class LAMMPS *lmp) : KSpace(lmp) {}
  virtual ~KSpaceDeprecated() {}

  virtual void init() {}
  virtual void settings(int, char **);
  virtual void setup() {}
  virtual void compute(int, int) {}
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

*/
