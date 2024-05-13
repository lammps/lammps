/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(wall/reflect/stochastic,FixWallReflectStochastic);
// clang-format on
#else

#ifndef LMP_FIX_WALL_REFLECT_STOCHASTIC_H
#define LMP_FIX_WALL_REFLECT_STOCHASTIC_H

#include "fix_wall_reflect.h"

namespace LAMMPS_NS {

class FixWallReflectStochastic : public FixWallReflect {
 public:
  FixWallReflectStochastic(class LAMMPS *, int, char **);
  ~FixWallReflectStochastic() override;

 private:
  int seedfix;
  double walltemp[6], wallvel[6][3], wallaccom[6][3];
  int rstyle;

  class RanMars *random;

  void wall_particle(int m, int which, double coord) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
