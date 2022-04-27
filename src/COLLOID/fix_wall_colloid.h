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
FixStyle(wall/colloid,FixWallColloid);
// clang-format on
#else

#ifndef LMP_FIX_WALL_COLLOID_H
#define LMP_FIX_WALL_COLLOID_H

#include "fix_wall.h"

namespace LAMMPS_NS {

class FixWallColloid : public FixWall {
 public:
  FixWallColloid(class LAMMPS *, int, char **);
  void init() override;
  void precompute(int) override;
  void wall_particle(int, int, double) override;

 private:
  double coeff1[6], coeff2[6], coeff3[6], coeff4[6];
};

}    // namespace LAMMPS_NS

#endif
#endif
