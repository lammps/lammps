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
FixStyle(wall/lj93,FixWallLJ93);
// clang-format on
#else

#ifndef LMP_FIX_WALL_LJ93_H
#define LMP_FIX_WALL_LJ93_H

#include "fix_wall.h"

namespace LAMMPS_NS {

class FixWallLJ93 : public FixWall {
 public:
  FixWallLJ93(class LAMMPS *, int, char **);
  void precompute(int) override;
  void wall_particle(int, int, double) override;

 protected:
  double coeff1[6], coeff2[6], coeff3[6], coeff4[6], offset[6];
};

}    // namespace LAMMPS_NS

#endif
#endif
