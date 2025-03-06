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
FixStyle(wall/flow,FixWallFlow);
// clang-format on
#else

#ifndef LMP_FIX_WALL_FLOW_H
#define LMP_FIX_WALL_FLOW_H

#include "fix.h"

namespace LAMMPS_NS {

class FixWallFlow : public Fix {
 public:
  enum FlowAxis { AX_X = 0, AX_Y = 1, AX_Z = 2 };

  FixWallFlow(class LAMMPS *, int, char **);
  ~FixWallFlow() override;
  int setmask() override;
  void init() override;
  void end_of_step() override;

  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;

  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;

 protected:
  FlowAxis flowax;
  double flowvel;
  double kT;
  std::vector<double> walls;

  int flowdir;
  int rndseed;
  class RanMars *random;
  int *current_segment;

  int compute_current_segment(double pos) const;
  void generate_velocity(int i);
};

}    // namespace LAMMPS_NS

#endif
#endif
