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
FixStyle(press/berendsen,FixPressBerendsen);
// clang-format on
#else

#ifndef LMP_FIX_PRESS_BERENDSEN_H
#define LMP_FIX_PRESS_BERENDSEN_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPressBerendsen : public Fix {
 public:
  FixPressBerendsen(class LAMMPS *, int, char **);
  ~FixPressBerendsen() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void end_of_step() override;
  int modify_param(int, char **) override;

 protected:
  int dimension, which;
  double bulkmodulus;

  int pstyle, pcouple, allremap;
  int p_flag[3];    // 1 if control P on this dim, 0 if not
  double p_start[3], p_stop[3];
  double p_period[3], p_target[3];
  double p_current[3], dilation[3];
  double factor[3];
  int kspace_flag;            // 1 if KSpace invoked, 0 if not
  std::vector<Fix *> rfix;    // indices of rigid fixes

  char *id_temp, *id_press;
  class Compute *temperature, *pressure;
  int tflag, pflag;

  void couple();
  void remap();
};

}    // namespace LAMMPS_NS

#endif
#endif
