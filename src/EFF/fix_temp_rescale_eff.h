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
FixStyle(temp/rescale/eff,FixTempRescaleEff);
// clang-format on
#else

#ifndef LMP_FIX_TEMP_RESCALE_EFF_H
#define LMP_FIX_TEMP_RESCALE_EFF_H

#include "fix.h"

namespace LAMMPS_NS {

class FixTempRescaleEff : public Fix {
 public:
  FixTempRescaleEff(class LAMMPS *, int, char **);
  ~FixTempRescaleEff() override;
  int setmask() override;
  void init() override;
  void end_of_step() override;
  int modify_param(int, char **) override;
  void reset_target(double) override;
  double compute_scalar() override;

 protected:
  int which;
  double t_start, t_stop, t_window;
  double fraction, energy, efactor;

  char *id_temp;
  class Compute *temperature;
  int tflag;
};

}    // namespace LAMMPS_NS

#endif
#endif
