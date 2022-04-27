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
FixStyle(temp/berendsen,FixTempBerendsen);
// clang-format on
#else

#ifndef LMP_FIX_TEMP_BERENDSEN_H
#define LMP_FIX_TEMP_BERENDSEN_H

#include "fix.h"

namespace LAMMPS_NS {

class FixTempBerendsen : public Fix {
 public:
  FixTempBerendsen(class LAMMPS *, int, char **);
  ~FixTempBerendsen() override;
  int setmask() override;
  void init() override;
  void end_of_step() override;
  int modify_param(int, char **) override;
  void reset_target(double) override;
  double compute_scalar() override;
  void write_restart(FILE *) override;
  void restart(char *buf) override;
  void *extract(const char *, int &) override;

 private:
  int which;
  double t_start, t_stop, t_period, t_target;
  double energy;
  int tstyle, tvar;
  char *tstr;

  char *id_temp;
  class Compute *temperature;
  int tflag;
};

}    // namespace LAMMPS_NS

#endif
#endif
