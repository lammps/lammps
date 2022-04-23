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
FixStyle(dt/reset,FixDtReset);
// clang-format on
#else

#ifndef LMP_FIX_DT_RESET_H
#define LMP_FIX_DT_RESET_H

#include "fix.h"

namespace LAMMPS_NS {

class FixDtReset : public Fix {
 public:
  FixDtReset(class LAMMPS *, int, char **);

  int setmask() override;
  void init() override;
  void setup(int) override;
  void end_of_step() override;
  double compute_scalar() override;

 private:
  bigint laststep;
  int minbound, maxbound;
  double tmin, tmax, xmax, emax;
  double ftm2v, mvv2e;
  double dt, t_laststep;
  int respaflag;
};

}    // namespace LAMMPS_NS

#endif
#endif
