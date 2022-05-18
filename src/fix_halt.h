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
FixStyle(halt,FixHalt);
// clang-format on
#else

#ifndef LMP_FIX_HALT_H
#define LMP_FIX_HALT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixHalt : public Fix {
 public:
  FixHalt(class LAMMPS *, int, char **);
  ~FixHalt() override;
  int setmask() override;
  void init() override;
  void end_of_step() override;
  void min_post_force(int) override;
  void post_run() override;

 private:
  int attribute, operation, eflag, msgflag, ivar;
  bigint nextstep, thisstep;
  double value, tratio;
  char *idvar;
  char *dlimit_path;

  double bondmax();
  double tlimit();
  double diskfree();
};

}    // namespace LAMMPS_NS

#endif
#endif
