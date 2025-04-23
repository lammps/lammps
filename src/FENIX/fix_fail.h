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
FixStyle(fail,FixFail);
// clang-format on
#else

#ifndef LMP_FIX_FAIL_H
#define LMP_FIX_FAIL_H

#include "fix.h"

namespace LAMMPS_NS {

class FixFail : public Fix  {
 public:
  FixFail(class LAMMPS *, int, char**);

  void check_fail();

  int fail_rank = -1;
  int fail_timestep = -1;
  int fail_step = FixConst::END_OF_STEP;

  std::string fail_var = "";
  double fail_var_val;

  bool skip_double_failure = true;

  int setmask() override { return fail_step; }

  void initial_integrate(int) override { check_fail(); }
  void post_integrate() override { check_fail(); }
  void pre_exchange() override { check_fail(); }
  void pre_neighbor() override { check_fail(); }
  void post_neighbor() override { check_fail(); }
  void pre_force(int) override { check_fail(); }
  void pre_reverse(int, int) override { check_fail(); }
  void post_force(int) override { check_fail(); }
  void final_integrate() override { check_fail(); }
  void end_of_step() override { check_fail(); }
  void post_run() override { check_fail(); }
  void initial_integrate_respa(int, int, int) override { check_fail(); }
  void post_integrate_respa(int, int) override { check_fail(); }
  void pre_force_respa(int, int, int) override { check_fail(); }
  void post_force_respa(int, int, int) override { check_fail(); }
  void final_integrate_respa(int, int) override { check_fail(); }
  void min_pre_exchange() override { check_fail(); }
  void min_pre_neighbor() override { check_fail(); }
  void min_post_neighbor() override { check_fail(); }
  void min_pre_force(int) override { check_fail(); }
  void min_pre_reverse(int, int) override { check_fail(); }
  void min_post_force(int) override { check_fail(); }
};

}

#endif
#endif
