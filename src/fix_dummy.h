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
FixStyle(DUMMY,FixDummy);
// clang-format on
#else

#ifndef LMP_FIX_DUMMY_H
#define LMP_FIX_DUMMY_H

#include "fix.h"

namespace LAMMPS_NS {

class FixDummy : public Fix {
 public:
  FixDummy(class LAMMPS *, int, char **);

  int setmask() override;

 protected:
  int initial_integrate_flag, final_integrate_flag;
  int pre_exchange_flag, pre_neighbor_flag;
  int pre_force_flag, post_force_flag;
  int end_of_step_flag;
};

}    // namespace LAMMPS_NS

#endif
#endif
