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
FixStyle(lb/momentum,FixLbMomentum);
// clang-format on
#else

#ifndef LMP_FIX_LB_MOMENTUM_H
#define LMP_FIX_LB_MOMENTUM_H

#include "fix.h"

namespace LAMMPS_NS {

class FixLbMomentum : public Fix {
 public:
  FixLbMomentum(class LAMMPS *, int, char **);
  int setmask() override;
  void init() override;
  void end_of_step() override;

 private:
  int linear;
  int xflag, yflag, zflag;
  double masstotal;
  int count;

  class FixLbFluid *fix_lb_fluid;
};

}    // namespace LAMMPS_NS

#endif
#endif
