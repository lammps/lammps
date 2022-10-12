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
FixStyle(momentum,FixMomentum);
// clang-format on
#else

#ifndef LMP_FIX_MOMENTUM_H
#define LMP_FIX_MOMENTUM_H

#include "fix.h"

namespace LAMMPS_NS {

class FixMomentum : public Fix {
 public:
  FixMomentum(class LAMMPS *, int, char **);
  int setmask() override;
  void init() override;
  void end_of_step() override;

 protected:
  int linear, angular, rescale;
  int xflag, yflag, zflag;
  double masstotal;
};

}    // namespace LAMMPS_NS

#endif
#endif
