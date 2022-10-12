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
FixStyle(bond/create/angle,FixBondCreateAngle);
// clang-format on
#else

#ifndef LMP_FIX_BOND_CREATE_ANGLE_H
#define LMP_FIX_BOND_CREATE_ANGLE_H

#include "fix_bond_create.h"

namespace LAMMPS_NS {

class FixBondCreateAngle : public FixBondCreate {
 public:
  FixBondCreateAngle(LAMMPS *_lmp, int narg, char **arg) : FixBondCreate(_lmp, narg, arg) {}

 private:
  int constrain(int, int, double, double) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
