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

#ifdef BOND_CLASS
// clang-format off
BondStyle(DEPRECATED,BondDeprecated);
// clang-format on
#else

#ifndef LMP_BOND_DEPRECATED_H
#define LMP_BOND_DEPRECATED_H

#include "bond.h"

namespace LAMMPS_NS {

class BondDeprecated : public Bond {
 public:
  BondDeprecated(class LAMMPS *lmp) : Bond(lmp) {}

  void compute(int, int) override {}
  void settings(int, char **) override;
  void coeff(int, char **) override {}
  double equilibrium_distance(int) override { return 0.0; }
  void write_restart(FILE *) override {}
  void read_restart(FILE *) override {}
  double single(int, double, int, int, double &) override { return 0.0; }
};

}    // namespace LAMMPS_NS

#endif
#endif
