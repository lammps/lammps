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
BondStyle(fene/expand, BondFENEExpand);
// clang-format on
#else

#ifndef LMP_BOND_FENE_EXPAND_H
#define LMP_BOND_FENE_EXPAND_H

#include "bond.h"

namespace LAMMPS_NS {

class BondFENEExpand : public Bond {
 public:
  BondFENEExpand(class LAMMPS *_lmp) : Bond(_lmp) {}
  ~BondFENEExpand() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;
  double equilibrium_distance(int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_data(FILE *) override;
  double single(int, double, int, int, double &) override;

 protected:
  double *k, *r0, *epsilon, *sigma, *shift;

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
