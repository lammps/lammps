/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.

   Contributing Author: David Nicholson (MIT)
------------------------------------------------------------------------- */

#ifdef BOND_CLASS
// clang-format off
BondStyle(special,BondSpecial);
// clang-format on
#else

#ifndef LMP_BOND_SPECIAL_H
#define LMP_BOND_SPECIAL_H

#include "bond.h"

namespace LAMMPS_NS {

class BondSpecial : public Bond {
 public:
  BondSpecial(class LAMMPS *);
  ~BondSpecial() override;
  void init_style() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  double equilibrium_distance(int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_data(FILE *) override;
  double single(int, double, int, int, double &) override;

 protected:
  double *factor_lj, *factor_coul;

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
