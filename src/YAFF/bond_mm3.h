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

#ifdef BOND_CLASS
// clang-format off
BondStyle(mm3,BondMM3);
// clang-format on
#else

#ifndef LMP_BOND_MM3_H
#define LMP_BOND_MM3_H

#include "bond.h"

namespace LAMMPS_NS {

class BondMM3 : public Bond {
 public:
  BondMM3(class LAMMPS *);
  ~BondMM3() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  double equilibrium_distance(int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_data(FILE *) override;
  double single(int, double, int, int, double &) override;

 protected:
  double *r0, *k2;

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
