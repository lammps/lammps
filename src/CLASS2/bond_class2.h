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
BondStyle(class2,BondClass2);
// clang-format on
#else

#ifndef LMP_BOND_CLASS2_H
#define LMP_BOND_CLASS2_H

#include "bond.h"

namespace LAMMPS_NS {

class BondClass2 : public Bond {
 public:
  BondClass2(class LAMMPS *);
  ~BondClass2() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  double equilibrium_distance(int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_data(FILE *) override;
  double single(int, double, int, int, double &) override;
  void *extract(const char *, int &) override;

 protected:
  double *r0, *k2, *k3, *k4;

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
