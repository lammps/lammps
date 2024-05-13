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
BondStyle(harmonic/restrain,BondHarmonicRestrain);
// clang-format on
#else

#ifndef LMP_BOND_HARMONIC_RESTRAIN_H
#define LMP_BOND_HARMONIC_RESTRAIN_H

#include "bond.h"

namespace LAMMPS_NS {

class BondHarmonicRestrain : public Bond {
 public:
  BondHarmonicRestrain(class LAMMPS *);
  ~BondHarmonicRestrain() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;
  double equilibrium_distance(int) override { return -1.0; };  // return invalid value since not applicable
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_data(FILE *) override;
  double single(int, double, int, int, double &) override;
  void *extract(const char *, int &) override;

 protected:
  double *k;
  bigint natoms;
  class FixStoreAtom *initial;

  virtual void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif
