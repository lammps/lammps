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
  virtual ~BondDeprecated() {}

  virtual void compute(int, int) {}
  virtual void settings(int, char **);
  virtual void coeff(int, char **) {}
  virtual double equilibrium_distance(int) { return 0.0; }
  virtual void write_restart(FILE *) {}
  virtual void read_restart(FILE *) {}
  virtual double single(int, double, int, int, double &) { return 0.0; }
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

*/
