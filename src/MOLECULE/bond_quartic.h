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
BondStyle(quartic,BondQuartic);
// clang-format on
#else

#ifndef LMP_BOND_QUARTIC_H
#define LMP_BOND_QUARTIC_H

#include "bond.h"

namespace LAMMPS_NS {

class BondQuartic : public Bond {
 public:
  BondQuartic(class LAMMPS *);
  virtual ~BondQuartic();
  virtual void compute(int, int);
  void coeff(int, char **);
  void init_style();
  double equilibrium_distance(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_data(FILE *);
  double single(int, double, int, int, double &);

 protected:
  double TWO_1_3;
  double *k, *b1, *b2, *rc, *u0;

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args for bond coefficients

Self-explanatory.  Check the input script or data file.

E: Pair style does not support bond_style quartic

The pair style does not have a single() function, so it can
not be invoked by bond_style quartic.

E: Bond style quartic cannot be used with 3,4-body interactions

No angle, dihedral, or improper styles can be defined when using
bond style quartic.

E: Bond style quartic cannot be used with atom style template

This bond style can change the bond topology which is not
allowed with this atom style.

E: Bond style quartic requires special_bonds = 1,1,1

This is a restriction of the current bond quartic implementation.

*/
