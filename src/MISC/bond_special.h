/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

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
  virtual ~BondSpecial();
  void init_style();
  void compute(int, int);
  void coeff(int, char **);
  double equilibrium_distance(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_data(FILE *);
  double single(int, double, int, int, double &);

 protected:
  double *factor_lj, *factor_coul;

  void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args for bond coefficients

Self-explanatory.  Check the input script or data file.

E: Invalid 1-2 setting for bond style special.

Bond style special must be used with zero factors for 1-2 special bonds.

E: Invalid 1-3 setting for bond style special.

Bond style special must be used with 1.0 factors for 1-3 special bonds or the
angle keyword set to yes.

E: Invalid 1-4 setting for bond style special.

Bond style special must be used with 1.0 factors for 1-4 special bonds or the
dihedral keyword set to yes.

E: Bond style special is not compatible with long range Coulombic interactions.

Self-explanatory.

*/
