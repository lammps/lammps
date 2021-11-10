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
BondStyle(gaussian,BondGaussian);
// clang-format on
#else

#ifndef LMP_BOND_GAUSSIAN_H
#define LMP_BOND_GAUSSIAN_H

#include "bond.h"

namespace LAMMPS_NS {

class BondGaussian : public Bond {
 public:
  BondGaussian(class LAMMPS *);
  virtual ~BondGaussian();
  virtual void compute(int, int);
  virtual void coeff(int, char **);
  double equilibrium_distance(int);
  void write_restart(FILE *);
  virtual void read_restart(FILE *);
  void write_data(FILE *);
  double single(int, double, int, int, double &);

 protected:
  int *nterms;
  double *bond_temperature;
  double **alpha, **width, **r0;

  virtual void allocate();
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args for bond coefficients

Self-explanatory.  Check the input script or data file.

*/
