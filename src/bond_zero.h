/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef BOND_CLASS

BondStyle(zero,BondZero)

#else

#ifndef LMP_BOND_ZERO_H
#define LMP_BOND_ZERO_H

#include "bond.h"

namespace LAMMPS_NS {

class BondZero : public Bond {
 public:
  BondZero(class LAMMPS *);
  virtual ~BondZero();
  virtual void compute(int, int);
  virtual void settings(int, char **);

  void coeff(int, char **);
  double equilibrium_distance(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_data(FILE *);

  double single(int, double, int, int, double &);

 protected:
  double *r0;
  int coeffflag;

  virtual void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

UNDOCUMENTED

E: Incorrect args for bond coefficients

Self-explanatory.  Check the input script or data file.

*/
