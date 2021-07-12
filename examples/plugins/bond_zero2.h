/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/ Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_BOND_ZERO2_H
#define LMP_BOND_ZERO2_H

#include "bond.h"

namespace LAMMPS_NS {

class BondZero2 : public Bond {
 public:
  BondZero2(class LAMMPS *);
  virtual ~BondZero2();
  virtual void compute(int, int);
  virtual void settings(int, char **);

  void coeff(int, char **);
  double equilibrium_distance(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_data(FILE *);

  double single(int, double, int, int, double &);
  virtual void *extract(const char *, int &);

 protected:
  double *r0;
  int coeffflag;

  virtual void allocate();
};

}    // namespace LAMMPS_NS

#endif

/* ERROR/WARNING messages:

E: Illegal ... command

UNDOCUMENTED

E: Incorrect args for bond coefficients

Self-explanatory.  Check the input script or data file.

*/
