/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef BOND_HARMONIC_H
#define BOND_HARMONIC_H

#include "stdio.h"
#include "bond.h"

namespace LAMMPS_NS {

class BondHarmonic : public Bond {
 public:
  BondHarmonic(class LAMMPS *);
  ~BondHarmonic();
  void compute(int, int);
  void coeff(int, char **);
  double equilibrium_distance(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  double single(int, double, int, int);

 private:
  double *k,*r0;

  void allocate();
};

}

#endif
