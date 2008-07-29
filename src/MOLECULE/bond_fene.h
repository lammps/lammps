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

#ifndef BOND_FENE_H
#define BOND_FENE_H

#include "stdio.h"
#include "bond.h"

namespace LAMMPS_NS {

class BondFENE : public Bond {
 public:
  BondFENE(class LAMMPS *);
  ~BondFENE();
  void compute(int, int);
  void coeff(int, char **);
  void init_style();
  double equilibrium_distance(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  double single(int, double, int, int);

 private:
  double TWO_1_3;
  double *k,*r0,*epsilon,*sigma;

  void allocate();
};

}

#endif
