/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   www.cs.sandia.gov/~sjplimp/lammps.html
   Steve Plimpton, sjplimp@sandia.gov, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef BOND_NONLINEAR_H
#define BOND_NONLINEAR_H

#include "stdio.h"
#include "bond.h"

class BondNonlinear : public Bond {
 public:
  BondNonlinear() {}
  ~BondNonlinear();
  void compute(int, int);
  void coeff(int, char **);
  double equilibrium_distance(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void single(int, double, int, int, double, int, double &, double &);

 private:
  double *epsilon,*r0,*lamda;

  void allocate();
};

#endif
