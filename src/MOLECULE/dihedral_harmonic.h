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

#ifndef DIHEDRAL_HARMONIC_H
#define DIHEDRAL_HARMONIC_H

#include "stdio.h"
#include "dihedral.h"

class DihedralHarmonic : public Dihedral {
 public:
  DihedralHarmonic() {}
  ~DihedralHarmonic();
  void compute(int, int);
  void coeff(int, int, char **);
  void write_restart(FILE *);
  void read_restart(FILE *);

 private:
  double *k,*cos_shift,*sin_shift;
  int *sign,*multiplicity;

  void allocate();
};

#endif
