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

#ifndef DIHEDRAL_HARMONIC_H
#define DIHEDRAL_HARMONIC_H

#include "stdio.h"
#include "dihedral.h"

namespace LAMMPS_NS {

class DihedralHarmonic : public Dihedral {
 public:
  DihedralHarmonic(class LAMMPS *);
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

}

#endif
