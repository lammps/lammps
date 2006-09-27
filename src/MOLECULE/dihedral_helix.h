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

#ifndef DIHEDRAL_HELIX_H
#define DIHEDRAL_HELIX_H

#include "stdio.h"
#include "dihedral.h"

class DihedralHelix : public Dihedral {
 public:
  DihedralHelix() {}
  ~DihedralHelix();
  void compute(int, int);
  void coeff(int, int, char **);
  void write_restart(FILE *);
  void read_restart(FILE *);

 private:
  double *aphi,*bphi,*cphi;

  void allocate();
};

#endif
