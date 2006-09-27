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

#ifndef IMPROPER_CLASS2_H
#define IMPROPER_CLASS2_H

#include "stdio.h"
#include "improper.h"

class ImproperClass2 : public Improper {
 public:
  ImproperClass2();
  ~ImproperClass2();
  void compute(int, int);
  void coeff(int, int, char **);
  void write_restart(FILE *);
  void read_restart(FILE *);

 private:
  double *k0,*chi0;
  double *aa_k1,*aa_k2,*aa_k3,*aa_theta0_1,*aa_theta0_2,*aa_theta0_3;
  int *setflag_i,*setflag_aa;
  double PI;

  void allocate();
  void angleangle(int, int);
  void cross(double *, double *, double *);
  double dot(double *, double *);
};

#endif
