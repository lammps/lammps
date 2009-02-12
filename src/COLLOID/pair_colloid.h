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

#ifndef PAIR_COLLOID_H
#define PAIR_COLLOID_H

#include "pair.h"

namespace LAMMPS_NS {

class PairColloid : public Pair {
 public:
  PairColloid(class LAMMPS *);
  ~PairColloid();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  double single(int, int, int, int, double, double, double, double &);

 private:
  double cut_global;
  double **cut;
  double **a12,**d1,**d2,**diameter,**a1,**a2,**offset;
  double **sigma,**sigma3,**sigma6;
  double **lj1,**lj2,**lj3,**lj4;
  int **form;

  void allocate();
};

}

#endif
