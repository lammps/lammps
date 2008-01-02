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

#ifndef PAIR_LJ_SMOOTH_H
#define PAIR_LJ_SMOOTH_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLJSmooth : public Pair {
 public:
  PairLJSmooth(class LAMMPS *);
  ~PairLJSmooth();
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
  double cut_inner_global,cut_global;
  double **cut,**cut_inner,**cut_inner_sq;
  double **epsilon,**sigma;
  double **lj1,**lj2,**lj3,**lj4;
  double **ljsw0,**ljsw1,**ljsw2,**ljsw3,**ljsw4;
  double **offset;

  void allocate();
};

}

#endif
