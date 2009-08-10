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

#ifndef PAIR_LJ96_CUT_H
#define PAIR_LJ96_CUT_H

#include "pair.h"

namespace LAMMPS_NS {

class PairLJ96Cut : public Pair {
 public:
  PairLJ96Cut(class LAMMPS *);
  ~PairLJ96Cut();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  void init_list(int, class NeighList *);
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  double single(int, int, int, int, double, double, double, double &);

  void compute_inner();
  void compute_middle();
  void compute_outer(int, int);

 protected:
  double cut_global;
  double **cut;
  double **epsilon,**sigma;
  double **lj1,**lj2,**lj3,**lj4,**offset;
  double *cut_respa;

  void allocate();
};

}

#endif
