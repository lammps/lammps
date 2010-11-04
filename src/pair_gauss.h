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

#ifdef PAIR_CLASS

PairStyle(gauss,PairGauss)

#else

#ifndef PAIR_GAUSS_H
#define PAIR_GAUSS_H

#include "pair.h"

namespace LAMMPS_NS {

class PairGauss : public Pair {
 public:
  PairGauss(class LAMMPS *);
  ~PairGauss();
  void compute(int,int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);  
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  double single(int, int, int, int, double, double, double, double &);
  void *extract(char *, int &);
 
 private: 
  double cut_global;
  double **cut;
  double **a,**b;
  double **offset;
  
  void allocate();
};

}

#endif
#endif
