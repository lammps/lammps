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

PairStyle(lj/charmm/coul/charmm/omp,PairLJCharmmCoulCharmmOMP)

#else

#ifndef LMP_PAIR_LJ_CHARMM_COUL_CHARMM_OMP_H
#define LMP_PAIR_LJ_CHARMM_COUL_CHARMM_OMP_H

#include "pair_omp.h"

namespace LAMMPS_NS {

class PairLJCharmmCoulCharmmOMP : public PairOMP {
 public:
  PairLJCharmmCoulCharmmOMP(class LAMMPS *);
  ~PairLJCharmmCoulCharmmOMP();
  virtual void compute(int, int);
  virtual void settings(int, char **);
  virtual void coeff(int, char **);
  virtual void init_style();
  virtual double init_one(int, int);
  virtual void write_restart(FILE *);
  virtual void read_restart(FILE *);
  virtual void write_restart_settings(FILE *);
  virtual void read_restart_settings(FILE *);
  virtual double single(int, int, int, int, double, double, double, double &);
  virtual void *extract(char *, int &);

  virtual double memory_usage();

 protected:
  template <int EVFLAG, int EFLAG, int NEWTON_PAIR> void eval();

 protected:
  int implicit;
  double cut_lj_inner,cut_lj,cut_coul_inner,cut_coul;
  double cut_lj_innersq,cut_ljsq,cut_coul_innersq,cut_coulsq,cut_bothsq;
  double denom_lj,denom_coul;
  double **epsilon,**sigma,**eps14,**sigma14;
  double **lj1,**lj2,**lj3,**lj4;
  double **lj14_1,**lj14_2,**lj14_3,**lj14_4;

  void allocate();
};

}

#endif
#endif
