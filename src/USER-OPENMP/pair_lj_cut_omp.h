/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(lj/cut/omp,PairLJCutOMP)

#else

#ifndef LMP_PAIR_LJ_CUT_OMP_H
#define LMP_PAIR_LJ_CUT_OMP_H

#include "pair_omp.h"

namespace LAMMPS_NS {

class PairLJCutOMP : public PairOMP {
 public:

  PairLJCutOMP(class LAMMPS *);
  ~PairLJCutOMP();

  virtual void compute(int, int);
  virtual void compute_inner();
  virtual void compute_middle();
  virtual void compute_outer(int, int);

  virtual double single(int, int, int, int, 
			double, double, double, double &);

  virtual void settings(int, char **);
  virtual void coeff(int, char **);

  virtual void init_style();
  virtual void init_list(int, class NeighList *);
  virtual double init_one(int, int);

  virtual void write_restart(FILE *);
  virtual void read_restart(FILE *);
  virtual void write_restart_settings(FILE *);
  virtual void read_restart_settings(FILE *);

  virtual double memory_usage();

 protected:
  template <int EVFLAG, int EFLAG, int NEWTON_PAIR> void eval();
  template <int NEWTON_PAIR> void eval_inner();
  template <int NEWTON_PAIR> void eval_middle();
  template <int EVFLAG, int EFLAG, int VFLAG, int NEWTON_PAIR> 
  void eval_outer();

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
#endif
