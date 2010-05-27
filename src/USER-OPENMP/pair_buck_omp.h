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

PairStyle(buck/omp,PairBuckOMP)

#else

#ifndef LMP_PAIR_BUCK_OMP_H
#define LMP_PAIR_BUCK_OMP_H

#include "pair_omp.h"

namespace LAMMPS_NS {

class PairBuckOMP : public PairOMP {
 public:

 PairBuckOMP(class LAMMPS *);
  ~PairBuckOMP();

  virtual void compute(int, int);

  virtual double single(int, int, int, int,
			double, double, double, double &);

  virtual void settings(int, char **);
  virtual void coeff(int, char **);

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

 private:
  double cut_global;
  double **cut;
  double **a,**rho,**c;
  double **rhoinv,**buck1,**buck2,**offset;

  void allocate();
};

}

#endif
#endif
