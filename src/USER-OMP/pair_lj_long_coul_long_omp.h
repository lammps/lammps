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

PairStyle(lj/long/coul/long/omp,PairLJLongCoulLongOMP)

#else

#ifndef LMP_PAIR_LJ_LONG_COUL_LONG_OMP_H
#define LMP_PAIR_LJ_LONG_COUL_LONG_OMP_H

#include "pair_lj_long_coul_long.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class PairLJLongCoulLongOMP : public PairLJLongCoulLong, public ThrOMP {

 public:
  PairLJLongCoulLongOMP(class LAMMPS *);

  virtual void compute(int, int);
  virtual void compute_inner();
  virtual void compute_middle();
  virtual void compute_outer(int, int);
  virtual double memory_usage();

 private:
  template <const int EVFLAG, const int EFLAG,
    const int NEWTON_PAIR, const int CTABLE, const int LJTABLE,
    const int ORDER1, const int ORDER6 >
  void eval(int, int, ThrData * const);

  template <const int EVFLAG, const int EFLAG,
    const int NEWTON_PAIR, const int CTABLE, const int LJTABLE,
    const int ORDER1, const int ORDER6 >
  void eval_outer(int, int, ThrData * const);


  void eval_inner(int, int, ThrData *const);
  void eval_middle(int, int, ThrData *const);

  

};

}

#endif
#endif
