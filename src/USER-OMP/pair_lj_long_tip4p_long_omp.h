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

PairStyle(lj/long/tip4p/long/omp,PairLJLongTIP4PLongOMP)

#else

#ifndef LMP_PAIR_LJ_LONG_TIP4P_LONG_OMP_H
#define LMP_PAIR_LJ_LONG_TIP4P_LONG_OMP_H

#include "pair_lj_long_tip4p_long.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class PairLJLongTIP4PLongOMP : public PairLJLongTIP4PLong, public ThrOMP {

 public:
  PairLJLongTIP4PLongOMP(class LAMMPS *);
  virtual ~PairLJLongTIP4PLongOMP();

  virtual void compute(int, int);
  virtual void compute_inner();
  virtual void compute_middle();
  virtual void compute_outer(int, int);
  virtual double memory_usage();

 private:
  dbl3_t *newsite_thr;
  int3_t *hneigh_thr;

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

  void compute_newsite_thr(const dbl3_t &, const dbl3_t &,
                           const dbl3_t &, dbl3_t &) const;

  

};

}

#endif
#endif
