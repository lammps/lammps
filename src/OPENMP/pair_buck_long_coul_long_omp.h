/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(buck/long/coul/long/omp,PairBuckLongCoulLongOMP);
// clang-format on
#else

#ifndef LMP_PAIR_BUCK_LONG_COUL_LONG_OMP_H
#define LMP_PAIR_BUCK_LONG_COUL_LONG_OMP_H

#include "pair_buck_long_coul_long.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class PairBuckLongCoulLongOMP : public PairBuckLongCoulLong, public ThrOMP {

 public:
  PairBuckLongCoulLongOMP(class LAMMPS *);

  void compute(int, int) override;
  void compute_inner() override;
  void compute_middle() override;
  void compute_outer(int, int) override;

 private:
  template <const int EVFLAG, const int EFLAG, const int NEWTON_PAIR, const int CTABLE,
            const int LJTABLE, const int ORDER1, const int ORDER6>
  void eval(int, int, ThrData *const);

  template <const int EVFLAG, const int EFLAG, const int NEWTON_PAIR, const int CTABLE,
            const int LJTABLE, const int ORDER1, const int ORDER6>
  void eval_outer(int, int, ThrData *const);

  void eval_inner(int, int, ThrData *const);
  void eval_middle(int, int, ThrData *const);
};

}    // namespace LAMMPS_NS

#endif
#endif
