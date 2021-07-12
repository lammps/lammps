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

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(gran/hertz/history/omp,PairGranHertzHistoryOMP);
// clang-format on
#else

#ifndef LMP_PAIR_GRAN_HERTZ_HISTORY_OMP_H
#define LMP_PAIR_GRAN_HERTZ_HISTORY_OMP_H

#include "pair_gran_hertz_history.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class PairGranHertzHistoryOMP : public PairGranHertzHistory, public ThrOMP {

 public:
  PairGranHertzHistoryOMP(class LAMMPS *);

  virtual void compute(int, int);
  virtual double memory_usage();

 private:
  template <int EVFLAG, int SHEARUPDATE, int NEWTON_PAIR>
  void eval(int ifrom, int ito, ThrData *const thr);
};

}    // namespace LAMMPS_NS

#endif
#endif
