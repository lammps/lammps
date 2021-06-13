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
PairStyle(tip4p/cut/omp,PairTIP4PCutOMP);
// clang-format on
#else

#ifndef LMP_PAIR_TIP4P_CUT_OMP_H
#define LMP_PAIR_TIP4P_CUT_OMP_H

#include "pair_tip4p_cut.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class PairTIP4PCutOMP : public PairTIP4PCut, public ThrOMP {

 public:
  PairTIP4PCutOMP(class LAMMPS *);
  virtual ~PairTIP4PCutOMP();

  virtual void compute(int, int);
  virtual double memory_usage();

 private:
  dbl3_t *newsite_thr;
  int3_t *hneigh_thr;

  template <int, int, int> void eval(int, int, ThrData *const);
  void compute_newsite_thr(const dbl3_t &, const dbl3_t &, const dbl3_t &, dbl3_t &) const;
};

}    // namespace LAMMPS_NS

#endif
#endif
