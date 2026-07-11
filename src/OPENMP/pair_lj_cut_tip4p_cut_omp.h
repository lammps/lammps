/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

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
PairStyle(lj/cut/tip4p/cut/omp,PairLJCutTIP4PCutOMP);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_CUT_TIP4P_CUT_OMP_H
#define LMP_PAIR_LJ_CUT_TIP4P_CUT_OMP_H

#include "pair_lj_cut_tip4p_cut.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class PairLJCutTIP4PCutOMP : public PairLJCutTIP4PCut, public ThrOMP {

 public:
  PairLJCutTIP4PCutOMP(class LAMMPS *);
  ~PairLJCutTIP4PCutOMP() override;

  void compute(int, int) override;
  double memory_usage() override;

 private:
  dbl3_t *newsite_thr;
  int3_t *hneigh_thr;

  template <int, int, int> void eval(int, int, ThrData *const);
  void compute_newsite_thr(const dbl3_t &, const dbl3_t &, const dbl3_t &, dbl3_t &) const;
};

}    // namespace LAMMPS_NS

#endif
#endif
