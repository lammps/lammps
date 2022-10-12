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
PairStyle(lj/cut/tip4p/long/opt,PairLJCutTIP4PLongOpt);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_CUT_TIP4P_LONG_OPT_H
#define LMP_PAIR_LJ_CUT_TIP4P_LONG_OPT_H

#include "pair_lj_cut_tip4p_long.h"

namespace LAMMPS_NS {

class PairLJCutTIP4PLongOpt : public PairLJCutTIP4PLong {
 public:
  PairLJCutTIP4PLongOpt(class LAMMPS *);

  void compute(int, int) override;
  double memory_usage() override;

 protected:
  template <const int, const int, const int, const int> void eval();
  void compute_newsite_opt(const double *, const double *, const double *, double *) const;
};

}    // namespace LAMMPS_NS

#endif
#endif
