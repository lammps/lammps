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
PairStyle(lj/spica/coul/msm,PairLJSPICACoulMSM);
PairStyle(lj/sdk/coul/msm,PairLJSPICACoulMSM);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_SPICA_COUL_MSM_H
#define LMP_PAIR_LJ_SPICA_COUL_MSM_H

#include "pair_lj_spica_coul_long.h"

namespace LAMMPS_NS {

class PairLJSPICACoulMSM : public PairLJSPICACoulLong {
 public:
  PairLJSPICACoulMSM(class LAMMPS *);
  void compute(int, int) override;
  double single(int, int, int, int, double, double, double, double &) override;
  void *extract(const char *, int &) override;

 private:
  template <int EVFLAG, int EFLAG, int NEWTON_PAIR> void eval_msm();
};

}    // namespace LAMMPS_NS

#endif
#endif
