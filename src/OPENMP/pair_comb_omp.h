/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(comb/omp,PairCombOMP);
// clang-format on
#else

#ifndef LMP_PAIR_COMB_OMP_H
#define LMP_PAIR_COMB_OMP_H

#include "pair_comb.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class PairCombOMP : public PairComb, public ThrOMP {

 public:
  PairCombOMP(class LAMMPS *);

  void compute(int, int) override;
  double memory_usage() override;

  double yasu_char(double *, int &) override;

 private:
  template <int EVFLAG, int EFLAG, int VFLAG_ATOM>
  void eval(int ifrom, int ito, ThrData *const thr);

  void Short_neigh_thr();
};

}    // namespace LAMMPS_NS

#endif
#endif
