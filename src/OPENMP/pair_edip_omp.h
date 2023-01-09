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
PairStyle(edip/omp,PairEDIPOMP);
// clang-format on
#else

#ifndef LMP_PAIR_EDIP_OMP_H
#define LMP_PAIR_EDIP_OMP_H

#include "pair_edip.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class PairEDIPOMP : public PairEDIP, public ThrOMP {

 public:
  PairEDIPOMP(class LAMMPS *);

  void compute(int, int) override;
  double memory_usage() override;

 private:
  template <int EVFLAG, int EFLAG, int VFLAG_ATOM>
  void eval(int ifrom, int ito, ThrData *const thr);
};

}    // namespace LAMMPS_NS

#endif
#endif
