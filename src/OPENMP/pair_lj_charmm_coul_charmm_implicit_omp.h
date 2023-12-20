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
PairStyle(lj/charmm/coul/charmm/implicit/omp,PairLJCharmmCoulCharmmImplicitOMP);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_CHARMM_COUL_CHARMM_IMPLICIT_OMP_H
#define LMP_PAIR_LJ_CHARMM_COUL_CHARMM_IMPLICIT_OMP_H

#include "pair_lj_charmm_coul_charmm_implicit.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class PairLJCharmmCoulCharmmImplicitOMP : public PairLJCharmmCoulCharmmImplicit, public ThrOMP {

 public:
  PairLJCharmmCoulCharmmImplicitOMP(class LAMMPS *);

  void compute(int, int) override;
  double memory_usage() override;

 private:
  template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
  void eval(int ifrom, int ito, ThrData *const thr);
};

}    // namespace LAMMPS_NS

#endif
#endif
