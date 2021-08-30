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
PairStyle(nm/cut/coul/cut/omp,PairNMCutCoulCutOMP);
// clang-format on
#else

#ifndef LMP_PAIR_NM_CUT_COUL_CUT_OMP_H
#define LMP_PAIR_NM_CUT_COUL_CUT_OMP_H

#include "pair_nm_cut_coul_cut.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class PairNMCutCoulCutOMP : public PairNMCutCoulCut, public ThrOMP {

 public:
  PairNMCutCoulCutOMP(class LAMMPS *);

  virtual void compute(int, int);
  virtual double memory_usage();

 private:
  template <int EVFLAG, int EFLAG, int NEWTON_PAIR>
  void eval(int ifrom, int ito, ThrData *const thr);
};

}    // namespace LAMMPS_NS

#endif
#endif
