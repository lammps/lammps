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
PairStyle(coul/cut/global/omp,PairCoulCutGlobalOMP);
// clang-format on
#else

#ifndef LMP_PAIR_COUL_CUT_GLOBAL_OMP_H
#define LMP_PAIR_COUL_CUT_GLOBAL_OMP_H

#include "pair_coul_cut_omp.h"

namespace LAMMPS_NS {

class PairCoulCutGlobalOMP : public PairCoulCutOMP {
 public:
  PairCoulCutGlobalOMP(class LAMMPS *lmp) : PairCoulCutOMP(lmp) {}
  void coeff(int, char **) override;
  void *extract(const char *, int &) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
