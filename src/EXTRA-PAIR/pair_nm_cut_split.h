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
PairStyle(nm/cut/split,PairNMCutSplit);
// clang-format on
#else

#ifndef LMP_PAIR_NM_CUT_SPLIT_H
#define LMP_PAIR_NM_CUT_SPLIT_H

#include "pair_nm_cut.h"
namespace LAMMPS_NS {

class PairNMCutSplit : public PairNMCut {
 public:
  PairNMCutSplit(class LAMMPS *);
  double single(int, int, int, int, double, double, double, double &) override;
  void compute(int, int) override;
};
}    // namespace LAMMPS_NS
#endif
#endif
