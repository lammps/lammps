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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(lj/cut/coul/long/cs,PairLJCutCoulLongCS);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_CUT_COUL_LONG_CS_H
#define LMP_PAIR_LJ_CUT_COUL_LONG_CS_H

#include "pair_lj_cut_coul_long.h"

namespace LAMMPS_NS {

class PairLJCutCoulLongCS : public PairLJCutCoulLong {

 public:
  PairLJCutCoulLongCS(class LAMMPS *);
  void compute(int, int) override;
  void compute_inner() override;
  void compute_middle() override;
  void compute_outer(int, int) override;
};

}    // namespace LAMMPS_NS

#endif
#endif
