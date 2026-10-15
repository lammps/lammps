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
   Contributing author: Claude Fable 5.1 (Anthropic), under the direction of
   Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(sw/angle/table/omp,PairSWAngleTableOMP);
// clang-format on
#else

#ifndef LMP_PAIR_SW_ANGLE_TABLE_OMP_H
#define LMP_PAIR_SW_ANGLE_TABLE_OMP_H

#include "pair_sw_angle_table.h"
#include "thr_omp.h"

namespace LAMMPS_NS {

class PairSWAngleTableOMP : public PairSWAngleTable, public ThrOMP {

 public:
  PairSWAngleTableOMP(class LAMMPS *);

  void compute(int, int) override;
  double memory_usage() override;

 private:
  template <int EVFLAG, int EFLAG> void eval(int ifrom, int ito, ThrData *const thr);
};

}    // namespace LAMMPS_NS

#endif
#endif
