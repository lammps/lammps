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

#ifdef FIX_CLASS
// clang-format off
FixStyle(NEIGH_HISTORY/omp,FixNeighHistoryOMP);
// clang-format on
#else

#ifndef LMP_FIX_NEIGH_HISTORY_OMP_H
#define LMP_FIX_NEIGH_HISTORY_OMP_H

#include "fix_neigh_history.h"

namespace LAMMPS_NS {

class FixNeighHistoryOMP : public FixNeighHistory {

 public:
  FixNeighHistoryOMP(class LAMMPS *lmp, int narg, char **argv);
  void pre_exchange_onesided() override;
  void pre_exchange_newton() override;
  void pre_exchange_no_newton() override;
  void post_neighbor() override;
};

}    // namespace LAMMPS_NS

#endif
#endif
