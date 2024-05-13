/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS Development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(lj/cut/coul/long/soft/gpu,PairLJCutCoulLongSoftGPU);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_CUT_COUL_LONG_SOFT_GPU_H
#define LMP_PAIR_LJ_CUT_COUL_LONG_SOFT_GPU_H

#include "pair_lj_cut_coul_long_soft.h"

namespace LAMMPS_NS {

class PairLJCutCoulLongSoftGPU : public PairLJCutCoulLongSoft {
 public:
  PairLJCutCoulLongSoftGPU(LAMMPS *lmp);
  ~PairLJCutCoulLongSoftGPU() override;
  void cpu_compute(int, int, int, int, int *, int *, int **);
  void compute(int, int) override;
  void init_style() override;
  void reinit() override;
  double memory_usage() override;

  enum { GPU_FORCE, GPU_NEIGH, GPU_HYB_NEIGH };

 private:
  int gpu_mode;
  double cpu_time;
};

}    // namespace LAMMPS_NS
#endif
#endif
