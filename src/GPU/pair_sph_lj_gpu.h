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
PairStyle(sph/lj/gpu,PairSPHLJGPU);
// clang-format on
#else

#ifndef LMP_PAIR_SPH_LJ_GPU_H
#define LMP_PAIR_SPH_LJ_GPU_H

#include "pair_sph_lj.h"

namespace LAMMPS_NS {

class PairSPHLJGPU : public PairSPHLJ {
 public:
  PairSPHLJGPU(LAMMPS *lmp);
  ~PairSPHLJGPU() override;
  void cpu_compute(int, int, int, int, int *, int *, int **);
  void compute(int, int) override;
  void init_style() override;
  double memory_usage() override;

  enum { GPU_FORCE, GPU_NEIGH, GPU_HYB_NEIGH };

  void *drhoE_pinned;
  bool acc_float;

 private:
  int gpu_mode;
  double cpu_time;
};

}    // namespace LAMMPS_NS
#endif
#endif
