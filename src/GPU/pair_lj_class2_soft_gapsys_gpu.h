/* ----------------------------------------------------------------------
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
PairStyle(lj/class2/soft/gapsys/gpu,PairLJClass2SoftGapsysGPU);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_CLASS2_SOFT_GAPSYS_GPU_H
#define LMP_PAIR_LJ_CLASS2_SOFT_GAPSYS_GPU_H

#include "pair_lj_class2_soft_gapsys.h"

namespace LAMMPS_NS {

class PairLJClass2SoftGapsysGPU : public PairLJClass2SoftGapsys {
 public:
  PairLJClass2SoftGapsysGPU(LAMMPS *lmp);
  ~PairLJClass2SoftGapsysGPU() override;
  void compute(int, int) override;
  void init_style() override;
  double init_one(int, int) override;
  double memory_usage() override;

  enum { GPU_FORCE, GPU_NEIGH, GPU_HYB_NEIGH };

 private:
  int gpu_mode;
  bool update_coeffs;    // per-type parameters changed since the last copy to the GPU
};

}    // namespace LAMMPS_NS
#endif
#endif
