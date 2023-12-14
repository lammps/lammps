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
PairStyle(tersoff/zbl/gpu,PairTersoffZBLGPU);
// clang-format on
#else

#ifndef LMP_PAIR_TERSOFF_ZBL_GPU_H
#define LMP_PAIR_TERSOFF_ZBL_GPU_H

#include "pair_tersoff_zbl.h"

namespace LAMMPS_NS {

class PairTersoffZBLGPU : public PairTersoffZBL {
 public:
  PairTersoffZBLGPU(class LAMMPS *);
  ~PairTersoffZBLGPU() override;
  void compute(int, int) override;
  double init_one(int, int) override;
  void init_style() override;

  enum { GPU_FORCE, GPU_NEIGH, GPU_HYB_NEIGH };

 protected:
  void allocate() override;

  int gpu_mode;
  double cpu_time;
  int *gpulist;
};

}    // namespace LAMMPS_NS

#endif
#endif
