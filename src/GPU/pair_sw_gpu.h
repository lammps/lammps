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
PairStyle(sw/gpu,PairSWGPU);
// clang-format on
#else

#ifndef LMP_PAIR_SW_GPU_H
#define LMP_PAIR_SW_GPU_H

#include "pair_sw.h"

namespace LAMMPS_NS {

class PairSWGPU : public PairSW {
 public:
  PairSWGPU(class LAMMPS *);
  ~PairSWGPU() override;
  void compute(int, int) override;
  double init_one(int, int) override;
  void init_style() override;

  enum { GPU_FORCE, GPU_NEIGH, GPU_HYB_NEIGH };

 protected:
  void allocate() override;

  int gpu_mode;
  double cpu_time;
};

}    // namespace LAMMPS_NS

#endif
#endif
