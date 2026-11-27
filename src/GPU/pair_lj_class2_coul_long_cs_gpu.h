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
PairStyle(lj/class2/coul/long/cs/gpu,PairLJClass2CoulLongCSGPU);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_CLASS2_COUL_LONG_CS_GPU_H
#define LMP_PAIR_LJ_CLASS2_COUL_LONG_CS_GPU_H

#include "pair_lj_class2_coul_long_cs.h"

namespace LAMMPS_NS {

class PairLJClass2CoulLongCSGPU : public PairLJClass2CoulLongCS {
 public:
  PairLJClass2CoulLongCSGPU(LAMMPS *lmp);
  ~PairLJClass2CoulLongCSGPU() override;
  void compute(int, int) override;
  void init_style() override;
  double memory_usage() override;

  enum { GPU_FORCE, GPU_NEIGH, GPU_HYB_NEIGH };

 private:
  int gpu_mode;
};

}    // namespace LAMMPS_NS
#endif
#endif
