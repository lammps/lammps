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
PairStyle(morse/soft/gpu,PairMorseSoftGPU);
// clang-format on
#else

#ifndef LMP_PAIR_MORSE_SOFT_GPU_H
#define LMP_PAIR_MORSE_SOFT_GPU_H

#include "pair_morse_soft.h"

namespace LAMMPS_NS {

class PairMorseSoftGPU : public PairMorseSoft {
 public:
  PairMorseSoftGPU(LAMMPS *lmp);
  ~PairMorseSoftGPU() override;
  void compute(int, int) override;
  void init_style() override;
  double init_one(int, int) override;
  double memory_usage() override;

  enum { GPU_FORCE, GPU_NEIGH, GPU_HYB_NEIGH };

 private:
  int gpu_mode;
  bool update_coeffs;    // per-type parameters changed since the last copy to the GPU
  double **bs, **vs;
  void setup_soft_factors();
};

}    // namespace LAMMPS_NS
#endif
#endif
