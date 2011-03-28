/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(lj/cut/tgpu,PairLJCutTGPU)

#else

#ifndef LMP_PAIR_LJ_LIGHT_TGPU_H
#define LMP_PAIR_LJ_LIGHT_TGPU_H

#include "pair_lj_cut.h"
#include "pair_omp_gpu.h"

namespace LAMMPS_NS {

class PairLJCutTGPU : public PairLJCut {
 public:
  PairLJCutTGPU(LAMMPS *lmp);
  ~PairLJCutTGPU();
  void cpu_compute(int, int, int, int, int *, int *, int **);
  void compute(int, int);
  void init_style();
  double memory_usage();

 enum { GPU_PAIR, GPU_NEIGH };

 private:
  int gpu_mode;
  double cpu_time;
  int *gpulist;
  
  PairOMPGPU *omp;
};

}
#endif
#endif

