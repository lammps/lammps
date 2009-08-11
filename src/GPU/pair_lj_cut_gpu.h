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

#ifndef PAIR_LJ_CUT_GPU_H
#define PAIR_LJ_CUT_GPU_H

#include "pair_lj_cut.h"

namespace LAMMPS_NS {

class PairLJCutGPU : public PairLJCut {
 public:
  PairLJCutGPU(LAMMPS *lmp);
  ~PairLJCutGPU();
  void compute(int, int);
  void settings(int, char **);
  void init_style();
  double memory_usage();

  enum { ONE_NODE, ONE_GPU, MULTI_GPU };

 private:  
  int ij_size;
  int *ij, *ij_new;
 
  int last_neighbor, multi_gpu_mode, multi_gpu_param;
};

}
#endif
