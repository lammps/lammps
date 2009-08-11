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

#ifndef PAIR_GPU_H
#define PAIR_GPU_H

#include "pair_gayberne.h"
#define MAX_GPU_THREADS 4

namespace LAMMPS_NS {

class PairGayBerneGPU : public PairGayBerne {
 public:
  PairGayBerneGPU(LAMMPS *lmp);
  ~PairGayBerneGPU();
  void compute(int, int);
  void settings(int, char **);
  void init_style();
  double memory_usage();
 
  enum { ONE_NODE, ONE_GPU, MULTI_GPU };

 private:  
  int ij_size;
  int *ij[MAX_GPU_THREADS], *ij_new[MAX_GPU_THREADS], *olist[MAX_GPU_THREADS];
 
  int my_thread, nthreads, thread_inum[MAX_GPU_THREADS], omp_chunk;
 
  int last_neighbor, multi_gpu_mode, multi_gpu_param;
};

}
#endif
