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

PairStyle(eam/lj/gpu,PairEAMLJGPU)

#else

#ifndef LMP_PAIR_EAM_LJ_GPU_H
#define LMP_PAIR_EAM_LJ_GPU_H

#include "stdio.h"
#include "pair_eam.h"

namespace LAMMPS_NS {

class PairEAMLJGPU : public PairEAM {
 public:

  PairEAMLJGPU(class LAMMPS *);
  virtual ~PairEAMLJGPU();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  double memory_usage();
  void cpu_compute(int, int, int, int, int *, int *, int **);
  void cpu_compute_energy(int, int, int, int, int *, int *, int **);
  
  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);

 enum { GPU_FORCE, GPU_NEIGH, GPU_HYB_NEIGH };

 protected:
  double cut_global;
  double **cut;
  double **epsilon,**sigma;
  double **lj1,**lj2,**lj3,**lj4,**offset;

  void allocate();

 private:
  int gpu_mode;
  double cpu_time;
  int *gpulist;
  void *fp_pinned;
  bool fp_single;

  
};

}

#endif
#endif
