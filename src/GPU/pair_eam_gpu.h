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

PairStyle(eam/gpu,PairEAMGPU)

#else

#ifndef LMP_PAIR_EAM_GPU_H
#define LMP_PAIR_EAM_GPU_H

#include "stdio.h"
#include "pair_eam.h"

namespace LAMMPS_NS {

class PairEAMGPU : public PairEAM {
 public:

  PairEAMGPU(class LAMMPS *);
  virtual ~PairEAMGPU();
  void compute(int, int);
  void init_style();
  double memory_usage();

  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);

 enum { GPU_FORCE, GPU_NEIGH, GPU_HYB_NEIGH };

 protected:
  int gpu_mode;
  double cpu_time;
  int *gpulist;
  void *fp_pinned;
  bool fp_single;  
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Insufficient memory on accelerator

There is insufficient memory on one of the devices specified for the gpu
package

E: Cannot use newton pair with eam/gpu pair style

Self-explanatory.

*/
