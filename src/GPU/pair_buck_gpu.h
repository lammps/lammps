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

PairStyle(buck/gpu,PairBuckGPU)

#else

#ifndef LMP_PAIR_BUCK_GPU_H
#define LMP_PAIR_BUCK_GPU_H

#include "pair_buck.h"

namespace LAMMPS_NS {

class PairBuckGPU : public PairBuck {
 public:
  PairBuckGPU(LAMMPS *lmp);
  ~PairBuckGPU();
  void cpu_compute(int, int, int, int, int *, int *, int **);
  void compute(int, int);
  void init_style();
  double memory_usage();

 enum { GPU_FORCE, GPU_NEIGH, GPU_HYB_NEIGH };

 private:
  int gpu_mode;
  double cpu_time;
  int *gpulist;
};

}
#endif
#endif

/* ERROR/WARNING messages:

E: Insufficient memory on accelerator

There is insufficient memory on one of the devices specified for the gpu
package

E: Cannot use newton pair with buck/gpu pair style

Self-explanatory.

*/
