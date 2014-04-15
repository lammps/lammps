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

PairStyle(resquared/gpu,PairRESquaredGPU)

#else

#ifndef LMP_PAIR_RESQUARED_GPU_H
#define LMP_PAIR_RESQUARED_GPU_H

#include "pair_resquared.h"

namespace LAMMPS_NS {

class PairRESquaredGPU : public PairRESquared {
 public:
  PairRESquaredGPU(LAMMPS *lmp);
  ~PairRESquaredGPU();
  void cpu_compute(int, int, int, int, int *, int *, int **);
  void compute(int, int);
  void init_style();
  double memory_usage();

  enum { GPU_FORCE, GPU_NEIGH, GPU_HYB_NEIGH };

 private:
  int gpu_mode;
  double cpu_time;
  int quat_nmax;
  double **quat;
};

}
#endif
#endif

/* ERROR/WARNING messages:

E: Pair resquared/gpu requires atom style ellipsoid

Self-explanatory.

E: Insufficient memory on accelerator

There is insufficient memory on one of the devices specified for the gpu
package

E: Cannot use newton pair with resquared/gpu pair style

Self-explanatory.

E: Pair resquared/gpu requires atoms with same type have same shape

Self-explanatory.

*/
