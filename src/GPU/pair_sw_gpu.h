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

PairStyle(sw/gpu,PairSWGPU)

#else

#ifndef LMP_PAIR_SW_GPU_H
#define LMP_PAIR_SW_GPU_H

#include "pair_sw.h"

namespace LAMMPS_NS {

class PairSWGPU : public PairSW {
 public:
  PairSWGPU(class LAMMPS *);
  ~PairSWGPU();
  void compute(int, int);
  double init_one(int, int);
  void init_style();

 enum { GPU_FORCE, GPU_NEIGH, GPU_HYB_NEIGH };

 protected:
  void allocate();

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

E: Pair style sw/gpu requires newton pair on

See the newton command.  This is a restriction to use the SW
potential.

E: Pair style sw/gpu is currently limited to one element.

Self-explanatory.

*/

