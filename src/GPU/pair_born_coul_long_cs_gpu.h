/* -*- c++ -*- ----------------------------------------------------------
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

PairStyle(born/coul/long/cs/gpu,PairBornCoulLongCSGPU)

#else

#ifndef LMP_PAIR_BORN_COUL_LONG_CS_GPU_H
#define LMP_PAIR_BORN_COUL_LONG_CS_GPU_H

#include "pair_born_coul_long_cs.h"

namespace LAMMPS_NS {

class PairBornCoulLongCSGPU : public PairBornCoulLongCS {
 public:
  PairBornCoulLongCSGPU(LAMMPS *lmp);
  ~PairBornCoulLongCSGPU();
  void cpu_compute(int, int, int, int, int *, int *, int **);
  void compute(int, int);
  void init_style();
  double memory_usage();

 enum { GPU_FORCE, GPU_NEIGH, GPU_HYB_NEIGH };

 private:
  int gpu_mode;
  double cpu_time;
};

}
#endif
#endif

/* ERROR/WARNING messages:

E: Insufficient memory on accelerator

There is insufficient memory on one of the devices specified for the gpu
package

E: Pair style born/coul/long/cs/gpu requires atom attribute q

The atom style defined does not have this attribute.

E: Cannot use newton pair with born/coul/long/cs/gpu pair style

Self-explanatory.

E: Pair style requires a KSpace style

No kspace style is defined.

*/
