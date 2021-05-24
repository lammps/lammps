/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(born/coul/wolf/cs/gpu,PairBornCoulWolfCSGPU);
// clang-format on
#else

#ifndef LMP_PAIR_BORN_COUL_WOLF_CS_GPU_H
#define LMP_PAIR_BORN_COUL_WOLF_CS_GPU_H

#include "pair_born_coul_wolf_cs.h"

namespace LAMMPS_NS {

class PairBornCoulWolfCSGPU : public PairBornCoulWolfCS {
 public:
  PairBornCoulWolfCSGPU(LAMMPS *lmp);
  ~PairBornCoulWolfCSGPU();
  void cpu_compute(int, int, int, int, int *, int *, int **);
  void compute(int, int);
  void init_style();
  double memory_usage();

  enum { GPU_FORCE, GPU_NEIGH, GPU_HYB_NEIGH };

 private:
  int gpu_mode;
  double cpu_time;
};

}    // namespace LAMMPS_NS
#endif
#endif

/* ERROR/WARNING messages:

E: Insufficient memory on accelerator

There is insufficient memory on one of the devices specified for the gpu
package

E: Cannot use newton pair with born/coul/wolf/cs/gpu pair style

Self-explanatory.

*/
