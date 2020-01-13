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

/* ----------------------------------------------------------------------
   Contributing author: Vsevolod Nikolskiy (HSE)
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(lj/cut/tip4p/long/gpu,PairLJCutTIP4PLongGPU)

#else

#ifndef LMP_PAIR_LJ_TIP4P_LONG_GPU_H
#define LMP_PAIR_LJ_TIP4P_LONG_GPU_H

#include "pair_lj_cut_tip4p_long.h"

namespace LAMMPS_NS {

class PairLJCutTIP4PLongGPU : public PairLJCutTIP4PLong {
 public:
  PairLJCutTIP4PLongGPU(LAMMPS *lmp);
  ~PairLJCutTIP4PLongGPU();
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
