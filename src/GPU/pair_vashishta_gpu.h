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
PairStyle(vashishta/gpu,PairVashishtaGPU);
// clang-format on
#else

#ifndef LMP_PAIR_VASHISHTA_GPU_H
#define LMP_PAIR_VASHISHTA_GPU_H

#include "pair_vashishta.h"

namespace LAMMPS_NS {

class PairVashishtaGPU : public PairVashishta {
 public:
  PairVashishtaGPU(class LAMMPS *);
  ~PairVashishtaGPU();
  void compute(int, int);
  double init_one(int, int);
  void init_style();

  enum { GPU_FORCE, GPU_NEIGH, GPU_HYB_NEIGH };

 protected:
  void allocate();
  int gpu_allocated;
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

E: Pair style vashishta/gpu requires atom IDs

This is a requirement to use this potential.

E: Pair style vashishta/gpu requires newton pair off

See the newton command.  This is a restriction to use this potential.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

*/
