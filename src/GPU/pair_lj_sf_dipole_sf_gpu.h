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

PairStyle(lj/sf/dipole/sf/gpu,PairLJSFDipoleSFGPU)

#else

#ifndef LMP_PAIR_LJ_SF_DIPOLE_SF_GPU_H
#define LMP_PAIR_LJ_SF_DIPOLE_SF_GPU_H

#include "pair_lj_sf_dipole_sf.h"

namespace LAMMPS_NS {

class PairLJSFDipoleSFGPU : public PairLJSFDipoleSF {
 public:
  PairLJSFDipoleSFGPU(LAMMPS *lmp);
  ~PairLJSFDipoleSFGPU();
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

E: Pair style dipole/cut/gpu requires atom attribute q

The atom style defined does not have this attribute.

E: Cannot use newton pair with dipole/cut/gpu pair style

Self-explanatory.

E: Cannot (yet) use 'electron' units with dipoles

This feature is not yet supported.

*/
