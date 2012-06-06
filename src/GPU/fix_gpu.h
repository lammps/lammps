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

#ifdef FIX_CLASS

FixStyle(GPU,FixGPU)

#else

#ifndef LMP_FIX_GPU_H
#define LMP_FIX_GPU_H

#include "fix.h"

namespace LAMMPS_NS {

class FixGPU : public Fix {
 public:
  FixGPU(class LAMMPS *, int, char **);
  ~FixGPU();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void min_post_force(int);
  double memory_usage();

 private:
  int _gpu_mode;
  double _particle_split;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Cannot use fix GPU with USER-CUDA mode enabled

You cannot use both the GPU and USER-CUDA packages
together.  Use one or the other.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot use force/neigh with triclinic box

This is a current limitation of the GPU implementation
in LAMMPS.

E: Cannot use force/hybrid_neigh with triclinic box

Self-explanatory.

E: No OpenMP support compiled in

An OpenMP flag is set, but LAMMPS was not built with
OpenMP support.

E: Cannot use pair hybrid with GPU neighbor builds

See documentation for fix gpu.

E: Fix GPU split must be positive for hybrid pair styles

Self-explanatory.

E: Cannot use neigh_modify exclude with GPU neighbor builds

This is a current limitation of the GPU implementation
in LAMMPS.

*/
