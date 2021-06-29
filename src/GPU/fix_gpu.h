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

#ifdef FIX_CLASS
// clang-format off
FixStyle(GPU,FixGPU);
// clang-format on
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
  void post_force_respa(int, int, int);
  double memory_usage();

  double binsize(const double subx, const double suby, const double subz, const int nlocal,
                 const double cut);

 private:
  int _gpu_mode;
  int _nlevels_respa;
  double _particle_split;
  double _binsize;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: No OpenMP support compiled in

An OpenMP flag is set, but LAMMPS was not built with
OpenMP support.

E: GPU package does not (yet) work with atom_style template

Self-explanatory.

E: Cannot use package gpu neigh yes with triclinic box

This is a current restriction in LAMMPS.

W: Using package gpu without any pair style defined

Self-explanatory.

E: Cannot use neigh_modify exclude with GPU neighbor builds

This is a current limitation of the GPU implementation
in LAMMPS.

E: Too many neighbors on GPU. Use neigh_modify one to increase limit.

The expected maximum number of neighbors is determined in the GPU package
automatically. This error means the actual number of neighbors is exceeding
the expected value. Use neigh_modify one command to increase GPU allocations
(e.g. doubling this value doubles the GPU allocation).

*/
