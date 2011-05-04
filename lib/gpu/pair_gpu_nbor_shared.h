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
   Contributing authors: Mike Brown (ORNL), brownw@ornl.gov
------------------------------------------------------------------------- */

#ifndef PAIR_GPU_NBOR_SHARED_H
#define PAIR_GPU_NBOR_SHARED_H

#ifdef USE_OPENCL

#include "geryon/ocl_kernel.h"
#include "geryon/ocl_texture.h"
using namespace ucl_opencl;

#else

#include "geryon/nvd_kernel.h"
#include "geryon/nvd_texture.h"
using namespace ucl_cudadr;

#endif

class PairGPUNborShared {
 public:
  PairGPUNborShared() : _compiled(false) {}
  ~PairGPUNborShared() { clear(); }
 
  /// Free all memory on host and device
  void clear();

  /// Texture for cached position/type access with CUDA
  UCL_Texture neigh_tex;

  /// Compile kernels for neighbor lists
  void compile_kernels(UCL_Device &dev, const bool gpu_nbor);

  // ----------------------------- Kernels
  UCL_Program *nbor_program, *build_program;
  UCL_Kernel k_nbor, k_cell_id, k_cell_counts, k_build_nbor;
  UCL_Kernel k_transpose, k_special;

 private:
  bool _compiled, _gpu_nbor;
};

#endif
