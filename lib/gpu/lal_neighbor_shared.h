/***************************************************************************
                              neighbor_shared.h
                             -------------------
                            W. Michael Brown (ORNL)

  Class for management of data shared by all neighbor lists

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : brownw@ornl.gov
 ***************************************************************************/

#ifndef LAL_NEIGHBOR_SHARED_H
#define LAL_NEIGHBOR_SHARED_H

#if defined(USE_OPENCL)
#include "geryon/ocl_kernel.h"
#include "geryon/ocl_texture.h"
using namespace ucl_opencl;
#elif defined(USE_CUDART)
#include "geryon/nvc_kernel.h"
#include "geryon/nvc_texture.h"
using namespace ucl_cudart;
#elif defined(USE_HIP)
#include "geryon/hip_kernel.h"
#include "geryon/hip_texture.h"
using namespace ucl_hip;
#else
#include "geryon/nvd_kernel.h"
#include "geryon/nvd_texture.h"
using namespace ucl_cudadr;
#endif

namespace LAMMPS_AL {

class NeighborShared {
 public:
  NeighborShared() : _compiled(false) {}
  ~NeighborShared() { clear(); }

  /// Free all memory on host and device
  void clear();

  /// Texture for cached position/type access with CUDA
  UCL_Texture neigh_tex;

  /// Compile kernels for neighbor lists
  void compile_kernels(UCL_Device &dev, const int gpu_nbor,
                       const std::string flags);

  // ----------------------------- Kernels
  UCL_Program *nbor_program, *build_program;
  UCL_Kernel k_nbor, k_cell_id, k_cell_counts, k_build_nbor;
  UCL_Kernel k_transpose, k_special;

 private:
  bool _compiled;
  int _gpu_nbor;
};

}

#endif
