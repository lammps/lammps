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

  /// Use a heuristic to approximate best bin size assuming uniform density
  /** This is only called by core LAMMPS for atom sort sizes **/
  inline double update_cell_size(const double subx, const double suby,
                                 const double subz, const int nlocal,
                                 const double cut) {
    if (_auto_cell_size==false || subz==0.0) return cut;
    else {
      _cell_size=best_cell_size(subx, suby, subz, nlocal, cut);
      _cached_cell_size=true;
      _cut_sort=cut;
      return _cell_size;
    }
  }

  /// Use a heuristic to approximate best bin size assuming uniform density
  double best_cell_size(const double subx, const double suby,
                        const double subz, const int nlocal,
                        const double cut);

  /// Current cutoff used for cell size determination
  inline double neighbor_cutoff() { return _neighbor_cutoff; }

  /// Current neighbor cell size
  inline double cell_size() { return _cell_size; }

  /// Return setting for auto cell size
  inline bool auto_cell_size() { return _auto_cell_size; }

  inline void setup_auto_cell_size(const bool autosize, const double cut,
                                   const int simd_size) {
    _auto_cell_size = autosize;
    _cached_cell_size = false;
    _neighbor_cutoff = cut;
    _cell_size = cut;
    _simd_size = simd_size;
    if (_simd_size < 2) _auto_cell_size = false;
  }

  /// Compile kernels for neighbor lists
  void compile_kernels(UCL_Device &dev, const int gpu_nbor,
                       const std::string &flags);

  // ----------------------------- Kernels
  UCL_Program *nbor_program, *build_program;
  UCL_Kernel k_nbor, k_cell_id, k_cell_counts, k_build_nbor;
  UCL_Kernel k_transpose, k_special;

 private:
  bool _compiled;
  int _gpu_nbor;
  bool _auto_cell_size, _cached_cell_size;
  double _neighbor_cutoff, _cell_size, _simd_size, _cut_sort;
};

}

#endif
