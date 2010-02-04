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
   Contributing authors: Mike Brown (SNL), wmbrown@sandia.gov
                         Peng Wang (Nvidia), penwang@nvidia.com
                         Paul Crozier (SNL), pscrozi@sandia.gov
------------------------------------------------------------------------- */

#ifndef LJ_GPU_MEMORY_H
#define LJ_GPU_MEMORY_H

#include "nvc_device.h"
#include "nvc_traits.h"
#include "pair_gpu_atom.h"
#include "pair_gpu_nbor.h"

#define BLOCK_1D 64
#define CELL_SIZE 64
#define MAX_SHARED_TYPES 8
#define PERCENT_GPU_MEMORY 0.7
#define BIG_NUMBER 100000000

template <class numtyp, class acctyp>
class LJ_GPU_Memory {
 public:
  LJ_GPU_Memory() : allocated(false) {}
  ~LJ_GPU_Memory() { clear(); }
 
  /// Allocate memory on host and device
  bool init(const int ij_size, const int ntypes, double **host_cutsq, 
            double **host_sigma, double **host_epsilon, 
            double **host_lj1, double **host_lj2, double **host_lj3, 
            double **host_lj4, double **host_offset, double *host_special_lj,
            const int max_nbors, const int me);
  /// Free any memory on host and device
  void clear();

  /// Returns memory usage on GPU per atom
  int bytes_per_atom(const int max_nbors) const;
  /// Maximum number of atoms that can be stored on GPU
  int get_max_atoms(const size_t gpu_bytes, const int max_nbors);
  /// Total host memory used by library
  double host_memory_usage() const;
  
  // -------------------------   DATA   -----------------------------

  // Device Properties
  NVCDevice gpu;
  // Device Error Flag
  NVC_VecI dev_error;
  // Stream for asynchronous work
  cudaStream_t pair_stream;
  
  // Atom Data
  PairGPUAtom<numtyp,acctyp> atom;
  // Neighbor Data
  PairGPUNbor nbor;
  
  // --------------- Const Data for Atoms
  NVC_ConstMatT sigma, epsilon, cutsq, offset;
  NVC_ConstMat< typename nvc_vec_traits<numtyp>::vec2 > lj1, lj3;
  NVC_VecT special_lj;
  
  size_t max_atoms;
  
  // Timing for pair calculation
  NVCTimer time_pair;
  
  // If atom type constants fit in shared memory, use fast kernels
  bool shared_types;
   
 protected:
  bool allocated;
};

#endif
