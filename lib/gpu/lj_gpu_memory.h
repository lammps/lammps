/***************************************************************************
                               lj_gpu_memory.h
                             -------------------
                               W. Michael Brown

  Global variables for GPU Lennard-Jones Library

 __________________________________________________________________________
    This file is part of the LAMMPS GPU Library
 __________________________________________________________________________

    begin                : Tue Aug 4 2009
    copyright            : (C) 2009 by W. Michael Brown
    email                : wmbrown@sandia.gov
 ***************************************************************************/

/* -----------------------------------------------------------------------
   Copyright (2009) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.
   ----------------------------------------------------------------------- */

#ifndef LJ_GPU_MEMORY_H
#define LJ_GPU_MEMORY_H

#include "nvc_device.h"
#include "nvc_traits.h"
#include "pair_gpu_atom.h"
#include "pair_gpu_nbor.h"

#define BLOCK_1D 64
#define MAX_SHARED_TYPES 8
#define PERCENT_GPU_MEMORY 0.7

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
