/***************************************************************************
                               pair_gpu_nbor.h
                             -------------------
                               W. Michael Brown

  Neighbor memory operations for LAMMPS GPU Library

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

#ifndef PAIR_GPU_NBOR_H
#define PAIR_GPU_NBOR_H

#include "nvc_macros.h"
#include "nvc_timer.h"
#include "nvc_memory.h"

#define BLOCK_1D 64
#define IJ_SIZE 131072

class PairGPUNbor {
 public:
  PairGPUNbor() : _use_packing(false), allocated(false) {}
  ~PairGPUNbor() { clear(); }
 
  /// Determine whether neighbor packing should be used
  /** If true, twice as much memory is reserved to allow packing neighbors by 
    * atom for coalesced access after cutoff evaluation. This can be used
    * for expensive potentials where it is more efficient to evaluate the 
    * cutoff separately from the potential in order to reduce thread divergence 
    * for expensive routines **/
  void packing(const bool use_packing) { _use_packing=use_packing; }
  
  /// Called once to allocate memory
  void init(const int ij_size, const int max_atoms, const int max_nbors);
  
  /// Free all memory on host and device
  void clear();
 
  /// Bytes per atom used on device
  int bytes_per_atom(const int max_nbors) const;
  /// Total host memory used by class
  double host_memory_usage() const;

  /// Reset neighbor data (first time or from a rebuild)  
  void reset(const int inum, int *ilist, const int *numj, cudaStream_t &s);
  /// Add neighbor data from host
  inline void add(const int num_ij, cudaStream_t &s)
    { host_ij.copy_to_device(ij.begin()+ij_total,num_ij,s); ij_total+=num_ij; }

  /// Pack neighbors satisfying cutoff by atom for coalesced access
  void pack_nbors(const int GX, const int BX, const int start, 
                  const int inum, const int form_low, const int form_high);

    
  // ------------------------------- Data -------------------------------

  // Store IJ interactions on device
  NVC_VecI ij;
  // Buffer for moving ij data to GPU
  NVC_HostI host_ij;

  // --------------- Atom neighbors
  // 3 x n
  // - 1st row is i
  // - 2nd row is numj (number of neighbors)
  // - 3rd row is starting address in host_ij of neighbors
  NVC_MatI dev_nbor;

  // --------------- Timing Stuff
  NVCTimer time_nbor;
  
  int ij_total;
 private:
  bool allocated, _use_packing;
};

#endif
