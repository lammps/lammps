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

#include "pair_gpu_nbor.h"

int PairGPUNbor::bytes_per_atom(const int max_nbors) const {
  if (_use_packing)
    return (max_nbors*2+4)*sizeof(int);
  else
    return (max_nbors+3)*sizeof(int);
}

bool PairGPUNbor::init(const int ij_size, const int max_atoms, 
                       const int max_nbors) {
  bool success=true;
  if (allocated)
    clear();
    
  // Initialize timers for the selected GPU
  time_nbor.init();

  if (_use_packing)
    success=success && dev_nbor.alloc((max_nbors+4)*max_atoms);
  else  
    success=success && dev_nbor.alloc(3*max_atoms);
  
  success=success && ij.alloc(max_nbors*max_atoms);
  success=success && host_ij.alloc_w(ij_size);
    
  allocated=true;
  
  return success;
}
  
void PairGPUNbor::resize(const int nlocal, const int max_nbor, bool &success) {
  dev_nbor.clear();
  ij.clear();
  if (_use_packing)
    success=success && dev_nbor.alloc((max_nbor+4)*nlocal);
  else  
    success=success && dev_nbor.alloc(3*nlocal);
  success=success && ij.alloc(max_nbor*nlocal);
  allocated=true;
}

void PairGPUNbor::clear() {
  if (!allocated)
    return;
  allocated=false;

  ij.clear();
  host_ij.clear();
  dev_nbor.clear();
}  

double PairGPUNbor::host_memory_usage() const {
  return IJ_SIZE*sizeof(int)+sizeof(PairGPUNbor);
}

void PairGPUNbor::reset(const int inum, int *ilist, const int *numj, 
                        cudaStream_t &s) {  
  ij_total=0;

  dev_nbor.copy_from_host(ilist,inum);
  int acc=0;
   
  int ij_size=host_ij.numel();
  if (inum*2<ij_size) {
    for (int i=0; i<inum; i++) {
      host_ij[i]=numj[ilist[i]];
      host_ij[i+inum]=acc;
      acc+=numj[ilist[i]];
    }
    host_ij.copy_to_device(dev_nbor.begin()+inum,2*inum, s);
  } else {
    int offset=0;
    int half=ij_size/2;
    int hi=0;
    for (int i=0; i<inum; i++) {
      host_ij[hi]=numj[ilist[i]];
      host_ij[hi+half]=acc;
      acc+=numj[ilist[i]];
      hi++;
      if (hi==half) {
        host_ij.copy_to_device(dev_nbor.begin()+inum+offset,half,s);
        host_ij.copy_to_device(half,dev_nbor.begin()+2*inum+offset,half,s);
        offset+=half;
        hi=0;
        CUDA_SAFE_CALL(cudaStreamSynchronize(s));
      }
   }
   if (hi>0) {
     host_ij.copy_to_device(dev_nbor.begin()+inum+offset,hi,s);
     host_ij.copy_to_device(half,dev_nbor.begin()+2*inum+offset,hi,s);
   }
 }
}
