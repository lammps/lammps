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

void PairGPUNbor::init(const int ij_size, const int max_atoms, 
                       const int max_nbors) {
  if (allocated)
    clear();
    
  // Initialize timers for the selected GPU
  time_nbor.init();

  if (_use_packing)
    dev_nbor.safe_alloc(max_nbors+4,max_atoms);
  else  
    dev_nbor.safe_alloc(3,max_atoms);
  
  ij.safe_alloc(max_nbors*max_atoms);
  host_ij.safe_alloc_w(ij_size);
    
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
  for (int i=0; i<inum; i++) {
    host_ij[i]=numj[ilist[i]];
    host_ij[i+inum]=acc;
    acc+=numj[ilist[i]];
  }
  
  host_ij.copy_to_2Ddevice(dev_nbor.begin()+dev_nbor.row_size(),
                           dev_nbor.row_size(),2,inum, s);
}
