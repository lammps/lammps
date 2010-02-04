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

#include "lj_gpu_memory.h"
#define LJ_GPU_MemoryT LJ_GPU_Memory<numtyp, acctyp>

template <class numtyp, class acctyp>
int LJ_GPU_MemoryT::bytes_per_atom(const int max_nbors) const {
  return atom.bytes_per_atom()+nbor.bytes_per_atom(max_nbors); 
}

template <class numtyp, class acctyp>
int LJ_GPU_MemoryT::get_max_atoms(const size_t gpu_bytes, const int max_nbors) {
  int matoms=static_cast<int>(PERCENT_GPU_MEMORY*gpu_bytes/
                              bytes_per_atom(max_nbors));
  if (matoms>MAX_ATOMS)
    matoms=MAX_ATOMS;
  return matoms;
}
  
template <class numtyp, class acctyp>
bool LJ_GPU_MemoryT::init(const int ij_size, const int ntypes, 
                          double **host_cutsq, double **host_sigma, 
                          double **host_epsilon, double **host_lj1, 
                          double **host_lj2, double **host_lj3, 
                          double **host_lj4, double **host_offset, 
                          double *host_special_lj, const int max_nbors, 
                          const int me) {
  if (allocated)
    clear();
    
  if (me>=gpu.num_devices())
    return false;
  gpu.set(me);
  if (gpu.revision()<1.0)
    return false;  
    
  // Initialize timers for the selected GPU
  time_pair.init();

  // Initialize atom and nbor data
  max_atoms=get_max_atoms(gpu.bytes(),max_nbors);
  atom.init(max_atoms);
  nbor.init(ij_size,max_atoms,max_nbors);
  
  // Get a stream for computing pair potentials
  CUDA_SAFE_CALL(cudaStreamCreate(&pair_stream));
    
  // Use the write buffer from atom for data initialization
  NVC_HostT &host_write=atom.host_write;
  assert(host_write.numel()>4 && host_write.numel()>ntypes*ntypes*2);

  // Copy data for bonded interactions
  special_lj.safe_alloc(4);
  special_lj.cast_copy(host_special_lj,host_write);

  // Copy sigma, epsilon, and cutsq onto GPU
  sigma.safe_alloc(ntypes,ntypes,sigma_get_texture<numtyp>());
  sigma.cast_copy(host_sigma[0],host_write);
  epsilon.safe_alloc(ntypes,ntypes,epsilon_get_texture<numtyp>());
  epsilon.cast_copy(host_epsilon[0],host_write);
  cutsq.safe_alloc(ntypes,ntypes,cutsq_get_texture<numtyp>());
  cutsq.cast_copy(host_cutsq[0],host_write);

  // If atom type constants fit in shared memory use fast kernel
  int lj_types=ntypes;
  shared_types=false;
  if (lj_types<=MAX_SHARED_TYPES) {
    lj_types=MAX_SHARED_TYPES;
    shared_types=true;
  }
  offset.safe_alloc(lj_types,lj_types,offset_get_texture<numtyp>());
  offset.cast_copy2D(host_offset[0],host_write,ntypes,ntypes);
  double *t1=host_lj1[0];
  double *t2=host_lj2[0];
  for (int i=0; i<ntypes*ntypes; i++) {
    host_write[i*2]=t1[i];
    host_write[i*2+1]=t2[i];
  }
  lj1.safe_alloc(lj_types,lj_types,lj1_get_texture<numtyp>());
  lj1.copy_2Dfrom_host(reinterpret_cast<typename nvc_vec_traits<numtyp>::vec2 *> (host_write.begin()),
                       ntypes,ntypes);
  t1=host_lj3[0];
  t2=host_lj4[0];
  for (int i=0; i<ntypes*ntypes; i++) {
    host_write[i*2]=t1[i];
    host_write[i*2+1]=t2[i];
  }
  lj3.safe_alloc(lj_types,lj_types,lj3_get_texture<numtyp>());
  lj3.copy_2Dfrom_host(reinterpret_cast<typename nvc_vec_traits<numtyp>::vec2 *> (host_write.begin()),
                       ntypes,ntypes);
        
  dev_error.safe_alloc(1);
  dev_error.zero();
    
  allocated=true;
  return true;
}
  
template <class numtyp, class acctyp>
void LJ_GPU_MemoryT::clear() {
  if (!allocated)
    return;
  allocated=false;
      
  // Check for any pair style specific errors here
  int err_flag;
  dev_error.copy_to_host(&err_flag);
 
  atom.clear();
  nbor.clear();
    
  CUDA_SAFE_CALL(cudaStreamDestroy(pair_stream));

  dev_error.clear();
  sigma.clear();
  epsilon.clear();
  special_lj.clear();
  cutsq.clear();
  offset.clear();
  lj1.clear();
  lj3.clear();
}  
 
template <class numtyp, class acctyp>
double LJ_GPU_MemoryT::host_memory_usage() const {
  return atom.host_memory_usage(max_atoms)+nbor.host_memory_usage()+
         sizeof(LJ_GPU_Memory<numtyp,acctyp>);
}

template class LJ_GPU_Memory<PRECISION,ACC_PRECISION>;
