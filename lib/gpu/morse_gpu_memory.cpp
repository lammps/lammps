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

#ifdef USE_OPENCL
#include "morse_gpu_cl.h"
#else
#include "morse_gpu_ptx.h"
#endif

#include "morse_gpu_memory.h"
#include <cassert>
#define MOR_GPU_MemoryT MOR_GPU_Memory<numtyp, acctyp>

extern PairGPUDevice<PRECISION,ACC_PRECISION> pair_gpu_device;

template <class numtyp, class acctyp>
MOR_GPU_MemoryT::MOR_GPU_Memory() : AtomicGPUMemory<numtyp,acctyp>(), _allocated(false) {
}

template <class numtyp, class acctyp>
MOR_GPU_MemoryT::~MOR_GPU_Memory() { 
  clear();
}
 
template <class numtyp, class acctyp>
int MOR_GPU_MemoryT::bytes_per_atom(const int max_nbors) const {
  return this->bytes_per_atom_atomic(max_nbors);
}

template <class numtyp, class acctyp>
int MOR_GPU_MemoryT::init(const int ntypes, 
                          double **host_cutsq, double **host_morse1, 
                          double **host_r0, double **host_alpha, 
                          double **host_d0, double **host_offset, 
                          double *host_special_lj, const int nlocal,
                          const int nall, const int max_nbors,
                          const int maxspecial, const double cell_size,
                          const double gpu_split, FILE *_screen) {
  int success;
  success=this->init_atomic(nlocal,nall,max_nbors,maxspecial,cell_size,gpu_split,
                            _screen,morse_gpu_kernel);
  if (success!=0)
    return success;

  // If atom type constants fit in shared memory use fast kernel
  int types=ntypes;
  shared_types=false;
  int max_shared_types=this->device->max_shared_types();
  if (types<=max_shared_types && this->_block_size>=max_shared_types) {
    types=max_shared_types;
    shared_types=true;
  }
  _types=types;

  // Allocate a host write buffer for data initialization
  UCL_H_Vec<numtyp> host_write(types*types*32,*(this->ucl_device),
                               UCL_WRITE_OPTIMIZED);

  for (int i=0; i<types*types; i++)
    host_write[i]=0.0;

  mor1.alloc(types*types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack4(ntypes,types,mor1,host_write,host_cutsq,host_morse1,
                         host_r0,host_alpha);

  mor2.alloc(types*types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack2(ntypes,types,mor2,host_write,host_d0,host_offset);

  UCL_H_Vec<double> dview;
  sp_lj.alloc(4,*(this->ucl_device),UCL_READ_ONLY);
  dview.view(host_special_lj,4,*(this->ucl_device));
  ucl_copy(sp_lj,dview,false);

  _allocated=true;
  this->_max_bytes=mor1.row_bytes()+mor2.row_bytes()+sp_lj.row_bytes();
  return 0;
}

template <class numtyp, class acctyp>
void MOR_GPU_MemoryT::clear() {
  if (!_allocated)
    return;
  _allocated=false;

  mor1.clear();
  mor2.clear();
  sp_lj.clear();
  this->clear_atomic();
}

template <class numtyp, class acctyp>
double MOR_GPU_MemoryT::host_memory_usage() const {
  return this->host_memory_usage_atomic()+sizeof(MOR_GPU_Memory<numtyp,acctyp>);
}

// ---------------------------------------------------------------------------
// Calculate energies, forces, and torques
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void MOR_GPU_MemoryT::loop(const bool _eflag, const bool _vflag) {
  // Compute the block size and grid size to keep all cores busy
  const int BX=this->block_size();
  int eflag, vflag;
  if (_eflag)
    eflag=1;
  else
    eflag=0;

  if (_vflag)
    vflag=1;
  else
    vflag=0;
  
  int GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/
                               (BX/this->_threads_per_atom)));

  int ainum=this->ans->inum();
  int anall=this->atom->nall();
  int nbor_pitch=this->nbor->nbor_pitch();
  this->time_pair.start();
  if (shared_types) {
    this->k_pair_fast.set_size(GX,BX);
    this->k_pair_fast.run(&this->atom->dev_x.begin(), &mor1.begin(),
                          &mor2.begin(), &sp_lj.begin(),
                          &this->nbor->dev_nbor.begin(),
                          &this->_nbor_data->begin(),
                          &this->ans->dev_ans.begin(),
                          &this->ans->dev_engv.begin(), &eflag, &vflag,
                          &ainum, &anall, &nbor_pitch, 
                          &this->_threads_per_atom);
  } else {
    this->k_pair.set_size(GX,BX);
    this->k_pair.run(&this->atom->dev_x.begin(), &mor1.begin(), &mor2.begin(),
                     &_types, &sp_lj.begin(), &this->nbor->dev_nbor.begin(),
                     &this->_nbor_data->begin(), &this->ans->dev_ans.begin(),
                     &this->ans->dev_engv.begin(), &eflag, &vflag, &ainum,
                     &anall, &nbor_pitch, &this->_threads_per_atom);
  }
  this->time_pair.stop();
}

template class MOR_GPU_Memory<PRECISION,ACC_PRECISION>;

