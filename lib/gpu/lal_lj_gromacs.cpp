/***************************************************************************
                                 lj_gromacs.cpp
                             -------------------
                            Trung Dac Nguyen (ORNL)

  Class for acceleration of the lj/gromacs pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : 
    email                : nguyentd@ornl.gov
 ***************************************************************************/

#if defined(USE_OPENCL)
#include "lj_gromacs_cl.h"
#elif defined(USE_CUDART)
const char *lj_gromacs=0;
#else
#include "lj_gromacs_cubin.h"
#endif

#include "lal_lj_gromacs.h"
#include <cassert>
using namespace LAMMPS_AL;
#define LJGROMACST LJGROMACS<numtyp, acctyp>

extern Device<PRECISION,ACC_PRECISION> device;

template <class numtyp, class acctyp>
LJGROMACST::LJGROMACS() : BaseAtomic<numtyp,acctyp>(),
                                    _allocated(false) {
}

template <class numtyp, class acctyp>
LJGROMACST::~LJGROMACS() {
  clear();
}
 
template <class numtyp, class acctyp>
int LJGROMACST::bytes_per_atom(const int max_nbors) const {
  return this->bytes_per_atom_atomic(max_nbors);
}

template <class numtyp, class acctyp>
int LJGROMACST::init(const int ntypes, double **host_cutsq,
                     double **host_lj1, double **host_lj2, double **host_lj3,
                     double **host_lj4, double *host_special_lj,
                     const int nlocal, const int nall, const int max_nbors, 
                     const int maxspecial, const double cell_size, 
                     const double gpu_split, FILE *_screen,
                     double **host_ljsw1, double **host_ljsw2, double **host_ljsw3,
                     double **host_ljsw4, double **host_ljsw5, 
                     double **cut_inner, double **cut_inner_sq) {
  int success;
  success=this->init_atomic(nlocal,nall,max_nbors,maxspecial,cell_size,gpu_split,
                            _screen,lj_gromacs,"k_lj_gromacs");
  if (success!=0)
    return success;

  // If atom type constants fit in shared memory use fast kernel
  int lj_types=ntypes;
  shared_types=false;
  int max_shared_types=this->device->max_shared_types();
  if (lj_types<=max_shared_types && this->_block_size>=max_shared_types) {
    lj_types=max_shared_types;
    shared_types=true;
  }
  _lj_types=lj_types;

  // Allocate a host write buffer for data initialization
  UCL_H_Vec<numtyp> host_write(lj_types*lj_types*32,*(this->ucl_device),
                               UCL_WRITE_ONLY);

  for (int i=0; i<lj_types*lj_types; i++)
    host_write[i]=0.0;

  lj1.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack4(ntypes,lj_types,lj1,host_write,host_lj1,host_lj2,
                         host_cutsq,cut_inner_sq);

  lj3.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack4(ntypes,lj_types,lj3,host_write,host_lj3,host_lj4,
                         cut_inner,host_ljsw5);

  ljsw.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack4(ntypes,lj_types,ljsw,host_write,host_ljsw1,host_ljsw2,
                         host_ljsw3,host_ljsw4);

  UCL_H_Vec<double> dview;
  sp_lj.alloc(4,*(this->ucl_device),UCL_READ_ONLY);
  dview.view(host_special_lj,4,*(this->ucl_device));
  ucl_copy(sp_lj,dview,false);

  _allocated=true;
  this->_max_bytes=lj1.row_bytes()+lj3.row_bytes()
    +ljsw.row_bytes()+sp_lj.row_bytes();
  return 0;
}

template <class numtyp, class acctyp>
void LJGROMACST::clear() {
  if (!_allocated)
    return;
  _allocated=false;

  lj1.clear();
  lj3.clear();
  ljsw.clear();
  sp_lj.clear();
  this->clear_atomic();
}

template <class numtyp, class acctyp>
double LJGROMACST::host_memory_usage() const {
  return this->host_memory_usage_atomic()+sizeof(LJGROMACS<numtyp,acctyp>);
}

// ---------------------------------------------------------------------------
// Calculate energies, forces, and torques
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void LJGROMACST::loop(const bool _eflag, const bool _vflag) {
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
  int nbor_pitch=this->nbor->nbor_pitch();
  this->time_pair.start();
  if (shared_types) {
    this->k_pair_fast.set_size(GX,BX);
    this->k_pair_fast.run(&this->atom->x, &lj1, &lj3, &ljsw,
                          &sp_lj, &this->nbor->dev_nbor,
                          &this->_nbor_data->begin(),
                          &this->ans->force, &this->ans->engv, 
                          &eflag, &vflag, &ainum, &nbor_pitch,
                          &this->_threads_per_atom);
  } else {
    this->k_pair.set_size(GX,BX);
    this->k_pair.run(&this->atom->x, &lj1, &lj3, &ljsw, &_lj_types, 
                     &sp_lj, &this->nbor->dev_nbor,
                     &this->_nbor_data->begin(), 
                     &this->ans->force, &this->ans->engv, 
                     &eflag, &vflag, &ainum, &nbor_pitch, 
                     &this->_threads_per_atom);
  }
  this->time_pair.stop();
}

template class LJGROMACS<PRECISION,ACC_PRECISION>;
