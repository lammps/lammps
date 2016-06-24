/***************************************************************************
                                 coul.cpp
                             -------------------
                               Trung Dac Nguyen

  Class for acceleration of the coul/cut pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : ndtrung@umich.edu
 ***************************************************************************/

#if defined(USE_OPENCL)
#include "coul_cl.h"
#elif defined(USE_CUDART)
const char *coul=0;
#else
#include "coul_cubin.h"
#endif

#include "lal_coul.h"
#include <cassert>
using namespace LAMMPS_AL;
#define CoulT Coul<numtyp, acctyp>

extern Device<PRECISION,ACC_PRECISION> device;

template <class numtyp, class acctyp>
CoulT::Coul() : BaseCharge<numtyp,acctyp>(),
  _allocated(false) {
}

template <class numtyp, class acctyp>
CoulT::~Coul() {
  clear();
}

template <class numtyp, class acctyp>
int CoulT::bytes_per_atom(const int max_nbors) const {
  return this->bytes_per_atom_atomic(max_nbors);
}

template <class numtyp, class acctyp>
int CoulT::init(const int ntypes, double **host_scale, double **host_cutsq,
                double *host_special_coul, const int nlocal,
                const int nall, const int max_nbors,
                const int maxspecial, const double cell_size,
                const double gpu_split, FILE *_screen,
                const double qqrd2e) {
  int success;
  success=this->init_atomic(nlocal,nall,max_nbors,maxspecial,cell_size,gpu_split,
                            _screen,coul,"k_coul");
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

  scale.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack1(ntypes,lj_types,scale,host_write,host_scale);

  cutsq.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack1(ntypes,lj_types,cutsq,host_write,host_cutsq);

  sp_cl.alloc(4,*(this->ucl_device),UCL_READ_ONLY);
  for (int i=0; i<4; i++) {
    host_write[i]=host_special_coul[i];
  }
  ucl_copy(sp_cl,host_write,4,false);

  _qqrd2e=qqrd2e;

  _allocated=true;
  this->_max_bytes=cutsq.row_bytes()+sp_cl.row_bytes();
  return 0;
}

template <class numtyp, class acctyp>
void CoulT::reinit(const int ntypes, double **host_scale) {
  // Allocate a host write buffer for data initialization
  UCL_H_Vec<numtyp> host_write(_lj_types*_lj_types*32,*(this->ucl_device),
                               UCL_WRITE_ONLY);

  for (int i=0; i<_lj_types*_lj_types; i++)
    host_write[i]=0.0;

  this->atom->type_pack1(ntypes,_lj_types,scale,host_write,host_scale);
}

template <class numtyp, class acctyp>
void CoulT::clear() {
  if (!_allocated)
    return;
  _allocated=false;

  scale.clear();
  cutsq.clear();
  sp_cl.clear();
  this->clear_atomic();
}

template <class numtyp, class acctyp>
double CoulT::host_memory_usage() const {
  return this->host_memory_usage_atomic()+sizeof(Coul<numtyp,acctyp>);
}

// ---------------------------------------------------------------------------
// Calculate energies, forces, and torques
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void CoulT::loop(const bool _eflag, const bool _vflag) {
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
    this->k_pair_fast.run(&this->atom->x, &scale, &sp_cl,
                          &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                          &this->ans->force, &this->ans->engv, &eflag,
                          &vflag, &ainum, &nbor_pitch, &this->atom->q,
                          &cutsq, &_qqrd2e, &this->_threads_per_atom);
  } else {
    this->k_pair.set_size(GX,BX);
    this->k_pair.run(&this->atom->x, &scale, &_lj_types, &sp_cl,
                     &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                     &this->ans->force, &this->ans->engv,
                     &eflag, &vflag, &ainum, &nbor_pitch, &this->atom->q,
                     &cutsq, &_qqrd2e, &this->_threads_per_atom);
  }
  this->time_pair.stop();
}

template class Coul<PRECISION,ACC_PRECISION>;
