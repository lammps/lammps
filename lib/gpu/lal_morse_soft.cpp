/***************************************************************************
                                 morse_soft.cpp
                             -------------------

  Class for acceleration of the morse/soft pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : October 2026
 ***************************************************************************/

#if defined(USE_OPENCL)
#include "morse_soft_cl.h"
#elif defined(USE_CUDART)
const char *morse_soft=0;
#else
#include "morse_soft_cubin.h"
#endif

#include "lal_morse_soft.h"
#include <cassert>
namespace LAMMPS_AL {
#define MorseSoftT MorseSoft<numtyp, acctyp>

extern Device<PRECISION,ACC_PRECISION> device;

template <class numtyp, class acctyp>
MorseSoftT::MorseSoft() : BaseAtomic<numtyp,acctyp>(), _allocated(false) {
}

template <class numtyp, class acctyp>
MorseSoftT::~MorseSoft() {
  clear();
}

template <class numtyp, class acctyp>
int MorseSoftT::bytes_per_atom(const int max_nbors) const {
  return this->bytes_per_atom_atomic(max_nbors);
}

template <class numtyp, class acctyp>
int MorseSoftT::init(const int ntypes,
                          double **host_cutsq, double **host_morse1,
                          double **host_r0, double **host_alpha,
                          double **host_d0, double **host_offset,
                          double **host_bs, double **host_vs,
                          double *host_special_lj, const int nlocal,
                          const int nall, const int max_nbors,
                          const int maxspecial, const double cell_size,
                           FILE *_screen) {
  int success;
  success=this->init_atomic(nlocal,nall,max_nbors,maxspecial,cell_size,
                            _screen,morse_soft,"k_morse_soft");
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
                               UCL_WRITE_ONLY);

  for (int i=0; i<types*types; i++)
    host_write[i]=0.0;

  mor1.alloc(types*types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack4(ntypes,types,mor1,host_write,host_cutsq,host_morse1,
                         host_r0,host_alpha);

  mor2.alloc(types*types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack4(ntypes,types,mor2,host_write,host_d0,host_offset,
                         host_bs,host_vs);

  UCL_H_Vec<double> dview;
  sp_lj.alloc(4,*(this->ucl_device),UCL_READ_ONLY);
  dview.view(host_special_lj,4,*(this->ucl_device));
  ucl_copy(sp_lj,dview,false);

  _allocated=true;
  this->_max_bytes=mor1.row_bytes()+mor2.row_bytes()+sp_lj.row_bytes();
  return 0;
}

template <class numtyp, class acctyp>
void MorseSoftT::reinit(const int ntypes, double **host_cutsq, double **host_morse1,
                        double **host_r0, double **host_alpha, double **host_d0,
                        double **host_offset, double **host_bs, double **host_vs) {
  // Allocate a host write buffer for data initialization
  UCL_H_Vec<numtyp> host_write(_types*_types*32,*(this->ucl_device),
                               UCL_WRITE_ONLY);

  for (int i=0; i<_types*_types; i++)
    host_write[i]=0.0;

  this->atom->type_pack4(ntypes,_types,mor1,host_write,host_cutsq,host_morse1,
                         host_r0,host_alpha);
  this->atom->type_pack4(ntypes,_types,mor2,host_write,host_d0,host_offset,
                         host_bs,host_vs);
}

template <class numtyp, class acctyp>
void MorseSoftT::clear() {
  if (!_allocated)
    return;
  _allocated=false;

  mor1.clear();
  mor2.clear();
  sp_lj.clear();
  this->clear_atomic();
}

template <class numtyp, class acctyp>
double MorseSoftT::host_memory_usage() const {
  return this->host_memory_usage_atomic()+sizeof(MorseSoft<numtyp,acctyp>);
}

// ---------------------------------------------------------------------------
// Calculate energies, forces, and torques
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int MorseSoftT::loop(const int eflag, const int vflag) {
  // Compute the block size and grid size to keep all cores busy
  const int BX=this->block_size();
  int GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/
                               (BX/this->_threads_per_atom)));

  int ainum=this->ans->inum();
  int nbor_pitch=this->nbor->nbor_pitch();
  this->time_pair.start();
  if (shared_types) {
    this->k_pair_sel->set_size(GX,BX);
    this->k_pair_sel->run(&this->atom->x, &mor1, &mor2, &sp_lj,
                          &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                          &this->ans->force, &this->ans->engv, &eflag,
                          &vflag, &ainum, &nbor_pitch,
                          &this->_threads_per_atom);
  } else {
    this->k_pair.set_size(GX,BX);
    this->k_pair.run(&this->atom->x, &mor1, &mor2, &_types, &sp_lj,
                     &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                     &this->ans->force, &this->ans->engv, &eflag, &vflag,
                     &ainum, &nbor_pitch, &this->_threads_per_atom);
  }
  this->time_pair.stop();
  return GX;
}

template class MorseSoft<PRECISION,ACC_PRECISION>;
}
