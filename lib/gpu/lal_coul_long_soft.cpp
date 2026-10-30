/***************************************************************************
                                 coul_long_soft.cpp
                             -------------------

  Class for acceleration of the coul/long/soft pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : October 2026
 ***************************************************************************/

#if defined(USE_OPENCL)
#include "coul_long_soft_cl.h"
#elif defined(USE_CUDART)
const char *coul_long_soft=0;
#else
#include "coul_long_soft_cubin.h"
#endif

#include "lal_coul_long_soft.h"
#include <cassert>
namespace LAMMPS_AL {
#define CoulLongSoftT CoulLongSoft<numtyp, acctyp>

extern Device<PRECISION,ACC_PRECISION> pair_gpu_device;

template <class numtyp, class acctyp>
CoulLongSoftT::CoulLongSoft() : BaseCharge<numtyp,acctyp>(), _allocated(false) {
}

template <class numtyp, class acctyp>
CoulLongSoftT::~CoulLongSoft() {
  clear();
}

template <class numtyp, class acctyp>
int CoulLongSoftT::bytes_per_atom(const int max_nbors) const {
  return this->bytes_per_atom_atomic(max_nbors);
}

template <class numtyp, class acctyp>
int CoulLongSoftT::init(const int ntypes, double **host_scale, double **host_lam1, double **host_lam2,
                    const int nlocal, const int nall, const int max_nbors,
                    const int maxspecial, const double cell_size,
                           FILE *_screen,
                    const double host_cut_coulsq, double *host_special_coul,
                    const double qqrd2e, const double g_ewald) {
  int success;
  success=this->init_atomic(nlocal,nall,max_nbors,maxspecial,cell_size,
                            _screen,coul_long_soft,"k_coul_long_soft");
  if (success!=0)
    return success;

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

  scl.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack4(ntypes,lj_types,scl,host_write,host_scale,host_lam1,host_lam2);

  sp_cl.alloc(4,*(this->ucl_device),UCL_READ_ONLY);
  for (int i=0; i<4; i++) {
    host_write[i]=host_special_coul[i];
  }
  ucl_copy(sp_cl,host_write,4,false);

  _cut_coulsq=host_cut_coulsq;
  _qqrd2e=qqrd2e;
  _g_ewald=g_ewald;

  _allocated=true;
  this->_max_bytes=scl.row_bytes()+sp_cl.row_bytes();
  return 0;
}

template <class numtyp, class acctyp>
void CoulLongSoftT::reinit(const int ntypes, double **host_scale, double **host_lam1,
                           double **host_lam2) {
  UCL_H_Vec<numtyp> hscale(_lj_types*_lj_types*4,*(this->ucl_device),
                           UCL_WRITE_ONLY);
  this->atom->type_pack4(ntypes,_lj_types,scl,hscale,host_scale,host_lam1,host_lam2);
}

template <class numtyp, class acctyp>
void CoulLongSoftT::clear() {
  if (!_allocated)
    return;
  _allocated=false;

  scl.clear();
  sp_cl.clear();
  this->clear_atomic();
}

template <class numtyp, class acctyp>
double CoulLongSoftT::host_memory_usage() const {
  return this->host_memory_usage_atomic()+sizeof(CoulLongSoft<numtyp,acctyp>);
}

// ---------------------------------------------------------------------------
// Calculate energies, forces, and torques
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int CoulLongSoftT::loop(const int eflag, const int vflag) {
  // Compute the block size and grid size to keep all cores busy
  const int BX=this->block_size();
  int GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/
                               (BX/this->_threads_per_atom)));

  int ainum=this->ans->inum();
  int nbor_pitch=this->nbor->nbor_pitch();
  this->time_pair.start();
  if (shared_types) {
    this->k_pair_sel->set_size(GX,BX);
    this->k_pair_sel->run(&this->atom->x, &scl, &sp_cl,
                          &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                          &this->ans->force, &this->ans->engv,
                          &eflag, &vflag, &ainum, &nbor_pitch,
                          &this->atom->q, &_cut_coulsq, &_qqrd2e, &_g_ewald,
                          &this->_threads_per_atom);
  } else {
    this->k_pair.set_size(GX,BX);
    this->k_pair.run(&this->atom->x, &scl, &_lj_types, &sp_cl,
                     &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                     &this->ans->force, &this->ans->engv, &eflag, &vflag,
                     &ainum, &nbor_pitch, &this->atom->q, &_cut_coulsq,
                     &_qqrd2e, &_g_ewald, &this->_threads_per_atom);
  }
  this->time_pair.stop();
  return GX;
}

template class CoulLongSoft<PRECISION,ACC_PRECISION>;
}
