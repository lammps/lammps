/***************************************************************************
                                 coul_cut_soft.cpp
                             -------------------
     Claude Opus 5.5 (Anthropic), under the direction of Axel Kohlmeyer

  Class for acceleration of the coul/cut/soft pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : October 2026
 ***************************************************************************/

#if defined(USE_OPENCL)
#include "coul_cut_soft_cl.h"
#elif defined(USE_CUDART)
const char *coul_cut_soft=0;
#else
#include "coul_cut_soft_cubin.h"
#endif

#include "lal_coul_cut_soft.h"
#include <cassert>
namespace LAMMPS_AL {
#define CoulCutSoftT CoulCutSoft<numtyp, acctyp>

extern Device<PRECISION,ACC_PRECISION> device;

template <class numtyp, class acctyp>
CoulCutSoftT::CoulCutSoft() : BaseCharge<numtyp,acctyp>(),
  _allocated(false) {
}

template <class numtyp, class acctyp>
CoulCutSoftT::~CoulCutSoft() {
  clear();
}

template <class numtyp, class acctyp>
int CoulCutSoftT::bytes_per_atom(const int max_nbors) const {
  return this->bytes_per_atom_atomic(max_nbors);
}

template <class numtyp, class acctyp>
int CoulCutSoftT::init(const int ntypes, double **host_lam1, double **host_lam2, double **host_cutsq,
                double *host_special_coul, const int nlocal,
                const int nall, const int max_nbors,
                const int maxspecial, const double cell_size,
                           FILE *_screen,
                const double qqrd2e) {
  int success;
  success=this->init_atomic(nlocal,nall,max_nbors,maxspecial,cell_size,
                            _screen,coul_cut_soft,"k_coul_cut_soft");
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

  lam.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack2(ntypes,lj_types,lam,host_write,host_lam1,host_lam2);

  cutsq.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack1(ntypes,lj_types,cutsq,host_write,host_cutsq);

  sp_cl.alloc(4,*(this->ucl_device),UCL_READ_ONLY);
  for (int i=0; i<4; i++) {
    host_write[i]=host_special_coul[i];
  }
  ucl_copy(sp_cl,host_write,4,false);

  _qqrd2e=qqrd2e;

  _allocated=true;
  this->_max_bytes=lam.row_bytes()+cutsq.row_bytes()+sp_cl.row_bytes();
  return 0;
}

template <class numtyp, class acctyp>
void CoulCutSoftT::reinit(const int ntypes, double **host_lam1, double **host_lam2) {
  // Allocate a host write buffer for data initialization
  UCL_H_Vec<numtyp> host_write(_lj_types*_lj_types*32,*(this->ucl_device),
                               UCL_WRITE_ONLY);

  for (int i=0; i<_lj_types*_lj_types; i++)
    host_write[i]=0.0;

  this->atom->type_pack2(ntypes,_lj_types,lam,host_write,host_lam1,host_lam2);
}

template <class numtyp, class acctyp>
void CoulCutSoftT::clear() {
  if (!_allocated)
    return;
  _allocated=false;

  lam.clear();
  cutsq.clear();
  sp_cl.clear();
  this->clear_atomic();
}

template <class numtyp, class acctyp>
double CoulCutSoftT::host_memory_usage() const {
  return this->host_memory_usage_atomic()+sizeof(CoulCutSoft<numtyp,acctyp>);
}

// ---------------------------------------------------------------------------
// Calculate energies, forces, and torques
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int CoulCutSoftT::loop(const int eflag, const int vflag) {
  // Compute the block size and grid size to keep all cores busy
  const int BX=this->block_size();
  int GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/
                               (BX/this->_threads_per_atom)));

  int ainum=this->ans->inum();
  int nbor_pitch=this->nbor->nbor_pitch();
  this->time_pair.start();
  if (shared_types) {
    this->k_pair_sel->set_size(GX,BX);
    this->k_pair_sel->run(&this->atom->x, &lam, &sp_cl,
                          &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                          &this->ans->force, &this->ans->engv, &eflag,
                          &vflag, &ainum, &nbor_pitch, &this->atom->q,
                          &cutsq, &_qqrd2e, &this->_threads_per_atom);
  } else {
    this->k_pair.set_size(GX,BX);
    this->k_pair.run(&this->atom->x, &lam, &_lj_types, &sp_cl,
                     &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                     &this->ans->force, &this->ans->engv,
                     &eflag, &vflag, &ainum, &nbor_pitch, &this->atom->q,
                     &cutsq, &_qqrd2e, &this->_threads_per_atom);
  }
  this->time_pair.stop();
  return GX;
}

template class CoulCutSoft<PRECISION,ACC_PRECISION>;
}
