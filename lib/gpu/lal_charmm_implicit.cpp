/***************************************************************************
                                 charmm_implicit.cpp
                             -------------------

  Class for acceleration of the lj/charmm/coul/charmm/implicit pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : October 2026
 ***************************************************************************/

#if defined(USE_OPENCL)
#include "charmm_implicit_cl.h"
#elif defined(USE_CUDART)
const char *charmm_implicit=0;
#else
#include "charmm_implicit_cubin.h"
#endif

#include "lal_charmm_implicit.h"
#include <cassert>
namespace LAMMPS_AL {
#define CHARMMImplicitT CHARMMImplicit<numtyp, acctyp>

extern Device<PRECISION,ACC_PRECISION> device;

template <class numtyp, class acctyp>
CHARMMImplicitT::CHARMMImplicit() : BaseCharge<numtyp,acctyp>(), _allocated(false) {
}

template <class numtyp, class acctyp>
CHARMMImplicitT::~CHARMMImplicit() {
  clear();
}

template <class numtyp, class acctyp>
int CHARMMImplicitT::bytes_per_atom(const int max_nbors) const {
  return this->bytes_per_atom_atomic(max_nbors);
}

template <class numtyp, class acctyp>
int CHARMMImplicitT::init(const int ntypes, double **host_lj1, double **host_lj2, double **host_lj3,
                          double **host_lj4, double *host_special_lj, double *host_special_coul,
                          const int nlocal, const int nall, const int max_nbors,
                          const int maxspecial, const double cell_size, FILE *_screen,
                          const double host_qqrd2e, const double host_cut_bothsq,
                          const double host_cut_ljsq, const double host_cut_lj_innersq,
                          const double host_cut_coulsq, const double host_cut_coul_innersq,
                          const double host_inv_denom_lj, const double host_inv_denom_coul) {
  int success;
  success=this->init_atomic(nlocal, nall, max_nbors, maxspecial, cell_size, _screen,
                            charmm_implicit, "k_charmm_implicit");
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
  UCL_H_Vec<numtyp> host_write(lj_types*lj_types*32,*(this->ucl_device), UCL_WRITE_ONLY);

  for (int i=0; i<lj_types*lj_types; i++)
    host_write[i]=0.0;

  lj1.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack4(ntypes, lj_types, lj1, host_write, host_lj1, host_lj2, host_lj3, host_lj4);

  sp_lj.alloc(8,*(this->ucl_device),UCL_READ_ONLY);
  for (int i=0; i<4; i++) {
    host_write[i]=host_special_lj[i];
    host_write[i+4]=host_special_coul[i];
  }
  ucl_copy(sp_lj,host_write,8,false);

  _qqrd2e=host_qqrd2e;
  _cut_bothsq=host_cut_bothsq;
  _cut_ljsq=host_cut_ljsq;
  _cut_lj_innersq=host_cut_lj_innersq;
  _cut_coulsq=host_cut_coulsq;
  _cut_coul_innersq=host_cut_coul_innersq;
  _inv_denom_lj=host_inv_denom_lj;
  _inv_denom_coul=host_inv_denom_coul;

  _allocated=true;
  this->_max_bytes=lj1.row_bytes()+sp_lj.row_bytes();
  return 0;
}

template <class numtyp, class acctyp>
void CHARMMImplicitT::reinit(const int ntypes, double **host_lj1, double **host_lj2,
                             double **host_lj3, double **host_lj4) {
  // Allocate a host write buffer for data initialization
  UCL_H_Vec<numtyp> host_write(_lj_types*_lj_types*32,*(this->ucl_device), UCL_WRITE_ONLY);

  for (int i=0; i<_lj_types*_lj_types; i++)
    host_write[i]=0.0;

  this->atom->type_pack4(ntypes, _lj_types, lj1, host_write, host_lj1, host_lj2, host_lj3,
                         host_lj4);
}

template <class numtyp, class acctyp>
void CHARMMImplicitT::clear() {
  if (!_allocated)
    return;
  _allocated=false;

  lj1.clear();
  sp_lj.clear();
  this->clear_atomic();
}

template <class numtyp, class acctyp>
double CHARMMImplicitT::host_memory_usage() const {
  return this->host_memory_usage_atomic()+sizeof(CHARMMImplicit<numtyp,acctyp>);
}

// ---------------------------------------------------------------------------
// Calculate energies, forces, and torques
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int CHARMMImplicitT::loop(const int eflag, const int vflag) {
  // Compute the block size and grid size to keep all cores busy
  const int BX=this->block_size();
  int GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/
                               (BX/this->_threads_per_atom)));

  int ainum=this->ans->inum();
  int nbor_pitch=this->nbor->nbor_pitch();
  this->time_pair.start();
  if (shared_types) {
    this->k_pair_sel->set_size(GX,BX);
    this->k_pair_sel->run(&this->atom->x, &lj1, &sp_lj, &this->nbor->dev_nbor,
                          &this->_nbor_data->begin(), &this->ans->force, &this->ans->engv, &eflag,
                          &vflag, &ainum, &nbor_pitch, &this->atom->q, &_qqrd2e, &_cut_bothsq,
                          &_cut_ljsq, &_cut_lj_innersq, &_cut_coulsq, &_cut_coul_innersq,
                          &_inv_denom_lj, &_inv_denom_coul, &this->_threads_per_atom);
  } else {
    this->k_pair.set_size(GX,BX);
    this->k_pair.run(&this->atom->x, &lj1, &_lj_types, &sp_lj, &this->nbor->dev_nbor,
                     &this->_nbor_data->begin(), &this->ans->force, &this->ans->engv, &eflag,
                     &vflag, &ainum, &nbor_pitch, &this->atom->q, &_qqrd2e, &_cut_bothsq,
                     &_cut_ljsq, &_cut_lj_innersq, &_cut_coulsq, &_cut_coul_innersq, &_inv_denom_lj,
                     &_inv_denom_coul, &this->_threads_per_atom);
  }
  this->time_pair.stop();
  return GX;
}

template class CHARMMImplicit<PRECISION,ACC_PRECISION>;
}
