/***************************************************************************
                                 nm_cut_coul_long.cpp
                             -------------------

  Class for acceleration of the nm/cut/coul/long pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : October 2026
 ***************************************************************************/

#if defined(USE_OPENCL)
#include "nm_cut_coul_long_cl.h"
#elif defined(USE_CUDART)
const char *nm_cut_coul_long=0;
#else
#include "nm_cut_coul_long_cubin.h"
#endif

#include "lal_nm_cut_coul_long.h"
#include <cassert>
namespace LAMMPS_AL {
#define NMCutCoulLongT NMCutCoulLong<numtyp, acctyp>

extern Device<PRECISION,ACC_PRECISION> device;

template <class numtyp, class acctyp>
NMCutCoulLongT::NMCutCoulLong() : BaseCharge<numtyp,acctyp>(), _allocated(false) {
}

template <class numtyp, class acctyp>
NMCutCoulLongT::~NMCutCoulLong() {
  clear();
}

template <class numtyp, class acctyp>
int NMCutCoulLongT::bytes_per_atom(const int max_nbors) const {
  return this->bytes_per_atom_atomic(max_nbors);
}

template <class numtyp, class acctyp>
int NMCutCoulLongT::init(const int ntypes, double **host_e0nm, double **host_nm, double **host_r0n,
                         double **host_r0m, double **host_nn, double **host_mm,
                         double **host_cut_ljsq, double **host_offset, double **host_cutsq,
                         double *host_special_lj, double *host_special_coul, const int nlocal,
                         const int nall, const int max_nbors, const int maxspecial,
                         const double cell_size, FILE *_screen, const double host_cut_coulsq,
                         const double host_qqrd2e, const double host_g_ewald) {
  int success;
  success=this->init_atomic(nlocal, nall, max_nbors, maxspecial, cell_size, _screen,
                            nm_cut_coul_long, "k_nm_cut_coul_long");
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

  nm1.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack4(ntypes, lj_types, nm1, host_write, host_e0nm, host_nm, host_r0n, host_r0m);

  nm2.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack4(ntypes, lj_types, nm2, host_write, host_nn, host_mm, host_cut_ljsq,
                         host_offset);

  cutsq.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack1(ntypes, lj_types, cutsq, host_write, host_cutsq);

  sp_lj.alloc(8,*(this->ucl_device),UCL_READ_ONLY);
  for (int i=0; i<4; i++) {
    host_write[i]=host_special_lj[i];
    host_write[i+4]=host_special_coul[i];
  }
  ucl_copy(sp_lj,host_write,8,false);

  _cut_coulsq=host_cut_coulsq;
  _qqrd2e=host_qqrd2e;
  _g_ewald=host_g_ewald;

  _allocated=true;
  this->_max_bytes=nm1.row_bytes()+nm2.row_bytes()+cutsq.row_bytes()+sp_lj.row_bytes();
  return 0;
}

template <class numtyp, class acctyp>
void NMCutCoulLongT::reinit(const int ntypes, double **host_e0nm, double **host_nm,
                            double **host_r0n, double **host_r0m, double **host_nn,
                            double **host_mm, double **host_cut_ljsq, double **host_offset,
                            double **host_cutsq) {
  // Allocate a host write buffer for data initialization
  UCL_H_Vec<numtyp> host_write(_lj_types*_lj_types*32,*(this->ucl_device), UCL_WRITE_ONLY);

  for (int i=0; i<_lj_types*_lj_types; i++)
    host_write[i]=0.0;

  this->atom->type_pack4(ntypes, _lj_types, nm1, host_write, host_e0nm, host_nm, host_r0n,
                         host_r0m);
  this->atom->type_pack4(ntypes, _lj_types, nm2, host_write, host_nn, host_mm, host_cut_ljsq,
                         host_offset);
  this->atom->type_pack1(ntypes, _lj_types, cutsq, host_write, host_cutsq);
}

template <class numtyp, class acctyp>
void NMCutCoulLongT::clear() {
  if (!_allocated)
    return;
  _allocated=false;

  nm1.clear();
  nm2.clear();
  cutsq.clear();
  sp_lj.clear();
  this->clear_atomic();
}

template <class numtyp, class acctyp>
double NMCutCoulLongT::host_memory_usage() const {
  return this->host_memory_usage_atomic()+sizeof(NMCutCoulLong<numtyp,acctyp>);
}

// ---------------------------------------------------------------------------
// Calculate energies, forces, and torques
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int NMCutCoulLongT::loop(const int eflag, const int vflag) {
  // Compute the block size and grid size to keep all cores busy
  const int BX=this->block_size();
  int GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/
                               (BX/this->_threads_per_atom)));

  int ainum=this->ans->inum();
  int nbor_pitch=this->nbor->nbor_pitch();
  this->time_pair.start();
  if (shared_types) {
    this->k_pair_sel->set_size(GX,BX);
    this->k_pair_sel->run(&this->atom->x, &nm1, &nm2, &cutsq, &sp_lj, &this->nbor->dev_nbor,
                          &this->_nbor_data->begin(), &this->ans->force, &this->ans->engv, &eflag,
                          &vflag, &ainum, &nbor_pitch, &this->atom->q, &_cut_coulsq, &_qqrd2e,
                          &_g_ewald, &this->_threads_per_atom);
  } else {
    this->k_pair.set_size(GX,BX);
    this->k_pair.run(&this->atom->x, &nm1, &nm2, &cutsq, &_lj_types, &sp_lj, &this->nbor->dev_nbor,
                     &this->_nbor_data->begin(), &this->ans->force, &this->ans->engv, &eflag,
                     &vflag, &ainum, &nbor_pitch, &this->atom->q, &_cut_coulsq, &_qqrd2e, &_g_ewald,
                     &this->_threads_per_atom);
  }
  this->time_pair.stop();
  return GX;
}

template class NMCutCoulLong<PRECISION,ACC_PRECISION>;
}
