/***************************************************************************
                           coul_slater_long_ext.cpp
                           ------------------------
                           Trung Nguyen (U Chicago)

  Class for acceleration of the coul/slater/long pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : September 2023
    email                : ndactrung@gmail.com
 ***************************************************************************/

#if defined(USE_OPENCL)
#include "coul_slater_long_cl.h"
#elif defined(USE_CUDART)
const char *coul_slater_long=0;
#else
#include "coul_slater_long_cubin.h"
#endif

#include "lal_coul_slater_long.h"
#include <cassert>
namespace LAMMPS_AL {
#define CoulSlaterLongT CoulSlaterLong<numtyp, acctyp>

extern Device<PRECISION,ACC_PRECISION> pair_gpu_device;

template <class numtyp, class acctyp>
CoulSlaterLongT::CoulSlaterLong() : BaseCharge<numtyp,acctyp>(), _allocated(false) {
}

template <class numtyp, class acctyp>
CoulSlaterLongT::~CoulSlaterLong() {
  clear();
}

template <class numtyp, class acctyp>
int CoulSlaterLongT::bytes_per_atom(const int max_nbors) const {
  return this->bytes_per_atom_atomic(max_nbors);
}

template <class numtyp, class acctyp>
int CoulSlaterLongT::init(const int ntypes, double **host_scale,
                    const int nlocal, const int nall, const int max_nbors,
                    const int maxspecial, const double cell_size,
                    const double gpu_split, FILE *_screen,
                    const double host_cut_coulsq, double *host_special_coul,
                    const double qqrd2e, const double g_ewald, double lamda) {
  int success;
  success=this->init_atomic(nlocal,nall,max_nbors,maxspecial,cell_size,
                            gpu_split,_screen,coul_slater_long,"k_coul_slater_long");
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

  scale.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack1(ntypes,lj_types,scale,host_write,host_scale);

  sp_cl.alloc(4,*(this->ucl_device),UCL_READ_ONLY);
  for (int i=0; i<4; i++) {
    host_write[i]=host_special_coul[i];
  }
  ucl_copy(sp_cl,host_write,4,false);

  _cut_coulsq=host_cut_coulsq;
  _qqrd2e=qqrd2e;
  _g_ewald=g_ewald;
  _lamda=lamda;

  _allocated=true;
  this->_max_bytes=scale.row_bytes()+sp_cl.row_bytes();
  return 0;
}

template <class numtyp, class acctyp>
void CoulSlaterLongT::reinit(const int ntypes, double **host_scale) {
  UCL_H_Vec<numtyp> hscale(_lj_types*_lj_types,*(this->ucl_device),
                           UCL_WRITE_ONLY);
  this->atom->type_pack1(ntypes,_lj_types,scale,hscale,host_scale);
}

template <class numtyp, class acctyp>
void CoulSlaterLongT::clear() {
  if (!_allocated)
    return;
  _allocated=false;

  scale.clear();
  sp_cl.clear();
  this->clear_atomic();
}

template <class numtyp, class acctyp>
double CoulSlaterLongT::host_memory_usage() const {
  return this->host_memory_usage_atomic()+sizeof(CoulSlaterLong<numtyp,acctyp>);
}

// ---------------------------------------------------------------------------
// Calculate energies, forces, and torques
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int CoulSlaterLongT::loop(const int eflag, const int vflag) {
  // Compute the block size and grid size to keep all cores busy
  const int BX=this->block_size();
  int GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/
                               (BX/this->_threads_per_atom)));

  int ainum=this->ans->inum();
  int nbor_pitch=this->nbor->nbor_pitch();
  this->time_pair.start();
  if (shared_types) {
    this->k_pair_sel->set_size(GX,BX);
    this->k_pair_sel->run(&this->atom->x, &scale, &sp_cl,
                          &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                          &this->ans->force, &this->ans->engv,
                          &eflag, &vflag, &ainum, &nbor_pitch,
                          &this->atom->q, &_cut_coulsq, &_qqrd2e, &_g_ewald,
                          &_lamda, &this->_threads_per_atom);
  } else {
    this->k_pair.set_size(GX,BX);
    this->k_pair.run(&this->atom->x, &scale, &_lj_types, &sp_cl,
                     &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                     &this->ans->force, &this->ans->engv, &eflag, &vflag,
                     &ainum, &nbor_pitch, &this->atom->q, &_cut_coulsq,
                     &_qqrd2e, &_g_ewald, &_lamda, &this->_threads_per_atom);
  }
  this->time_pair.stop();
  return GX;
}

template class CoulSlaterLong<PRECISION,ACC_PRECISION>;
}
