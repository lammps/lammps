/***************************************************************************
                                 coul_wolf_cs.cpp
                             -------------------

  Class for acceleration of the coul/wolf/cs pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : October 2026
 ***************************************************************************/

#if defined(USE_OPENCL)
#include "coul_wolf_cs_cl.h"
#elif defined(USE_CUDART)
const char *coul_wolf_cs=0;
#else
#include "coul_wolf_cs_cubin.h"
#endif

#include "lal_coul_wolf_cs.h"
#include <cassert>
namespace LAMMPS_AL {
#define CoulWolfCST CoulWolfCS<numtyp, acctyp>

extern Device<PRECISION,ACC_PRECISION> device;

template <class numtyp, class acctyp>
CoulWolfCST::CoulWolfCS() : BaseCharge<numtyp,acctyp>(),
                                    _allocated(false) {
}

template <class numtyp, class acctyp>
CoulWolfCST::~CoulWolfCS() {
  clear();
}

template <class numtyp, class acctyp>
int CoulWolfCST::bytes_per_atom(const int max_nbors) const {
  return this->bytes_per_atom_atomic(max_nbors);
}

template <class numtyp, class acctyp>
int CoulWolfCST::init(const int ntypes, const int nlocal, const int nall,
                   const int max_nbors, const int maxspecial,
                   const double cell_size, FILE *_screen,
                   const double host_cut_coulsq, double *host_special_coul,
                   const double qqrd2e, const double e_shift, const double f_shift,
                   const double alpha) {
  int success;
  success=this->init_atomic(nlocal,nall,max_nbors,maxspecial,cell_size,
                            _screen,coul_wolf_cs,"k_coul_wolf_cs");
  if (success!=0)
    return success;

  _cut_coulsq=host_cut_coulsq;
  _e_shift=e_shift;
  _f_shift=f_shift;
  _alpha=alpha;

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

  sp_lj.alloc(4,*(this->ucl_device),UCL_READ_ONLY);
  for (int i=0; i<4; i++) {
    host_write[i]=host_special_coul[i];
  }
  ucl_copy(sp_lj,host_write,4,false);

  _qqrd2e=qqrd2e;

  _allocated=true;
  this->_max_bytes=sp_lj.row_bytes();
  return 0;
}

template <class numtyp, class acctyp>
void CoulWolfCST::clear() {
  if (!_allocated)
    return;
  _allocated=false;

  sp_lj.clear();
  this->clear_atomic();
}

template <class numtyp, class acctyp>
double CoulWolfCST::host_memory_usage() const {
  return this->host_memory_usage_atomic()+sizeof(CoulWolfCS<numtyp,acctyp>);
}

// ---------------------------------------------------------------------------
// Calculate energies, forces, and torques
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int CoulWolfCST::loop(const int eflag, const int vflag) {
  // Compute the block size and grid size to keep all cores busy
  const int BX=this->block_size();
  int GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/
                               (BX/this->_threads_per_atom)));

  int ainum=this->ans->inum();
  int nbor_pitch=this->nbor->nbor_pitch();
  this->time_pair.start();
  if (shared_types) {
    this->k_pair_sel->set_size(GX,BX);
    this->k_pair_sel->run(&this->atom->x, &sp_lj,
                          &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                          &this->ans->force, &this->ans->engv, &eflag,
                          &vflag, &ainum, &nbor_pitch, &this->atom->q,
                          &_cut_coulsq, &_qqrd2e, &_e_shift, &_f_shift, &_alpha,
                          &this->_threads_per_atom);
  } else {
    this->k_pair.set_size(GX,BX);
    this->k_pair.run(&this->atom->x, &_lj_types, &sp_lj,
                     &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                     &this->ans->force, &this->ans->engv,
                     &eflag, &vflag, &ainum, &nbor_pitch, &this->atom->q,
                     &_cut_coulsq, &_qqrd2e, &_e_shift, &_f_shift, &_alpha,
                     &this->_threads_per_atom);
  }
  this->time_pair.stop();
  return GX;
}

template class CoulWolfCS<PRECISION,ACC_PRECISION>;
}
