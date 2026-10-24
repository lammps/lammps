/***************************************************************************
                                 coul_diel.cpp
                             -------------------

  Class for acceleration of the coul/diel pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : October 2026
 ***************************************************************************/

#if defined(USE_OPENCL)
#include "coul_diel_cl.h"
#elif defined(USE_CUDART)
const char *coul_diel=0;
#else
#include "coul_diel_cubin.h"
#endif

#include "lal_coul_diel.h"
#include <cassert>
namespace LAMMPS_AL {
#define CoulDielT CoulDiel<numtyp, acctyp>

extern Device<PRECISION,ACC_PRECISION> device;

template <class numtyp, class acctyp>
CoulDielT::CoulDiel() : BaseCharge<numtyp,acctyp>(), _allocated(false) {
}

template <class numtyp, class acctyp>
CoulDielT::~CoulDiel() {
  clear();
}

template <class numtyp, class acctyp>
int CoulDielT::bytes_per_atom(const int max_nbors) const {
  return this->bytes_per_atom_atomic(max_nbors);
}

template <class numtyp, class acctyp>
int CoulDielT::init(const int ntypes, double **host_cutsq, double **host_rme, double **host_sigmae,
                    double **host_offset, double *host_special_lj, double *host_special_coul,
                    const int nlocal, const int nall, const int max_nbors, const int maxspecial,
                    const double cell_size, FILE *_screen, const double host_qqrd2e,
                    const double host_eps_s, const double host_a_eps, const double host_b_eps) {
  int success;
  success=this->init_atomic(nlocal, nall, max_nbors, maxspecial, cell_size, _screen, coul_diel,
                            "k_coul_diel");
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

  cutsq.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack4(ntypes, lj_types, cutsq, host_write, host_cutsq, host_rme, host_sigmae,
                         host_offset);

  sp_lj.alloc(8,*(this->ucl_device),UCL_READ_ONLY);
  for (int i=0; i<4; i++) {
    host_write[i]=host_special_lj[i];
    host_write[i+4]=host_special_coul[i];
  }
  ucl_copy(sp_lj,host_write,8,false);

  _qqrd2e=host_qqrd2e;
  _eps_s=host_eps_s;
  _a_eps=host_a_eps;
  _b_eps=host_b_eps;

  _allocated=true;
  this->_max_bytes=cutsq.row_bytes()+sp_lj.row_bytes();
  return 0;
}

template <class numtyp, class acctyp>
void CoulDielT::clear() {
  if (!_allocated)
    return;
  _allocated=false;

  cutsq.clear();
  sp_lj.clear();
  this->clear_atomic();
}

template <class numtyp, class acctyp>
double CoulDielT::host_memory_usage() const {
  return this->host_memory_usage_atomic()+sizeof(CoulDiel<numtyp,acctyp>);
}

// ---------------------------------------------------------------------------
// Calculate energies, forces, and torques
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int CoulDielT::loop(const int eflag, const int vflag) {
  // Compute the block size and grid size to keep all cores busy
  const int BX=this->block_size();
  int GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/
                               (BX/this->_threads_per_atom)));

  int ainum=this->ans->inum();
  int nbor_pitch=this->nbor->nbor_pitch();
  this->time_pair.start();
  if (shared_types) {
    this->k_pair_sel->set_size(GX,BX);
    this->k_pair_sel->run(&this->atom->x, &cutsq, &sp_lj, &this->nbor->dev_nbor,
                          &this->_nbor_data->begin(), &this->ans->force, &this->ans->engv, &eflag,
                          &vflag, &ainum, &nbor_pitch, &this->atom->q, &_qqrd2e, &_eps_s, &_a_eps,
                          &_b_eps, &this->_threads_per_atom);
  } else {
    this->k_pair.set_size(GX,BX);
    this->k_pair.run(&this->atom->x, &cutsq, &_lj_types, &sp_lj, &this->nbor->dev_nbor,
                     &this->_nbor_data->begin(), &this->ans->force, &this->ans->engv, &eflag,
                     &vflag, &ainum, &nbor_pitch, &this->atom->q, &_qqrd2e, &_eps_s, &_a_eps,
                     &_b_eps, &this->_threads_per_atom);
  }
  this->time_pair.stop();
  return GX;
}

template class CoulDiel<PRECISION,ACC_PRECISION>;
}
