/***************************************************************************
                                   colloid.cpp
                             -------------------
                            Trung Dac Nguyen (ORNL)

  Class for acceleration of the colloid pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : nguyentd@ornl.gov
 ***************************************************************************/

#ifdef USE_OPENCL
#include "colloid_cl.h"
#elif defined(USE_CUDART)
const char *colloid=0;
#else
#include "colloid_cubin.h"
#endif

#include "lal_colloid.h"
#include <cassert>
namespace LAMMPS_AL {
#define ColloidT Colloid<numtyp, acctyp>

extern Device<PRECISION,ACC_PRECISION> device;

template <class numtyp, class acctyp>
ColloidT::Colloid() : BaseAtomic<numtyp,acctyp>(), _allocated(false) {
}

template <class numtyp, class acctyp>
ColloidT::~Colloid() {
  clear();
}

template <class numtyp, class acctyp>
int ColloidT::bytes_per_atom(const int max_nbors) const {
  return this->bytes_per_atom_atomic(max_nbors);
}

template <class numtyp, class acctyp>
int ColloidT::init(const int ntypes,
                   double **host_cutsq, double **host_lj1,
                   double **host_lj2, double **host_lj3,
                   double **host_lj4, double **host_offset,
                   double *host_special_lj, double **host_a12,
                   double **host_a1, double **host_a2,
                   double **host_d1, double **host_d2,
                   double **host_sigma3, double **host_sigma6,
                   int **host_form, const int nlocal,
                   const int nall, const int max_nbors,
                   const int maxspecial, const double cell_size,
                   const double gpu_split, FILE *_screen) {
  int success;
  success=this->init_atomic(nlocal,nall,max_nbors,maxspecial,cell_size,gpu_split,
                            _screen,colloid,"k_colloid");
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
                         host_cutsq);

  lj3.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack4(ntypes,lj_types,lj3,host_write,host_lj3,host_lj4,
                         host_offset);

  colloid1.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack4(ntypes,lj_types,colloid1,host_write,host_a12,host_a1,
                         host_a2);

  colloid2.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack4(ntypes,lj_types,colloid2,host_write,host_d1,host_d2,
                         host_sigma3,host_sigma6);

  UCL_H_Vec<int> dview_form(lj_types*lj_types,*(this->ucl_device),
                             UCL_WRITE_ONLY);
  for (int i=0; i<lj_types*lj_types; i++) dview_form[i]=0;

  form.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  for (int i=0; i<ntypes; i++)
    for (int j=0; j<ntypes; j++) {
      dview_form[i*lj_types+j]=host_form[i][j];
  }
  ucl_copy(form,dview_form,false);

  UCL_H_Vec<double> dview;
  sp_lj.alloc(4,*(this->ucl_device),UCL_READ_ONLY);
  dview.view(host_special_lj,4,*(this->ucl_device));
  ucl_copy(sp_lj,dview,false);

  _allocated=true;
  this->_max_bytes=lj1.row_bytes()+lj3.row_bytes()
        +colloid1.row_bytes()+colloid2.row_bytes()+sp_lj.row_bytes();
  return 0;
}

template <class numtyp, class acctyp>
void ColloidT::clear() {
  if (!_allocated)
    return;
  _allocated=false;

  lj1.clear();
  lj3.clear();
  colloid1.clear();
  colloid2.clear();
  form.clear();
  sp_lj.clear();
  this->clear_atomic();
}

template <class numtyp, class acctyp>
double ColloidT::host_memory_usage() const {
  return this->host_memory_usage_atomic()+sizeof(Colloid<numtyp,acctyp>);
}

// ---------------------------------------------------------------------------
// Calculate energies, forces, and torques
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void ColloidT::loop(const bool _eflag, const bool _vflag) {
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
    this->k_pair_fast.run(&this->atom->x, &lj1, &lj3, &sp_lj,
                          &colloid1, &colloid2, &form,
                          &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                          &this->ans->force, &this->ans->engv, &eflag, &vflag,
                          &ainum, &nbor_pitch, &this->_threads_per_atom);
  } else {
    this->k_pair.set_size(GX,BX);
    this->k_pair.run(&this->atom->x, &lj1, &lj3, &_lj_types, &sp_lj,
                     &colloid1, &colloid2, &form,
                     &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                     &this->ans->force, &this->ans->engv, &eflag, &vflag,
                     &ainum, &nbor_pitch, &this->_threads_per_atom);
  }
  this->time_pair.stop();
}

template class Colloid<PRECISION,ACC_PRECISION>;
}
