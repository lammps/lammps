/***************************************************************************
                                 lj_spica.cpp
                             -------------------
                            W. Michael Brown (ORNL)

  Class for acceleration of the lj/spica/cut pair style

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : brownw@ornl.gov
 ***************************************************************************/

#if defined(USE_OPENCL)
#include "lj_spica_cl.h"
#elif defined(USE_CUDART)
const char *lj_spica=0;
#else
#include "lj_spica_cubin.h"
#endif

#include "lal_lj_spica.h"
#include <cassert>
namespace LAMMPS_AL {
#define CGCMMT CGCMM<numtyp, acctyp>

extern Device<PRECISION,ACC_PRECISION> device;

template <class numtyp, class acctyp>
CGCMMT::CGCMM() : BaseAtomic<numtyp,acctyp>(), _allocated(false) {
}

template <class numtyp, class acctyp>
CGCMMT::~CGCMM() {
  clear();
}

template <class numtyp, class acctyp>
int CGCMMT::bytes_per_atom(const int max_nbors) const {
  return this->bytes_per_atom_atomic(max_nbors);
}

template <class numtyp, class acctyp>
int CGCMMT::init(const int ntypes, double **host_cutsq,
                          int **host_cg_type, double **host_lj1,
                          double **host_lj2, double **host_lj3,
                          double **host_lj4, double **host_offset,
                          double *host_special_lj, const int nlocal,
                          const int nall, const int max_nbors,
                          const int maxspecial, const double cell_size,
                          const double gpu_split, FILE *_screen) {
  int success;
  success=this->init_atomic(nlocal,nall,max_nbors,maxspecial,cell_size,gpu_split,
                            _screen,lj_spica,"k_lj_spica");
  if (success!=0)
    return success;

  // If atom type constants fit in shared memory use fast kernel
  int spica_types=ntypes;
  shared_types=false;
  int max_shared_types=this->device->max_shared_types();
  if (spica_types<=max_shared_types && this->_block_size>=max_shared_types) {
    spica_types=max_shared_types;
    shared_types=true;
  }
  _spica_types=spica_types;

  // Allocate a host write buffer for data initialization
  UCL_H_Vec<numtyp> host_write(spica_types*spica_types*32,*(this->ucl_device),
                               UCL_WRITE_ONLY);

  for (int i=0; i<spica_types*spica_types; i++)
    host_write[i]=0.0;

  lj1.alloc(spica_types*spica_types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack4(ntypes,spica_types,lj1,host_write,host_cutsq,
                         host_cg_type,host_lj1,host_lj2);

  lj3.alloc(spica_types*spica_types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack4(ntypes,spica_types,lj3,host_write,host_lj3,host_lj4,
                         host_offset);

  UCL_H_Vec<double> dview;
  sp_lj.alloc(4,*(this->ucl_device),UCL_READ_ONLY);
  dview.view(host_special_lj,4,*(this->ucl_device));
  ucl_copy(sp_lj,dview,false);

  _allocated=true;
  this->_max_bytes=lj1.row_bytes()+lj3.row_bytes()+sp_lj.row_bytes();
  return 0;
}

template <class numtyp, class acctyp>
void CGCMMT::clear() {
  if (!_allocated)
    return;
  _allocated=false;

  lj1.clear();
  lj3.clear();
  sp_lj.clear();
  this->clear_atomic();
}

template <class numtyp, class acctyp>
double CGCMMT::host_memory_usage() const {
  return this->host_memory_usage_atomic()+sizeof(CGCMM<numtyp,acctyp>);
}

// ---------------------------------------------------------------------------
// Calculate energies, forces, and torques
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int CGCMMT::loop(const int eflag, const int vflag) {
  // Compute the block size and grid size to keep all cores busy
  const int BX=this->block_size();
  int GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/
                               (BX/this->_threads_per_atom)));

  int ainum=this->ans->inum();
  int nbor_pitch=this->nbor->nbor_pitch();
  this->time_pair.start();
  if (shared_types) {
    this->k_pair_sel->set_size(GX,BX);
    this->k_pair_sel->run(&this->atom->x, &lj1, &lj3, &sp_lj,
                          &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                          &this->ans->force, &this->ans->engv, &eflag,
                          &vflag, &ainum, &nbor_pitch,
                          &this->_threads_per_atom);
  } else {
    this->k_pair.set_size(GX,BX);
    this->k_pair.run(&this->atom->x, &lj1, &lj3,
                     &_spica_types, &sp_lj, &this->nbor->dev_nbor,
                     &this->_nbor_data->begin(), &this->ans->force,
                     &this->ans->engv, &eflag, &vflag, &ainum,
                     &nbor_pitch, &this->_threads_per_atom);
  }
  this->time_pair.stop();
  return GX;
}

template class CGCMM<PRECISION,ACC_PRECISION>;
}
