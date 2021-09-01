/***************************************************************************
                                 amoeba.cpp
                             -------------------
                          Trung Dac Nguyen (Northwestern)

  Class for acceleration of the amoeba pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : trung.nguyen@northwestern.edu
 ***************************************************************************/

#if defined(USE_OPENCL)
#include "amoeba_cl.h"
#elif defined(USE_CUDART)
const char *amoeba=0;
#else
#include "amoeba_cubin.h"
#endif

#include "lal_amoeba.h"
#include <cassert>
namespace LAMMPS_AL {
#define AmoebaT Amoeba<numtyp, acctyp>

extern Device<PRECISION,ACC_PRECISION> device;

template <class numtyp, class acctyp>
AmoebaT::Amoeba() : BaseAmoeba<numtyp,acctyp>(),
  _allocated(false) {
}

template <class numtyp, class acctyp>
AmoebaT::~Amoeba() {
  clear();
}

template <class numtyp, class acctyp>
int AmoebaT::bytes_per_atom(const int max_nbors) const {
  return this->bytes_per_atom_atomic(max_nbors);
}

template <class numtyp, class acctyp>
int AmoebaT::init(const int ntypes, const int max_amtype, const double *host_pdamp,
                  const double *host_thole, const double *host_dirdamp, 
                  const double *host_special_polar_wscale,
                  const double *host_special_polar_piscale,
                  const double *host_special_polar_pscale,
                  const int nlocal, const int nall, const int max_nbors,
                  const int maxspecial, const int maxspecial15,
                  const double cell_size, const double gpu_split, FILE *_screen,
                  const double aewald, const double felec,
                  const double off2, const double polar_dscale,
                  const double polar_uscale) {
  int success;
  success=this->init_atomic(nlocal,nall,max_nbors,maxspecial,maxspecial15,
                            cell_size,gpu_split,_screen,amoeba,"k_amoeba_polar");
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

  UCL_H_Vec<numtyp4> host_write(max_amtype, *(this->ucl_device), UCL_WRITE_ONLY);
  for (int i = 0; i < max_amtype; i++) {
    host_write[i].x = host_pdamp[i];
    host_write[i].y = host_thole[i];
    host_write[i].z = host_dirdamp[i];
    host_write[i].w = (numtyp)0;
  }

  damping.alloc(max_amtype,*(this->ucl_device), UCL_READ_ONLY);
  ucl_copy(damping,host_write,false);

  UCL_H_Vec<numtyp4> dview(5, *(this->ucl_device), UCL_WRITE_ONLY);
  sp_polar.alloc(5,*(this->ucl_device),UCL_READ_ONLY);
  for (int i=0; i<5; i++) {
    dview[i].x=host_special_polar_wscale[i];
    dview[i].y=host_special_polar_piscale[i];
    dview[i].z=host_special_polar_pscale[i];
    dview[i].w=(numtyp)0;
  }
  ucl_copy(sp_polar,dview,5,false);

  _aewald = aewald;
  _felec = felec;
  _off2 = off2;
  _polar_dscale = polar_dscale;
  _polar_uscale = polar_uscale;

  _allocated=true;
  this->_max_bytes=damping.row_bytes()
    + sp_polar.row_bytes()
    + this->_tep.row_bytes();
  return 0;
}

template <class numtyp, class acctyp>
void AmoebaT::clear() {
  if (!_allocated)
    return;
  _allocated=false;

  damping.clear();
  sp_polar.clear();
  
  this->clear_atomic();
}

template <class numtyp, class acctyp>
double AmoebaT::host_memory_usage() const {
  return this->host_memory_usage_atomic()+sizeof(Amoeba<numtyp,acctyp>);
}

// ---------------------------------------------------------------------------
// Calculate energies, forces, and torques
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int AmoebaT::loop(const int eflag, const int vflag) {
  // Compute the block size and grid size to keep all cores busy
  const int BX=this->block_size();
  int GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/
                               (BX/this->_threads_per_atom)));

  int _nall=this->atom->nall();
  int ainum=this->ans->inum();
  int nbor_pitch=this->nbor->nbor_pitch();
  this->time_pair.start();

  this->k_polar.set_size(GX,BX);

  this->k_polar.run(&this->atom->x, &this->atom->extra,
                    &damping, &sp_polar,
                    &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                    &this->ans->force, &this->ans->engv, &this->_tep,
                    &eflag, &vflag, &ainum, &_nall, &nbor_pitch,
                    &this->_threads_per_atom,
                    &_aewald, &_felec, &_off2, &_polar_dscale, &_polar_uscale);
  this->time_pair.stop();
  return GX;
}

template class Amoeba<PRECISION,ACC_PRECISION>;
}
