/***************************************************************************
                                   dpd.cpp
                             -------------------
                            Trung Dac Nguyen (ORNL)

  Class for acceleration of the dpd pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : Jan 15, 2014
    email                : nguyentd@ornl.gov
 ***************************************************************************/

#if defined(USE_OPENCL)
#include "dpd_cl.h"
#elif defined(USE_CUDART)
const char *dpd=0;
#else
#include "dpd_cubin.h"
#endif

#include "lal_dpd.h"
#include <cassert>
namespace LAMMPS_AL {
#define DPDT DPD<numtyp, acctyp>

extern Device<PRECISION,ACC_PRECISION> device;

template <class numtyp, class acctyp>
DPDT::DPD() : BaseDPD<numtyp,acctyp>(), _allocated(false) {
}

template <class numtyp, class acctyp>
DPDT::~DPD() {
  clear();
}

template <class numtyp, class acctyp>
int DPDT::bytes_per_atom(const int max_nbors) const {
  return this->bytes_per_atom_atomic(max_nbors);
}

template <class numtyp, class acctyp>
int DPDT::init(const int ntypes,
               double **host_cutsq, double **host_a0,
               double **host_gamma, double **host_sigma,
               double **host_cut, double *host_special_lj,
               const bool tstat_only,
               const int nlocal, const int nall,
               const int max_nbors, const int maxspecial,
               const double cell_size,
               const double gpu_split, FILE *_screen) {
  const int max_shared_types=this->device->max_shared_types();

  int onetype=0;
  #ifdef USE_OPENCL
  if (maxspecial==0)
    for (int i=1; i<ntypes; i++)
      for (int j=i; j<ntypes; j++)
        if (host_cutsq[i][j]>0) {
          if (onetype>0)
            onetype=-1;
          else if (onetype==0)
            onetype=i*max_shared_types+j;
        }
  if (onetype<0) onetype=0;
  #endif

  int success;
  success=this->init_atomic(nlocal,nall,max_nbors,maxspecial,cell_size,
                            gpu_split,_screen,dpd,"k_dpd",onetype);
  if (success!=0)
    return success;

  // If atom type constants fit in shared memory use fast kernel
  int lj_types=ntypes;
  shared_types=false;
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

  coeff.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack4(ntypes,lj_types,coeff,host_write,host_a0,host_gamma,
                         host_sigma,host_cut);

  UCL_H_Vec<numtyp> host_rsq(lj_types*lj_types,*(this->ucl_device),
                             UCL_WRITE_ONLY);
  cutsq.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack1(ntypes,lj_types,cutsq,host_rsq,host_cutsq);

  double special_sqrt[4];
  special_sqrt[0] = sqrt(host_special_lj[0]);
  special_sqrt[1] = sqrt(host_special_lj[1]);
  special_sqrt[2] = sqrt(host_special_lj[2]);
  special_sqrt[3] = sqrt(host_special_lj[3]);

  UCL_H_Vec<double> dview;
  sp_lj.alloc(4,*(this->ucl_device),UCL_READ_ONLY);
  dview.view(host_special_lj,4,*(this->ucl_device));
  ucl_copy(sp_lj,dview,false);
  sp_sqrt.alloc(4,*(this->ucl_device),UCL_READ_ONLY);
  dview.view(special_sqrt,4,*(this->ucl_device));
  ucl_copy(sp_sqrt,dview,false);

  _tstat_only = 0;
  if (tstat_only) _tstat_only=1;

  _allocated=true;
  this->_max_bytes=coeff.row_bytes()+cutsq.row_bytes()+sp_lj.row_bytes()+sp_sqrt.row_bytes();
  return 0;
}

template <class numtyp, class acctyp>
void DPDT::clear() {
  if (!_allocated)
    return;
  _allocated=false;

  coeff.clear();
  cutsq.clear();
  sp_lj.clear();
  sp_sqrt.clear();
  this->clear_atomic();
}

template <class numtyp, class acctyp>
double DPDT::host_memory_usage() const {
  return this->host_memory_usage_atomic()+sizeof(DPD<numtyp,acctyp>);
}

// ---------------------------------------------------------------------------
// Calculate energies, forces, and torques
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int DPDT::loop(const int eflag, const int vflag) {
  // Compute the block size and grid size to keep all cores busy
  const int BX=this->block_size();
  int GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/
                               (BX/this->_threads_per_atom)));

  int ainum=this->ans->inum();
  int nbor_pitch=this->nbor->nbor_pitch();
  this->time_pair.start();
  if (shared_types) {
    this->k_pair_sel->set_size(GX,BX);
    this->k_pair_sel->run(&this->atom->x, &coeff, &sp_lj, &sp_sqrt,
                          &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                          &this->ans->force, &this->ans->engv, &eflag,
                          &vflag, &ainum, &nbor_pitch, &this->atom->v, &cutsq,
                          &this->_dtinvsqrt, &this->_seed, &this->_timestep,
                          &this->_tstat_only, &this->_threads_per_atom);
  } else {
    this->k_pair.set_size(GX,BX);
    this->k_pair.run(&this->atom->x, &coeff, &_lj_types, &sp_lj, &sp_sqrt,
                     &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                     &this->ans->force, &this->ans->engv, &eflag, &vflag,
                     &ainum, &nbor_pitch, &this->atom->v, &cutsq, &this->_dtinvsqrt,
                     &this->_seed, &this->_timestep, &this->_tstat_only,
                     &this->_threads_per_atom);
  }
  this->time_pair.stop();
  return GX;
}

template <class numtyp, class acctyp>
void DPDT::update_coeff(int ntypes, double **host_a0, double **host_gamma,
                        double **host_sigma, double **host_cut)
{
  UCL_H_Vec<numtyp> host_write(_lj_types*_lj_types*32,*(this->ucl_device),
                               UCL_WRITE_ONLY);
  this->atom->type_pack4(ntypes,_lj_types,coeff,host_write,host_a0,host_gamma,
                         host_sigma,host_cut);
}

template class DPD<PRECISION,ACC_PRECISION>;
}
