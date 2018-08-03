/***************************************************************************
                                   sw.cpp
                             -------------------
                            W. Michael Brown (ORNL)

  Class for acceleration of the sw pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : Tue March 26, 2013
    email                : brownw@ornl.gov
 ***************************************************************************/

#if defined(USE_OPENCL)
#include "sw_cl.h"
#elif defined(USE_CUDART)
const char *sw=0;
#else
#include "sw_cubin.h"
#endif

#include "lal_sw.h"
#include <cassert>
using namespace LAMMPS_AL;
#define SWT SW<numtyp, acctyp>

extern Device<PRECISION,ACC_PRECISION> device;

template <class numtyp, class acctyp>
SWT::SW() : BaseThree<numtyp,acctyp>(), _allocated(false) {
}

template <class numtyp, class acctyp>
SWT::~SW() {
  clear();
}

template <class numtyp, class acctyp>
int SWT::bytes_per_atom(const int max_nbors) const {
  return this->bytes_per_atom_atomic(max_nbors);
}

template <class numtyp, class acctyp>
int SWT::init(const int ntypes, const int nlocal, const int nall, const int max_nbors,
           const double cell_size, const double gpu_split, FILE *_screen,
           int* host_map, const int nelements, int*** host_elem2param, const int nparams,
           const double* epsilon, const double* sigma,
           const double* lambda, const double* gamma,
           const double* costheta, const double* biga,
           const double* bigb, const double* powerp,
           const double* powerq, const double* cut, const double* cutsq)
{
  int success;
  success=this->init_three(nlocal,nall,max_nbors,0,cell_size,gpu_split,
                           _screen,sw,"k_sw","k_sw_three_center",
                           "k_sw_three_end","k_sw_short_nbor");
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

  _nparams = nparams;
  _nelements = nelements;

  UCL_H_Vec<numtyp4> dview(nparams,*(this->ucl_device),
                             UCL_WRITE_ONLY);

  for (int i=0; i<nparams; i++) {
    dview[i].x=(numtyp)0;
    dview[i].y=(numtyp)0;
    dview[i].z=(numtyp)0;
    dview[i].w=(numtyp)0;
  }

  // pack coefficients into arrays
  sw1.alloc(nparams,*(this->ucl_device),UCL_READ_ONLY);

  for (int i=0; i<nparams; i++) {
    dview[i].x=static_cast<numtyp>(epsilon[i]);
    dview[i].y=static_cast<numtyp>(sigma[i]);
    dview[i].z=static_cast<numtyp>(lambda[i]);
    dview[i].w=static_cast<numtyp>(gamma[i]);
  }

  ucl_copy(sw1,dview,false);
  sw1_tex.get_texture(*(this->pair_program),"sw1_tex");
  sw1_tex.bind_float(sw1,4);

  sw2.alloc(nparams,*(this->ucl_device),UCL_READ_ONLY);

  for (int i=0; i<nparams; i++) {
    dview[i].x=static_cast<numtyp>(biga[i]);
    dview[i].y=static_cast<numtyp>(bigb[i]);
    dview[i].z=static_cast<numtyp>(powerp[i]);
    dview[i].w=static_cast<numtyp>(powerq[i]);
  }

  ucl_copy(sw2,dview,false);
  sw2_tex.get_texture(*(this->pair_program),"sw2_tex");
  sw2_tex.bind_float(sw2,4);

  sw3.alloc(nparams,*(this->ucl_device),UCL_READ_ONLY);

  for (int i=0; i<nparams; i++) {
    double sw_cut = cut[i];
    double sw_cutsq = cutsq[i];
    if (sw_cutsq>=sw_cut*sw_cut)
      sw_cutsq=sw_cut*sw_cut-1e-4;
    dview[i].x=static_cast<numtyp>(sw_cut);
    dview[i].y=static_cast<numtyp>(sw_cutsq);
    dview[i].z=static_cast<numtyp>(costheta[i]);
    dview[i].w=(numtyp)0;
  }

  ucl_copy(sw3,dview,false);
  sw3_tex.get_texture(*(this->pair_program),"sw3_tex");
  sw3_tex.bind_float(sw3,4);

  UCL_H_Vec<int> dview_elem2param(nelements*nelements*nelements,
                           *(this->ucl_device), UCL_WRITE_ONLY);

  elem2param.alloc(nelements*nelements*nelements,*(this->ucl_device),
                   UCL_READ_ONLY);

  for (int i = 0; i < nelements; i++)
    for (int j = 0; j < nelements; j++)
      for (int k = 0; k < nelements; k++) {
         int idx = i*nelements*nelements+j*nelements+k;
         dview_elem2param[idx] = host_elem2param[i][j][k];
      }

  ucl_copy(elem2param,dview_elem2param,false);

  UCL_H_Vec<int> dview_map(lj_types, *(this->ucl_device), UCL_WRITE_ONLY);
  for (int i = 0; i < ntypes; i++)
    dview_map[i] = host_map[i];

  map.alloc(lj_types,*(this->ucl_device), UCL_READ_ONLY);
  ucl_copy(map,dview_map,false);

  _allocated=true;
  this->_max_bytes=sw1.row_bytes()+sw2.row_bytes()+sw3.row_bytes()+
    map.row_bytes()+elem2param.row_bytes();
  return 0;
}

template <class numtyp, class acctyp>
void SWT::clear() {
  if (!_allocated)
    return;
  _allocated=false;

  sw1.clear();
  sw2.clear();
  sw3.clear();
  map.clear();
  elem2param.clear();
  this->clear_atomic();
}

template <class numtyp, class acctyp>
double SWT::host_memory_usage() const {
  return this->host_memory_usage_atomic()+sizeof(SW<numtyp,acctyp>);
}

#define KTHREADS this->_threads_per_atom
#define JTHREADS this->_threads_per_atom
// ---------------------------------------------------------------------------
// Calculate energies, forces, and torques
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void SWT::loop(const bool _eflag, const bool _vflag, const int evatom) {
  // Compute the block size and grid size to keep all cores busy
  int BX=this->block_pair();
  int eflag, vflag;
  if (_eflag)
    eflag=1;
  else
    eflag=0;

  if (_vflag)
    vflag=1;
  else
    vflag=0;

  // build the short neighbor list
  int ainum=this->_ainum;
  int nbor_pitch=this->nbor->nbor_pitch();
  int GX=static_cast<int>(ceil(static_cast<double>(ainum)/
                               (BX/this->_threads_per_atom)));
  this->k_short_nbor.set_size(GX,BX);
  this->k_short_nbor.run(&this->atom->x, &sw3, &map, &elem2param, &_nelements,
                 &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                 &this->dev_short_nbor, &ainum,
                 &nbor_pitch, &this->_threads_per_atom);

  // this->_nbor_data == nbor->dev_packed for gpu_nbor == 0 and tpa > 1
  // this->_nbor_data == nbor->dev_nbor for gpu_nbor == 1 or tpa == 1
  ainum=this->ans->inum();
  nbor_pitch=this->nbor->nbor_pitch();
  GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/
                               (BX/this->_threads_per_atom)));
  this->time_pair.start();
  
  this->k_pair.set_size(GX,BX);
  this->k_pair.run(&this->atom->x, &sw1, &sw2, &sw3,
                   &map, &elem2param, &_nelements,
                   &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                   &this->dev_short_nbor,
                   &this->ans->force, &this->ans->engv,
                   &eflag, &vflag, &ainum, &nbor_pitch,
                   &this->_threads_per_atom);

  BX=this->block_size();
  GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/
                           (BX/(KTHREADS*JTHREADS))));
  this->k_three_center.set_size(GX,BX);
  this->k_three_center.run(&this->atom->x, &sw1, &sw2, &sw3,
                           &map, &elem2param, &_nelements,
                           &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                           &this->dev_short_nbor,
                           &this->ans->force, &this->ans->engv, &eflag, &vflag, &ainum,
                           &nbor_pitch, &this->_threads_per_atom, &evatom);

  Answer<numtyp,acctyp> *end_ans;
  #ifdef THREE_CONCURRENT
  end_ans=this->ans2;
  #else
  end_ans=this->ans;
  #endif
  if (evatom!=0) {
    this->k_three_end_vatom.set_size(GX,BX);
    this->k_three_end_vatom.run(&this->atom->x, &sw1, &sw2, &sw3,
                          &map, &elem2param, &_nelements,
                          &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                          &this->nbor->dev_ilist, &this->dev_short_nbor,
                          &end_ans->force, &end_ans->engv, &eflag, &vflag, &ainum,
                          &nbor_pitch, &this->_threads_per_atom, &this->_gpu_nbor);

  } else {
    this->k_three_end.set_size(GX,BX);
    this->k_three_end.run(&this->atom->x, &sw1, &sw2, &sw3,
                          &map, &elem2param, &_nelements,
                          &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                          &this->nbor->dev_ilist, &this->dev_short_nbor,
                          &end_ans->force, &end_ans->engv, &eflag, &vflag, &ainum,
                          &nbor_pitch, &this->_threads_per_atom, &this->_gpu_nbor);

  }

  this->time_pair.stop();
}

template class SW<PRECISION,ACC_PRECISION>;

