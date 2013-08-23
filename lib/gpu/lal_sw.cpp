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
const char *lj=0;
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
int SWT::init(const int nlocal, const int nall, const int max_nbors, 
              const double cell_size, const double gpu_split, FILE *_screen,
              const double epsilon, const double sigma,
              const double lambda, const double gamma,
              const double costheta, const double biga,
              const double bigb, const double powerp,
              const double powerq, const double cut, const double cutsq) {

  sw_epsilon=static_cast<numtyp>(epsilon);
  sw_sigma=static_cast<numtyp>(sigma);
  sw_lambda=static_cast<numtyp>(lambda);
  sw_gamma=static_cast<numtyp>(gamma);
  sw_costheta=static_cast<numtyp>(costheta);
  sw_biga=static_cast<numtyp>(biga);
  sw_bigb=static_cast<numtyp>(bigb);
  sw_powerp=static_cast<numtyp>(powerp);
  sw_powerq=static_cast<numtyp>(powerq);
  sw_cut=static_cast<numtyp>(cut);
  sw_cutsq=static_cast<numtyp>(cutsq);
  if (sw_cutsq>=sw_cut*sw_cut) 
    sw_cutsq=sw_cut*sw_cut-1e-4;

  int success;
  success=this->init_three(nlocal,nall,max_nbors,0,cell_size,gpu_split,
                           _screen,sw,"k_sw","k_sw_three_center",
                           "k_sw_three_end");
  if (success!=0)
    return success;

  // If atom type constants fit in shared memory use fast kernel
  shared_types=true;

  _allocated=true;
  this->_max_bytes=0;
  return 0;
}

template <class numtyp, class acctyp>
void SWT::clear() {
  if (!_allocated)
    return;
  _allocated=false;

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
  
  int GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/
                               (BX/this->_threads_per_atom)));

  int ainum=this->ans->inum();
  int nbor_pitch=this->nbor->nbor_pitch();
  this->time_pair.start();
  this->k_pair.set_size(GX,BX);
  this->k_pair.run(&this->atom->x, &this->nbor->dev_nbor, 
                   &this->_nbor_data->begin(), &this->ans->force,
                   &this->ans->engv, &eflag, &vflag, &ainum, &nbor_pitch, 
                   &this->_threads_per_atom, &sw_cut, &sw_epsilon, &sw_sigma,
                   &sw_biga, &sw_bigb, &sw_powerp, &sw_powerq, &sw_cutsq);

  BX=this->block_size();
  GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/
                           (BX/(KTHREADS*JTHREADS)))); 
  this->k_three_center.set_size(GX,BX);
  this->k_three_center.run(&this->atom->x, &this->nbor->dev_nbor, 
                   &this->_nbor_data->begin(), &this->ans->force,
                   &this->ans->engv, &eflag, &vflag, &ainum, 
                   &nbor_pitch, &this->_threads_per_atom, &evatom,
                   &sw_cut, &sw_epsilon, &sw_sigma, &sw_lambda, &sw_gamma,
                   &sw_costheta, &sw_cutsq);
  Answer<numtyp,acctyp> *end_ans;
  #ifdef THREE_CONCURRENT
  end_ans=this->ans2;
  #else
  end_ans=this->ans;
  #endif
  if (evatom!=0) {
    this->k_three_end_vatom.set_size(GX,BX);
    this->k_three_end_vatom.run(&this->atom->x, &this->nbor->dev_nbor, 
                          &this->_nbor_data->begin(), &end_ans->force,
                          &end_ans->engv, &eflag, &vflag, &ainum, 
                          &nbor_pitch, &this->_threads_per_atom, &sw_cut, 
                          &sw_epsilon, &sw_sigma, &sw_lambda, &sw_gamma,
                          &sw_costheta, &sw_cutsq);
  } else {
    this->k_three_end.set_size(GX,BX);
    this->k_three_end.run(&this->atom->x, &this->nbor->dev_nbor, 
                          &this->_nbor_data->begin(), &end_ans->force,
                          &end_ans->engv, &eflag, &vflag, &ainum, 
                          &nbor_pitch, &this->_threads_per_atom, &sw_cut, 
                          &sw_epsilon, &sw_sigma, &sw_lambda, &sw_gamma,
                          &sw_costheta, &sw_cutsq);
  }
  this->time_pair.stop();
}

template class SW<PRECISION,ACC_PRECISION>;

