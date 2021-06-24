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
namespace LAMMPS_AL {
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
int SWT::init(const int ntypes, const int nlocal, const int nall,
              const int max_nbors, const double cell_size,
              const double gpu_split, FILE *_screen, double **ncutsq,
              double **ncut, double **sigma, double **powerp, double **powerq,
              double **sigma_gamma, double **c1, double **c2, double **c3,
              double **c4, double **c5, double **c6, double ***lambda_epsilon,
              double ***costheta, const int *map, int ***e2param) {
  _lj_types=ntypes;

  int oldparam=-1;
  int onetype=-1;
  int onetype3=0;
  int spq=1;
  int mtypes=0;
  #ifdef USE_OPENCL
  for (int ii=1; ii<ntypes; ii++) {
    int i=map[ii];
    if (i<0) continue;
    for (int jj=1; jj<ntypes; jj++) {
      int j=map[jj];
      if (j<0) continue;
      if (powerp[ii][jj] != 4.0 || powerq[ii][jj] != 0.0)
        spq=0;
      for (int kk=1; kk<ntypes; kk++) {
        int k=map[kk];
        if (k<0) continue;
        int param=e2param[i][j][k];
        if (oldparam!=param) {
          oldparam=param;
          onetype=ntypes*ii+jj;
          onetype3=ntypes*ntypes*ii+ntypes*jj+kk;
          mtypes++;
        }
      }
    }
  }
  if (mtypes>1) onetype=-1;
  #endif

  int success;
  success=this->init_three(nlocal,nall,max_nbors,0,cell_size,gpu_split,
                           _screen,sw,"k_sw","k_sw_three_center",
                           "k_sw_three_end","k_sw_short_nbor",onetype,
                           onetype3,spq);
  if (success!=0)
    return success;

  UCL_H_Vec<numtyp> host_write(ntypes*ntypes*ntypes*4,*(this->ucl_device),
                               UCL_WRITE_ONLY);
  host_write.zero();

  for (int i=1; i<ntypes; i++)
    for (int j=1; j<ntypes; j++) {
      double ccutsq = ncut[i][j]*ncut[i][j];
      if (ccutsq > 0.0 && ncutsq[i][j]>=ccutsq)
        ncutsq[i][j]=ccutsq*0.98;
  }

  // pack coefficients into arrays
  cutsq.alloc(ntypes*ntypes,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack1(ntypes,ntypes,cutsq,host_write,ncutsq);
  sw_pre.alloc(ntypes*ntypes,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack4(ntypes,ntypes,sw_pre,host_write,ncut,sigma,
                         powerp,powerq);
  c_14.alloc(ntypes*ntypes,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack4(ntypes,ntypes,c_14,host_write,c1,c2,c3,c4);
  c_56.alloc(ntypes*ntypes,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack2(ntypes,ntypes,c_56,host_write,c5,c6);
  cut_sigma_gamma.alloc(ntypes*ntypes,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack2(ntypes,ntypes,cut_sigma_gamma,host_write,ncut,
                         sigma_gamma);
  sw_pre3.alloc(ntypes*ntypes*ntypes,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack2(ntypes,sw_pre3,host_write,lambda_epsilon,costheta);

  _allocated=true;
  this->_max_bytes=cutsq.row_bytes()+sw_pre.row_bytes()+c_14.row_bytes()+
    c_56.row_bytes()+cut_sigma_gamma.row_bytes()+sw_pre3.row_bytes();
  return 0;
}

template <class numtyp, class acctyp>
void SWT::clear() {
  if (!_allocated)
    return;
  _allocated=false;

  cutsq.clear();
  sw_pre.clear();
  c_14.clear();
  c_56.clear();
  cut_sigma_gamma.clear();
  sw_pre3.clear();
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
int SWT::loop(const int eflag, const int vflag, const int evatom,
              bool &success) {
  const int nbor_pitch=this->nbor->nbor_pitch();

  // build the short neighbor list
  int ainum=this->_ainum;
  this->time_pair.start();

  int BX=this->block_pair();
  int GX=static_cast<int>(ceil(static_cast<double>(ainum)/BX));
  this->k_short_nbor.set_size(GX,BX);
  this->k_short_nbor.run(&this->atom->x, &cutsq, &_lj_types,
                         &this->nbor->dev_nbor, &this->nbor->dev_packed,
                         &ainum, &nbor_pitch, &this->_threads_per_atom);

  // this->_nbor_data == nbor->dev_packed for gpu_nbor == 0 and tpa > 1
  // this->_nbor_data == nbor->dev_nbor for gpu_nbor == 1 or tpa == 1
  ainum=this->ans->inum();
  BX=this->block_size();
  GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/
                           (BX/(KTHREADS*JTHREADS))));
  this->k_3center_sel->set_size(GX,BX);
  this->k_3center_sel->run(&this->atom->x, &cut_sigma_gamma, &sw_pre3,
                           &_lj_types, &this->nbor->dev_nbor,
                           &this->ans->force, &this->ans->engv, &eflag,
                           &vflag, &ainum, &nbor_pitch,
                           &this->_threads_per_atom, &evatom);

  Answer<numtyp,acctyp> *end_ans;
  #ifdef THREE_CONCURRENT
  end_ans=this->ans2;
  #else
  end_ans=this->ans;
  #endif
  if (evatom!=0) {
    this->k_three_end_vatom.set_size(GX,BX);
    this->k_three_end_vatom.run(&this->atom->x, &cut_sigma_gamma,
                                &sw_pre3, &_lj_types, &this->nbor->dev_nbor,
                                &this->nbor->three_ilist, &end_ans->force,
                                &end_ans->engv, &eflag, &vflag, &ainum,
                                &nbor_pitch,&this->_threads_per_atom,
                                &this->_gpu_nbor);
  } else {
    this->k_3end_sel->set_size(GX,BX);
    this->k_3end_sel->run(&this->atom->x, &cut_sigma_gamma, &sw_pre3,
                          &_lj_types, &this->nbor->dev_nbor,
                          &this->nbor->three_ilist, &end_ans->force,
                          &end_ans->engv, &eflag, &vflag, &ainum, &nbor_pitch,
                          &this->_threads_per_atom, &this->_gpu_nbor);
  }

  BX=this->block_pair();
  int GXT=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/
                               (BX/this->_threads_per_atom)));
  this->k_sel->set_size(GXT,BX);
  this->k_sel->run(&this->atom->x, &sw_pre, &c_14, &c_56,
                   &_lj_types, &this->nbor->dev_nbor,
                   &this->ans->force, &this->ans->engv, &eflag, &vflag,
                   &ainum, &nbor_pitch, &this->_threads_per_atom, &GX);

  this->time_pair.stop();
  return GX;
}

template class SW<PRECISION,ACC_PRECISION>;
}
