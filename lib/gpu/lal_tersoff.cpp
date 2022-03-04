/***************************************************************************
                                 tersoff.cpp
                             -------------------
                               Trung Dac Nguyen

  Class for acceleration of the tersoff pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : Thu April 17, 2014
    email                : ndactrung@gmail.com
 ***************************************************************************/

#if defined(USE_OPENCL)
#include "tersoff_cl.h"
#elif defined(USE_CUDART)
const char *tersoff=0;
#else
#include "tersoff_cubin.h"
#endif

#include "lal_tersoff.h"
#include <cassert>
namespace LAMMPS_AL {
#define TersoffT Tersoff<numtyp, acctyp>

extern Device<PRECISION,ACC_PRECISION> device;

template <class numtyp, class acctyp>
TersoffT::Tersoff() : BaseThree<numtyp,acctyp>(), _allocated(false) {
}

template <class numtyp, class acctyp>
TersoffT::~Tersoff() {
  clear();
}

template <class numtyp, class acctyp>
int TersoffT::bytes_per_atom(const int max_nbors) const {
  return this->bytes_per_atom_atomic(max_nbors)+max_nbors*sizeof(acctyp)*4;
}

template <class numtyp, class acctyp>
int TersoffT::init(const int ntypes, const int nlocal, const int nall, const int max_nbors,
                   const double cell_size, const double gpu_split, FILE *_screen,
                   int* host_map, const int nelements, int*** host_elem2param, const int nparams,
                   const double* lam1, const double* lam2, const double* lam3,const double* powermint,
                   const double* biga, const double* bigb, const double* bigr, const double* bigd,
                   const double* c1, const double* c2, const double* c3, const double* c4,
                   const double* c, const double* d, const double* h, const double* gamma,
                   const double* beta, const double* powern, const double* host_cutsq)
{
  int oldparam=-1;
  int onetype=-1;
  int onetype3=0;
  int spq=0;
  int mtypes=0;
  #ifdef USE_OPENCL
  for (int ii=1; ii<ntypes; ii++) {
    const int i=host_map[ii];
    for (int jj=1; jj<ntypes; jj++) {
      const int j=host_map[jj];
      for (int kk=1; kk<ntypes; kk++) {
        const int k=host_map[kk];
        if (i<0 || j<0 || k<0) continue;
        const int ijkparam = host_elem2param[i][j][k];
        if (oldparam!=ijkparam) {
          oldparam=ijkparam;
          onetype=ntypes*ii+jj;
          onetype3=ijkparam;
          mtypes++;
        }
      }
    }
  }
  if (mtypes>1) onetype=-1;
  if (onetype>=0) spq=powermint[onetype3];
  #endif

  int success;
  success=this->init_three(nlocal,nall,max_nbors,0,cell_size,gpu_split,
                           _screen,tersoff,"k_tersoff_repulsive",
                           "k_tersoff_three_center", "k_tersoff_three_end",
                           "k_tersoff_short_nbor",onetype,onetype3,spq,1);
  if (success!=0)
    return success;

  int ef_nall=nall;
  if (ef_nall==0)
    ef_nall=2000;
  if (this->nbor->max_nbors()) {
    _zetaij.alloc(ef_nall*this->nbor->max_nbors(),*(this->ucl_device),
                  UCL_READ_WRITE);
    _zetaij_eng.alloc(ef_nall*this->nbor->max_nbors(),*(this->ucl_device),
                      UCL_READ_WRITE);
  }

  k_zeta.set_function(*(this->pair_program),"k_tersoff_zeta");
  #if defined(LAL_OCL_EV_JIT)
  k_zeta_noev.set_function(*(this->pair_program_noev),"k_tersoff_zeta");
  #else
  k_zeta_selt = &k_zeta;
  #endif

  _ntypes=ntypes;
  _nparams = nparams;
  _nelements = nelements;

  _cutsq_max=0.0;
  for (int ii=1; ii<ntypes; ii++) {
    const int i=host_map[ii];
    for (int jj=1; jj<ntypes; jj++) {
      const int j=host_map[jj];
      for (int kk=1; kk<ntypes; kk++) {
        const int k=host_map[kk];
        if (i<0 || j<0 || k<0) continue;
        const int ijkparam = host_elem2param[i][j][k];
        if (host_cutsq[ijkparam]>_cutsq_max) _cutsq_max=host_cutsq[ijkparam];
      }
    }
  }

  // --------------------------------------------------------------------
  UCL_H_Vec<numtyp4> dview(nparams,*(this->ucl_device),
                           UCL_WRITE_ONLY);

  for (int i=0; i<nparams; i++) {
    dview[i].x=(numtyp)0;
    dview[i].y=(numtyp)0;
    dview[i].z=(numtyp)0;
    dview[i].w=(numtyp)0;
  }

  // pack coefficients into arrays
  // pack coefficients into arrays
  ts1.alloc(nparams,*(this->ucl_device),UCL_READ_ONLY);

  for (int i=0; i<nparams; i++) {
    dview[i].x=static_cast<numtyp>(lam3[i]);
    dview[i].y=static_cast<numtyp>(powermint[i]);
    dview[i].z=static_cast<numtyp>(bigr[i]);
    dview[i].w=static_cast<numtyp>(bigd[i]);
  }

  ucl_copy(ts1,dview,false);

  ts2.alloc(nparams,*(this->ucl_device),UCL_READ_ONLY);

  for (int i=0; i<nparams; i++) {
    dview[i].x=static_cast<numtyp>(biga[i]);
    dview[i].y=static_cast<numtyp>(lam1[i]);
    dview[i].z=static_cast<numtyp>(bigr[i]);
    dview[i].w=static_cast<numtyp>(bigd[i]);
  }

  ucl_copy(ts2,dview,false);

  ts3.alloc(nparams,*(this->ucl_device),UCL_READ_ONLY);

  for (int i=0; i<nparams; i++) {
    dview[i].x=static_cast<numtyp>(c1[i]);
    dview[i].y=static_cast<numtyp>(c2[i]);
    dview[i].z=static_cast<numtyp>(c3[i]);
    dview[i].w=static_cast<numtyp>(c4[i]);
  }

  ucl_copy(ts3,dview,false);

  ts4.alloc(nparams,*(this->ucl_device),UCL_READ_ONLY);

  for (int i=0; i<nparams; i++) {
    dview[i].x=static_cast<numtyp>(c[i]*c[i]);
    dview[i].y=static_cast<numtyp>(d[i]*d[i]);
    dview[i].z=static_cast<numtyp>(h[i]);
    dview[i].w=static_cast<numtyp>(gamma[i]);
  }

  ucl_copy(ts4,dview,false);

  ts5.alloc(nparams,*(this->ucl_device),UCL_READ_ONLY);

  for (int i=0; i<nparams; i++) {
    dview[i].x=static_cast<numtyp>(beta[i]);
    dview[i].y=static_cast<numtyp>(powern[i]);
    dview[i].z=static_cast<numtyp>(lam2[i]);
    dview[i].w=static_cast<numtyp>(bigb[i]);
  }

  ucl_copy(ts5,dview,false);

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

  UCL_H_Vec<int> dview_map(ntypes, *(this->ucl_device), UCL_WRITE_ONLY);
  for (int i = 0; i < ntypes; i++)
    dview_map[i] = host_map[i];

  map.alloc(ntypes,*(this->ucl_device), UCL_READ_ONLY);
  ucl_copy(map,dview_map,false);

  _allocated=true;
  this->_max_bytes=ts1.row_bytes()+ts2.row_bytes()+ts3.row_bytes()+
    ts4.row_bytes()+ts5.row_bytes()+map.row_bytes()+
    elem2param.row_bytes()+_zetaij.row_bytes()+_zetaij_eng.row_bytes();
  return 0;
}

template <class numtyp, class acctyp>
void TersoffT::clear() {
  if (!_allocated)
    return;
  _allocated=false;

  ts1.clear();
  ts2.clear();
  ts3.clear();
  ts4.clear();
  ts5.clear();
  map.clear();
  elem2param.clear();
  _zetaij.clear();
  _zetaij_eng.clear();

  k_zeta.clear();
  #if defined(LAL_OCL_EV_JIT)
  k_zeta_noev.clear();
  #endif

  this->clear_atomic();
}

template <class numtyp, class acctyp>
double TersoffT::host_memory_usage() const {
  return this->host_memory_usage_atomic()+sizeof(Tersoff<numtyp,acctyp>);
}

#define KTHREADS this->_threads_per_atom
#define JTHREADS this->_threads_per_atom
// ---------------------------------------------------------------------------
// Calculate energies, forces, and torques
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int TersoffT::loop(const int eflag, const int vflag, const int evatom,
                   bool &success) {
  const int nbor_pitch=this->nbor->nbor_pitch();

  // re-allocate zetaij if necessary
  int nall = this->_nall;
  if (nall*this->nbor->max_nbors() > _zetaij.cols()) {
    int _nmax=static_cast<int>(static_cast<double>(nall)*1.10);
    _zetaij.clear();
    _zetaij_eng.clear();
    success = success && (_zetaij.alloc(this->nbor->max_nbors()*_nmax,
                                        *(this->ucl_device),
                                        UCL_READ_WRITE) == UCL_SUCCESS);
    success = success && (_zetaij_eng.alloc(this->nbor->max_nbors()*_nmax,
                                            *(this->ucl_device),
                                            UCL_READ_WRITE) == UCL_SUCCESS);
    if (!success) return 0;
  }

  // build the short neighbor list
  int ainum=this->_ainum;
  this->time_pair.start();

  int BX=this->block_pair();
  int GX=static_cast<int>(ceil(static_cast<double>(ainum)/BX));
  this->k_short_nbor.set_size(GX,BX);
  this->k_short_nbor.run(&this->atom->x, &_cutsq_max, &_ntypes,
                         &this->nbor->dev_nbor, &this->nbor->dev_packed,
                         &ainum, &nbor_pitch, &this->_threads_per_atom);

  #if defined(LAL_OCL_EV_JIT)
  if (eflag || vflag) k_zeta_selt = &k_zeta;
  else k_zeta_selt = &k_zeta_noev;
  #endif

  GX=static_cast<int>(ceil(static_cast<double>(this->_ainum)/
                               (BX/(JTHREADS*KTHREADS))));
  k_zeta_selt->set_size(GX,BX);
  k_zeta_selt->run(&this->atom->x, &ts1, &ts3, &ts4, &ts5,
                   &map, &elem2param, &_nelements, &_nparams, &_zetaij,
                   &_zetaij_eng, &this->nbor->dev_nbor, &eflag, &this->_ainum,
                   &nbor_pitch, &this->_threads_per_atom);

  ainum=this->ans->inum();
  BX=this->block_size();
  GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/
                           (BX/(KTHREADS*JTHREADS))));
  this->k_3center_sel->set_size(GX,BX);
  this->k_3center_sel->run(&this->atom->x, &ts1, &ts4, &map,
                           &elem2param, &_nelements, &_nparams, &_zetaij,
                           &_zetaij_eng, &this->nbor->dev_nbor,
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
    this->k_three_end_vatom.run(&this->atom->x, &ts1, &ts4, &map, &elem2param,
                          &_nelements, &_nparams, &_zetaij, &_zetaij_eng,
                          &this->nbor->dev_nbor, &this->nbor->three_ilist,
                          &end_ans->force, &end_ans->engv, &eflag, &vflag,
                          &ainum, &nbor_pitch, &this->_threads_per_atom,
                          &this->_gpu_nbor);

  } else {
    this->k_3end_sel->set_size(GX,BX);
    this->k_3end_sel->run(&this->atom->x, &ts1, &ts4, &map, &elem2param,
                          &_nelements, &_nparams, &_zetaij, &_zetaij_eng,
                          &this->nbor->dev_nbor, &this->nbor->three_ilist,
                          &end_ans->force, &end_ans->engv, &eflag, &vflag,
                          &ainum, &nbor_pitch, &this->_threads_per_atom,
                          &this->_gpu_nbor);
  }

  BX=this->block_pair();
  int GXT=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/
                               (BX/this->_threads_per_atom)));
  this->k_sel->set_size(GXT,BX);
  this->k_sel->run(&this->atom->x, &ts2, &map, &elem2param, &_nelements,
                   &_nparams, &this->nbor->dev_nbor, &this->ans->force,
                   &this->ans->engv, &eflag, &vflag, &ainum, &nbor_pitch,
                   &this->_threads_per_atom, &GX);

  this->time_pair.stop();
  return GX;
}

template class Tersoff<PRECISION,ACC_PRECISION>;
}
