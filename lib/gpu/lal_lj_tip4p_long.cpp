/**************************************************************************
                               lj_tip4p_long.cpp
                             -------------------
                              V. Nikolskiy (HSE)

  Class for acceleration of the lj/tip4p/long pair style

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : thevsevak@gmail.com
***************************************************************************/

#if defined(USE_OPENCL)
#include "lj_tip4p_long_cl.h"
#elif defined(USE_CUDART)
const char *lj_tip4p=0;
#else
#include "lj_tip4p_long_cubin.h"
#endif

#include "lal_lj_tip4p_long.h"
#include <cassert>
using namespace LAMMPS_AL;
#define LJTIP4PLongT LJ_TIP4PLong<numtyp, acctyp>

extern Device<PRECISION,ACC_PRECISION> device;

template <class numtyp, class acctyp>
LJTIP4PLongT::LJ_TIP4PLong(): BaseCharge<numtyp,acctyp>(), _allocated(false) {
}

template <class numtyp, class acctyp>
LJTIP4PLongT::~LJ_TIP4PLong() {
  clear();
}

template <class numtyp, class acctyp>
int LJTIP4PLongT::bytes_per_atom(const int max_nbors) const {
  return this->bytes_per_atom_atomic(max_nbors);
}

template <class numtyp, class acctyp>
int LJTIP4PLongT::init(const int ntypes,
    double **host_cutsq, double **host_lj1,
    double **host_lj2, double **host_lj3,
    double **host_lj4, double **host_offset,
    double *host_special_lj, const int nlocal,
    const int tH, const int tO,
    const double a, const double qd,
    const int nall, const int max_nbors,
    const int maxspecial, const double cell_size,
    const double gpu_split, FILE *_screen,
    double **host_cut_ljsq,
    const double host_cut_coulsq, const double host_cut_coulsqplus,
    double *host_special_coul, const double qqrd2e,
    const double g_ewald, int map_size, int max_same) {
  int success;
  success=this->init_atomic(nlocal,nall,max_nbors,maxspecial,cell_size,gpu_split,
                            _screen,lj_tip4p_long,"k_lj_tip4p_long");
  if (success!=0)
    return success;
  k_pair_distrib.set_function(*this->pair_program,"k_lj_tip4p_long_distrib");
  k_pair_reneigh.set_function(*this->pair_program,"k_lj_tip4p_reneigh");
  k_pair_newsite.set_function(*this->pair_program,"k_lj_tip4p_newsite");

  TypeH = tH;
  TypeO = tO;
  alpha = a;
  qdist = qd;

  // If atom type constants fit in shared memory use fast kernel
  int lj_types=ntypes;
  shared_types=false;
//  int max_shared_types=this->device->max_shared_types();
//  if (lj_types<=max_shared_types && this->_block_size>=max_shared_types) {
//    lj_types=max_shared_types;
//    shared_types=true;
//  }
  _lj_types=lj_types;

  // Allocate a host write buffer for data initialization
  UCL_H_Vec<numtyp> host_write(lj_types*lj_types*32,*(this->ucl_device),
                               UCL_WRITE_ONLY);

  for (int i=0; i<lj_types*lj_types; i++)
    host_write[i]=0.0;

  lj1.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack4(ntypes,lj_types,lj1,host_write,host_lj1,host_lj2,
                         host_cut_ljsq);

  lj3.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack4(ntypes,lj_types,lj3,host_write,host_lj3,host_lj4,
                         host_offset);

  cutsq.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack1(ntypes,lj_types,cutsq,host_write,host_cutsq);

  sp_lj.alloc(8,*(this->ucl_device),UCL_READ_ONLY);
  for (int i=0; i<4; i++) {
    host_write[i]=host_special_lj[i];
    host_write[i+4]=host_special_coul[i];
  }
  ucl_copy(sp_lj,host_write,8,false);

  //force_comp.alloc(72*72, *(this->ucl_device), UCL_READ_WRITE);

  _qqrd2e=qqrd2e;
  _g_ewald=g_ewald;
  cut_coulsq = host_cut_coulsq;
  cut_coulsqplus = host_cut_coulsqplus;

  hneight.alloc(nall*4,*(this->ucl_device), UCL_READ_WRITE);
  m.alloc(nall,*(this->ucl_device), UCL_READ_WRITE);
  ansO.alloc(nall,*(this->ucl_device), UCL_READ_WRITE);

  this->tag.alloc(nall,*(this->ucl_device), UCL_READ_ONLY);
  this->atom_sametag.alloc(max_same, *(this->ucl_device), UCL_READ_ONLY);
  this->map_array.alloc(map_size,*(this->ucl_device), UCL_READ_ONLY);

  _allocated=true;
  this->_max_bytes=lj1.row_bytes()+lj3.row_bytes()+cutsq.row_bytes()+
      sp_lj.row_bytes() + hneight.row_bytes()+m.row_bytes()+
      this->tag.row_bytes()+this->atom_sametag.row_bytes() +
      this->map_array.row_bytes();
  return 0;
}


template <class numtyp, class acctyp>
void LJTIP4PLongT::clear() {
  if (!_allocated)
    return;
  _allocated=false;

  lj1.clear();
  lj3.clear();
  sp_lj.clear();
  cutsq.clear();
  hneight.clear();
  m.clear();
  tag.clear();
  atom_sametag.clear();
  map_array.clear();
  ansO.clear();
  //force_comp.clear();

  k_pair_distrib.clear();
  k_pair_reneigh.clear();
  k_pair_newsite.clear();

  this->clear_atomic();
}

template <class numtyp, class acctyp>
double LJTIP4PLongT::host_memory_usage() const {
  return this->host_memory_usage_atomic()+sizeof(LJ_TIP4PLong<numtyp,acctyp>);
}

// ---------------------------------------------------------------------------
// Calculate energies, forces, and torques
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void LJTIP4PLongT::loop(const bool _eflag, const bool _vflag) {
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

  int ainum=this->ans->inum();
  const int nall = this->atom->nall();
  int nbor_pitch=this->nbor->nbor_pitch();
  this->time_pair.start();
  int GX;
  GX=static_cast<int>(ceil(static_cast<double>(nall)/BX));
  if (t_ago == 0) {
    this->k_pair_reneigh.set_size(GX,BX);
    this->k_pair_reneigh.run(&this->atom->x,
        &this->nbor->dev_nbor, &this->_nbor_data->begin(),
        &nall, &ainum,&nbor_pitch, &this->_threads_per_atom,
        &hneight, &m, &TypeO, &TypeH,
        &tag, &map_array, &atom_sametag);
  }
  this->k_pair_newsite.set_size(GX,BX);
  this->k_pair_newsite.run(&this->atom->x,
          &this->nbor->dev_nbor, &this->_nbor_data->begin(),
          &nall, &ainum,
          &nbor_pitch, &this->_threads_per_atom,
          &hneight, &m, &TypeO, &TypeH, &alpha,
          &this->atom->q, &tag, &map_array,
          &atom_sametag);

  GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/
                               (BX/this->_threads_per_atom)));
  this->k_pair.set_size(GX,BX);
  if (vflag){
	  this->ansO.resize_ib(ainum*3);
  } else {
	  this->ansO.resize_ib(ainum);
  }
  this->ansO.zero();
  this->device->gpu->sync();
  this->k_pair.run(&this->atom->x, &lj1, &lj3, &_lj_types, &sp_lj,
          &this->nbor->dev_nbor, &this->_nbor_data->begin(),
          &this->ans->force, &this->ans->engv, &eflag, &vflag,
          &ainum, &nbor_pitch, &this->_threads_per_atom,
          &hneight, &m, &TypeO, &TypeH, &alpha,
          &this->atom->q, &cutsq, &_qqrd2e, &_g_ewald,
          &cut_coulsq, &cut_coulsqplus, &tag, &map_array,
          &atom_sametag, &this->ansO);
  GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/BX));
  this->k_pair_distrib.set_size(GX,BX);
  this->k_pair_distrib.run(&this->atom->x, &this->ans->force, &this->ans->engv,
      &eflag, &vflag, &ainum, &nbor_pitch, &this->_threads_per_atom,
      &hneight, &m, &TypeO, &TypeH, &alpha,&this->atom->q,  &this->ansO);
  this->time_pair.stop();
}


template <class numtyp, class acctyp>
void LJTIP4PLongT::copy_relations_data(int n, int* tag, int *map_array,
                      int map_size, int *sametag, int max_same, int ago){
  int nall = n;
  const int hn_sz = n*4; // matrix size = col size * col number
  hneight.resize_ib(hn_sz);
  if (ago == 0)
    hneight.zero();
  m.resize_ib(n);
  m.zero();

  UCL_H_Vec<int> host_tag_write(nall,*(this->ucl_device),UCL_WRITE_ONLY);
  this->tag.resize_ib(nall);
  for(int i=0; i<nall; ++i) host_tag_write[i] = tag[i];
  ucl_copy(this->tag, host_tag_write, nall, false);

  host_tag_write.resize_ib(max_same);
  this->atom_sametag.resize_ib(max_same);
  for(int i=0; i<max_same; ++i) host_tag_write[i] = sametag[i];
  ucl_copy(this->atom_sametag, host_tag_write, max_same, false);

  host_tag_write.resize_ib(map_size);
  this->map_array.resize_ib(map_size);
  for(int i=0; i<map_size; ++i) host_tag_write[i] = map_array[i];
  ucl_copy(this->map_array, host_tag_write, map_size, false);
}




// ---------------------------------------------------------------------------
// Copy nbor list from host if necessary and then calculate forces, virials,..
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void LJTIP4PLongT::compute(const int f_ago, const int inum_full,
                               const int nall, double **host_x, int *host_type,
                               int *ilist, int *numj, int **firstneigh,
                               const bool eflag, const bool vflag,
                               const bool eatom, const bool vatom,
                               int &host_start, const double cpu_time,
                               bool &success, double *host_q,
                               const int nlocal, double *boxlo, double *prd) {
  this->acc_timers();
  if (inum_full==0) {
    host_start=0;
    // Make sure textures are correct if realloc by a different hybrid style
    this->resize_atom(0,nall,success);
    this->zero_timers();
    return;
  }

  int ago=this->hd_balancer.ago_first(f_ago);
  int inum=this->hd_balancer.balance(ago,inum_full,cpu_time);
  this->ans->inum(inum);
  host_start=inum;

  if (ago==0) {
    this->reset_nbors(nall, inum, ilist, numj, firstneigh, success);
    if (!success)
      return;
  }

  this->atom->cast_x_data(host_x,host_type);
  this->atom->cast_q_data(host_q);
  this->hd_balancer.start_timer();
  this->atom->add_x_data(host_x,host_type);
  this->atom->add_q_data();

  this->device->precompute(f_ago,nlocal,nall,host_x,host_type,success,host_q,
                     boxlo, prd);

  t_ago = ago;
  loop(eflag,vflag);
  this->ans->copy_answers(eflag,vflag,eatom,vatom,ilist);
  this->device->add_ans_object(this->ans);
  this->hd_balancer.stop_timer();
}

// ---------------------------------------------------------------------------
// Reneighbor on GPU if necessary and then compute forces, virials, energies
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int** LJTIP4PLongT::compute(const int ago, const int inum_full,
                                const int nall, double **host_x, int *host_type,
                                double *sublo, double *subhi, tagint *tag,
                                int *map_array, int map_size, int *sametag, int max_same,
                                int **nspecial, tagint **special, const bool eflag,
                                const bool vflag, const bool eatom,
                                const bool vatom, int &host_start,
                                int **ilist, int **jnum,
                                const double cpu_time, bool &success,
                                double *host_q, double *boxlo, double *prd) {
  this->acc_timers();
  if (inum_full==0) {
    host_start=0;
    // Make sure textures are correct if realloc by a different hybrid style
    this->resize_atom(0,nall,success);
    this->zero_timers();
    return NULL;
  }

  this->hd_balancer.balance(cpu_time);
  int inum=this->hd_balancer.get_gpu_count(ago,inum_full);
  this->ans->inum(inum);
  host_start=inum;

  // Build neighbor list on GPU if necessary
  if (ago==0) {
    this->build_nbor_list(inum, inum_full-inum, nall, host_x, host_type,
                    sublo, subhi, tag, nspecial, special, success);
    if (!success)
      return NULL;
    this->atom->cast_q_data(host_q);
    this->hd_balancer.start_timer();
  } else {
    this->atom->cast_x_data(host_x,host_type);
    this->atom->cast_q_data(host_q);
    this->hd_balancer.start_timer();
    this->atom->add_x_data(host_x,host_type);
  }
  this->atom->add_q_data();
  *ilist=this->nbor->host_ilist.begin();
  *jnum=this->nbor->host_acc.begin();

  copy_relations_data(nall, tag, map_array, map_size, sametag, max_same, ago);

  this->device->precompute(ago,inum_full,nall,host_x,host_type,success,host_q,
                     boxlo, prd);

  t_ago = ago;
  loop(eflag,vflag);
  this->ans->copy_answers(eflag,vflag,eatom,vatom);
  this->device->add_ans_object(this->ans);
  this->hd_balancer.stop_timer();

  return this->nbor->host_jlist.begin()-host_start;
}




template class LJ_TIP4PLong<PRECISION,ACC_PRECISION>;
