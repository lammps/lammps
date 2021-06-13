/***************************************************************************
                              yukawa_colloid.cpp
                             -------------------
                            Trung Dac Nguyen (ORNL)

  Class for acceleration of the yukawa/colloid pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : nguyentd@ornl.gov
 ***************************************************************************/

#ifdef USE_OPENCL
#include "yukawa_colloid_cl.h"
#elif defined(USE_CUDART)
const char *yukawa_colloid=0;
#else
#include "yukawa_colloid_cubin.h"
#endif

#include "lal_yukawa_colloid.h"
#include <cassert>
namespace LAMMPS_AL {
#define YukawaColloidT YukawaColloid<numtyp, acctyp>

extern Device<PRECISION,ACC_PRECISION> device;

template <class numtyp, class acctyp>
YukawaColloidT::YukawaColloid() : BaseAtomic<numtyp,acctyp>(),
_max_rad_size(0), _allocated(false) {
}

template <class numtyp, class acctyp>
YukawaColloidT::~YukawaColloid() {
  clear();
}

template <class numtyp, class acctyp>
int YukawaColloidT::bytes_per_atom(const int max_nbors) const {
  return this->bytes_per_atom_atomic(max_nbors);
}

template <class numtyp, class acctyp>
int YukawaColloidT::init(const int ntypes,
                   double **host_cutsq, double **host_a,
                   double **host_offset, double *host_special_lj, const int nlocal,
                   const int nall, const int max_nbors,
                   const int maxspecial, const double cell_size,
                   const double gpu_split, FILE *_screen, const double kappa) {
  int success;
  success=this->init_atomic(nlocal,nall,max_nbors,maxspecial,cell_size,gpu_split,
                            _screen,yukawa_colloid,"k_yukawa_colloid");
  if (success!=0)
    return success;

  if (this->ucl_device->shared_memory() && sizeof(numtyp)==sizeof(double))
    _shared_view=true;
  else
    _shared_view=false;

  // allocate rad

  int ef_nall=nall;
  if (ef_nall==0)
    ef_nall=2000;

  _max_rad_size=static_cast<int>(static_cast<double>(ef_nall)*1.10);

  if (_shared_view==false)
    c_rad.alloc(_max_rad_size,*(this->ucl_device),UCL_WRITE_ONLY,UCL_READ_ONLY);

  rad_tex.get_texture(*(this->pair_program),"rad_tex");
  rad_tex.bind_float(c_rad,1);

  // If atom type constants fit in shared memory use fast kernel
  int lj_types=ntypes;
  shared_types=false;
  int max_shared_types=this->device->max_shared_types();
  if (lj_types<=max_shared_types && this->_block_size>=max_shared_types) {
    lj_types=max_shared_types;
    shared_types=true;
  }
  _lj_types=lj_types;

  _kappa = kappa;

  // Allocate a host write buffer for data initialization
  UCL_H_Vec<numtyp> host_write(lj_types*lj_types*32,*(this->ucl_device),
                               UCL_WRITE_ONLY);

  for (int i=0; i<lj_types*lj_types*32; i++)
    host_write[i]=(numtyp)0.0;

  coeff.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack4(ntypes,lj_types,coeff,host_write,host_a,
                         host_offset,host_cutsq);

  UCL_H_Vec<double> dview;
  sp_lj.alloc(4,*(this->ucl_device),UCL_READ_ONLY);
  dview.view(host_special_lj,4,*(this->ucl_device));
  ucl_copy(sp_lj,dview,false);

  _allocated=true;
  this->_max_bytes=coeff.row_bytes()+sp_lj.row_bytes();
  return 0;
}

template <class numtyp, class acctyp>
void YukawaColloidT::clear() {
  if (!_allocated)
    return;
  _allocated=false;

  coeff.clear();
  sp_lj.clear();

  c_rad.clear();

  this->clear_atomic();
}

template <class numtyp, class acctyp>
double YukawaColloidT::host_memory_usage() const {
  return this->host_memory_usage_atomic()+sizeof(YukawaColloid<numtyp,acctyp>);
}

// ---------------------------------------------------------------------------
// Copy nbor list from host if necessary and then compute atom energies/forces
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void YukawaColloidT::compute(const int f_ago, const int inum_full,
               const int nall, double **host_x, int *host_type, int *ilist,
               int *numj, int **firstneigh, const bool eflag_in,
               const bool vflag_in, const bool eatom, const bool vatom,
               int &host_start, const double cpu_time, bool &success,
               double *rad) {
  this->acc_timers();
  int eflag, vflag;
  if (eatom) eflag=2;
  else if (eflag_in) eflag=1;
  else eflag=0;
  if (vatom) vflag=2;
  else if (vflag_in) vflag=1;
  else vflag=0;

  #ifdef LAL_NO_BLOCK_REDUCE
  if (eflag) eflag=2;
  if (vflag) vflag=2;
  #endif

  this->set_kernel(eflag,vflag);

  // ------------------- Resize rad array --------------------------

  if (nall>_max_rad_size) {
    _max_rad_size=static_cast<int>(static_cast<double>(nall)*1.10);
    if (_shared_view==false) {
      c_rad.resize(_max_rad_size);
      rad_tex.bind_float(c_rad,1);
    }
  }

  // ----------------------------------------------------------------

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

  // -----------------------------------------------------------------

  if (ago==0) {
    this->reset_nbors(nall, inum, ilist, numj, firstneigh, success);
    if (!success)
      return;
  }

  this->atom->cast_x_data(host_x,host_type);
  this->cast_rad_data(rad);
  this->hd_balancer.start_timer();
  this->atom->add_x_data(host_x,host_type);
  this->add_rad_data();

  const int red_blocks=this->loop(eflag,vflag);
  this->ans->copy_answers(eflag_in,vflag_in,eatom,vatom,ilist,red_blocks);
  this->device->add_ans_object(this->ans);
  this->hd_balancer.stop_timer();
}

// ---------------------------------------------------------------------------
// Reneighbor on GPU and then compute per-atom densities
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int** YukawaColloidT::compute(const int ago, const int inum_full,
                const int nall, double **host_x, int *host_type, double *sublo,
                double *subhi, tagint *tag, int **nspecial,
                tagint **special, const bool eflag_in, const bool vflag_in,
                const bool eatom, const bool vatom, int &host_start,
                int **ilist, int **jnum, const double cpu_time, bool &success,
                double *rad) {
  this->acc_timers();
  int eflag, vflag;
  if (eatom) eflag=2;
  else if (eflag_in) eflag=1;
  else eflag=0;
  if (vatom) vflag=2;
  else if (vflag_in) vflag=1;
  else vflag=0;

  #ifdef LAL_NO_BLOCK_REDUCE
  if (eflag) eflag=2;
  if (vflag) vflag=2;
  #endif

  this->set_kernel(eflag,vflag);

  // ------------------- Resize rad array ----------------------------

  if (nall>_max_rad_size) {
    _max_rad_size=static_cast<int>(static_cast<double>(nall)*1.10);
    if (_shared_view==false) {
      c_rad.resize(_max_rad_size);
      rad_tex.bind_float(c_rad,1);
    }
  }

  // -----------------------------------------------------------------

  if (inum_full==0) {
    host_start=0;
    // Make sure textures are correct if realloc by a different hybrid style
    this->resize_atom(0,nall,success);
    this->zero_timers();
    return nullptr;
  }

  // load balance, returning the atom count on the device (inum)
  this->hd_balancer.balance(cpu_time);
  int inum=this->hd_balancer.get_gpu_count(ago,inum_full);
  this->ans->inum(inum);
  host_start=inum;

  // Build neighbor list on GPU if necessary
  if (ago==0) {
    this->build_nbor_list(inum, inum_full-inum, nall, host_x, host_type,
                          sublo, subhi, tag, nspecial, special, success);
    if (!success)
      return nullptr;
    this->cast_rad_data(rad);
    this->hd_balancer.start_timer();
  } else {
    this->atom->cast_x_data(host_x,host_type);
    this->cast_rad_data(rad);
    this->hd_balancer.start_timer();
    this->atom->add_x_data(host_x,host_type);
  }
  this->add_rad_data();
  *ilist=this->nbor->host_ilist.begin();
  *jnum=this->nbor->host_acc.begin();

  const int red_blocks=this->loop(eflag,vflag);
  this->ans->copy_answers(eflag_in,vflag_in,eatom,vatom,red_blocks);
  this->device->add_ans_object(this->ans);
  this->hd_balancer.stop_timer();

  return this->nbor->host_jlist.begin()-host_start;
}

// ---------------------------------------------------------------------------
// Calculate per-atom energies and forces
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int YukawaColloidT::loop(const int eflag, const int vflag) {
  // Compute the block size and grid size to keep all cores busy
  const int BX=this->block_size();
  int GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/
                               (BX/this->_threads_per_atom)));

  int ainum=this->ans->inum();
  int nbor_pitch=this->nbor->nbor_pitch();
  this->time_pair.start();
  if (shared_types) {
    this->k_pair_sel->set_size(GX,BX);
    this->k_pair_sel->run(&this->atom->x, &c_rad, &coeff, &sp_lj,
                          &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                          &this->ans->force, &this->ans->engv, &eflag, &vflag,
                          &ainum, &nbor_pitch, &this->_threads_per_atom, &_kappa);
  } else {
    this->k_pair.set_size(GX,BX);
    this->k_pair.run(&this->atom->x, &c_rad, &coeff, &_lj_types, &sp_lj,
                     &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                     &this->ans->force, &this->ans->engv, &eflag, &vflag,
                     &ainum, &nbor_pitch, &this->_threads_per_atom, &_kappa);
  }
  this->time_pair.stop();
  return GX;
}

template class YukawaColloid<PRECISION,ACC_PRECISION>;
}
