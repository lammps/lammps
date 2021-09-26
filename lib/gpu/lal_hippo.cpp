/***************************************************************************
                                 hippo.cpp
                             -------------------
                          Trung Dac Nguyen (Northwestern)

  Class for acceleration of the hippo pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : trung.nguyen@northwestern.edu
 ***************************************************************************/

#if defined(USE_OPENCL)
#include "hippo_cl.h"
#elif defined(USE_CUDART)
const char *hippo=0;
#else
#include "hippo_cubin.h"
#endif

#include "lal_hippo.h"
#include <cassert>
namespace LAMMPS_AL {
#define HippoT Hippo<numtyp, acctyp>

extern Device<PRECISION,ACC_PRECISION> device;

template <class numtyp, class acctyp>
HippoT::Hippo() : BaseAmoeba<numtyp,acctyp>(),
  _allocated(false) {
}

template <class numtyp, class acctyp>
HippoT::~Hippo() {
  clear();
  k_dispersion.clear();
}

template <class numtyp, class acctyp>
int HippoT::bytes_per_atom(const int max_nbors) const {
  return this->bytes_per_atom_atomic(max_nbors);
}

template <class numtyp, class acctyp>
int HippoT::init(const int ntypes, const int max_amtype, const int max_amclass,
                  const double *host_pdamp, const double *host_thole,
                  const double *host_dirdamp, const int *host_amtype2class,
                  const double *host_special_hal,
                  const double *host_special_repel,
                  const double *host_special_disp,
                  const double *host_special_mpole,
                  const double *host_special_polar_wscale,
                  const double *host_special_polar_piscale,
                  const double *host_special_polar_pscale,
                  const double *host_csix, const double *host_adisp,
                  const double *host_pcore, const double *host_palpha,
                  const int nlocal, const int nall, const int max_nbors,
                  const int maxspecial, const int maxspecial15,
                  const double cell_size, const double gpu_split, FILE *_screen,
                  const double polar_dscale, const double polar_uscale) {
  int success;
  success=this->init_atomic(nlocal,nall,max_nbors,maxspecial,maxspecial15,
                            cell_size,gpu_split,_screen,hippo,
                            "k_hippo_multipole",
                            "k_hippo_udirect2b", "k_hippo_umutual2b",
                            "k_hippo_polar", "k_hippo_short_nbor");
  if (success!=0)
    return success;

  // specific to HIPPO
  k_dispersion.set_function(*(this->pair_program),"k_hippo_dispersion");

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
    host_write[i].w = host_amtype2class[i];
  }

  coeff_amtype.alloc(max_amtype,*(this->ucl_device), UCL_READ_ONLY);
  ucl_copy(coeff_amtype,host_write,false);

  UCL_H_Vec<numtyp4> host_write2(max_amclass, *(this->ucl_device), UCL_WRITE_ONLY);
  for (int i = 0; i < max_amclass; i++) {
    host_write2[i].x = host_csix[i];
    host_write2[i].y = host_adisp[i];
    host_write2[i].z = host_pcore[i];
    host_write2[i].w = host_palpha[i];
  }

  coeff_amclass.alloc(max_amclass,*(this->ucl_device), UCL_READ_ONLY);
  ucl_copy(coeff_amclass,host_write2,false);

  UCL_H_Vec<numtyp4> dview(5, *(this->ucl_device), UCL_WRITE_ONLY);
  sp_polar.alloc(5,*(this->ucl_device),UCL_READ_ONLY);
  for (int i=0; i<5; i++) {
    dview[i].x=host_special_polar_wscale[i];
    dview[i].y=host_special_polar_piscale[i];
    dview[i].z=host_special_polar_pscale[i];
    dview[i].w=host_special_mpole[i];
  }
  ucl_copy(sp_polar,dview,5,false);

  sp_nonpolar.alloc(5,*(this->ucl_device),UCL_READ_ONLY);
  for (int i=0; i<5; i++) {
    dview[i].x=host_special_hal[i];
    dview[i].y=host_special_repel[i];
    dview[i].z=host_special_disp[i];
    dview[i].w=(numtyp)0;
  }
  ucl_copy(sp_nonpolar,dview,5,false);

  _polar_dscale = polar_dscale;
  _polar_uscale = polar_uscale;

  _allocated=true;
  this->_max_bytes=coeff_amtype.row_bytes() + coeff_amclass.row_bytes()
    + sp_polar.row_bytes() + sp_nonpolar.row_bytes() + this->_tep.row_bytes();
  return 0;
}

template <class numtyp, class acctyp>
void HippoT::clear() {
  if (!_allocated)
    return;
  _allocated=false;

  coeff_amtype.clear();
  coeff_amclass.clear();
  sp_polar.clear();
  sp_nonpolar.clear();
  
  this->clear_atomic();
}

template <class numtyp, class acctyp>
double HippoT::host_memory_usage() const {
  return this->host_memory_usage_atomic()+sizeof(Hippo<numtyp,acctyp>);
}

// ---------------------------------------------------------------------------
// Prepare for multiple kernel calls in a time step:
//   - reallocate per-atom arrays, if needed
//   - transfer extra data from host to device
//   - build the full neighbor lists for use by different kernels
// ---------------------------------------------------------------------------

template <class numtyp, class acctyp>
int** HippoT::precompute(const int ago, const int inum_full, const int nall,
                              double **host_x, int *host_type, int *host_amtype,
                              int *host_amgroup, double **host_rpole,
                              double **host_uind, double **host_uinp, double *host_pval,
                              double *sublo, double *subhi, tagint *tag,
                              int **nspecial, tagint **special,
                              int *nspecial15, tagint **special15,
                              const bool eflag_in, const bool vflag_in,
                              const bool eatom, const bool vatom, int &host_start,
                              int **&ilist, int **&jnum, const double cpu_time,
                              bool &success, double *host_q, double *boxlo,
                              double *prd) {
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

  // ------------------- Resize 1-5 neighbor arrays ------------------------

  if (nall>this->_nmax) {
    this->_nmax = nall;
    this->dev_nspecial15.clear();
    this->dev_special15.clear();
    this->dev_special15_t.clear();
    this->dev_nspecial15.alloc(nall,*(this->ucl_device),UCL_READ_ONLY);
    this->dev_special15.alloc(this->_maxspecial15*nall,*(this->ucl_device),UCL_READ_ONLY);
    this->dev_special15_t.alloc(nall*this->_maxspecial15,*(this->ucl_device),UCL_READ_ONLY);   
  }

  if (inum_full==0) {
    host_start=0;
    // Make sure textures are correct if realloc by a different hybrid style
    this->resize_atom(0,nall,success);
    this->zero_timers();
    return nullptr;
  }

  this->hd_balancer.balance(cpu_time);
  int inum=this->hd_balancer.get_gpu_count(ago,inum_full);
  this->ans->inum(inum);
  host_start=inum;

  // Build neighbor list on GPU if necessary
  if (ago==0) {
    this->_max_nbors = this->build_nbor_list(inum, inum_full-inum, nall, host_x, host_type,
                    sublo, subhi, tag, nspecial, special, nspecial15, special15,
                    success);
    if (!success)
      return nullptr;
    this->atom->cast_q_data(host_q);
    this->cast_extra_data(host_amtype, host_amgroup, host_rpole, host_uind, host_uinp, host_pval);
    this->hd_balancer.start_timer();
  } else {
    this->atom->cast_x_data(host_x,host_type);
    this->atom->cast_q_data(host_q);
    this->cast_extra_data(host_amtype, host_amgroup, host_rpole, host_uind, host_uinp, host_pval);
    this->hd_balancer.start_timer();
    this->atom->add_x_data(host_x,host_type);
  }
  this->atom->add_q_data();
  this->atom->add_extra_data();

  *ilist=this->nbor->host_ilist.begin();
  *jnum=this->nbor->host_acc.begin();

  this->device->precompute(ago,inum_full,nall,host_x,host_type,success,host_q,
                     boxlo, prd);

  // re-allocate dev_short_nbor if necessary
  if (inum_full*(2+this->_max_nbors) > this->dev_short_nbor.cols()) {
    int _nmax=static_cast<int>(static_cast<double>(inum_full)*1.10);
    this->dev_short_nbor.resize((2+this->_max_nbors)*this->_nmax);
  }

  return this->nbor->host_jlist.begin()-host_start;
}

// ---------------------------------------------------------------------------
// Reneighbor on GPU if necessary, and then compute dispersion real-space
// ---------------------------------------------------------------------------

template <class numtyp, class acctyp>
int** HippoT::compute_dispersion_real(const int ago, const int inum_full,
                                           const int nall, double **host_x,
                                           int *host_type, int *host_amtype,
                                           int *host_amgroup, double **host_rpole,
                                           double *sublo, double *subhi, tagint *tag,
                                           int **nspecial, tagint **special,
                                           int *nspecial15, tagint **special15,
                                           const bool eflag_in, const bool vflag_in,
                                           const bool eatom, const bool vatom,
                                           int &host_start, int **ilist, int **jnum,
                                           const double cpu_time, bool &success,
                                           const double aewald, const double off2_disp,
                                           double *host_q, double *boxlo, double *prd) {
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

  // reallocate per-atom arrays, transfer data from the host
  //   and build the neighbor lists if needed
  // NOTE: 
  //   For now we invoke precompute() again here,
  //     to be able to turn on/off the udirect2b kernel (which comes before this)
  //   Once all the kernels are ready, precompute() is needed only once
  //     in the first kernel in a time step.
  //   We only need to cast uind and uinp from host to device here
  //     if the neighbor lists are rebuilt and other per-atom arrays
  //     (x, type, amtype, amgroup, rpole) are ready on the device.

  int** firstneigh = nullptr;
  firstneigh = precompute(ago, inum_full, nall, host_x, host_type,
                          host_amtype, host_amgroup, host_rpole,
                          nullptr, nullptr, nullptr, sublo, subhi, tag,
                          nspecial, special, nspecial15, special15,
                          eflag_in, vflag_in, eatom, vatom,
                          host_start, ilist, jnum, cpu_time,
                          success, host_q, boxlo, prd);

  this->_off2_disp = off2_disp;
  this->_aewald = aewald;
  const int red_blocks=dispersion_real(eflag,vflag);

  // only copy them back if this is the last kernel
  //   otherwise, commenting out these two lines to leave the answers
  //   (forces, energies and virial) on the device until the last kernel
  //this->ans->copy_answers(eflag_in,vflag_in,eatom,vatom,red_blocks);
  //this->device->add_ans_object(this->ans);

  this->hd_balancer.stop_timer();

  return firstneigh; // nbor->host_jlist.begin()-host_start;
}

// ---------------------------------------------------------------------------
// Calculate the dispersion real-space term, returning tep
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int HippoT::dispersion_real(const int eflag, const int vflag) {
  int ainum=this->ans->inum();
  if (ainum == 0)
    return 0;

  int _nall=this->atom->nall();
  int nbor_pitch=this->nbor->nbor_pitch();

  // Compute the block size and grid size to keep all cores busy
  const int BX=this->block_size();
  int GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/
                               (BX/this->_threads_per_atom)));
  this->time_pair.start();

  // Build the short neighbor list for the cutoff off2_disp,
  //   at this point mpole is the first kernel in a time step
  
  this->k_short_nbor.set_size(GX,BX);
  this->k_short_nbor.run(&this->atom->x, &this->nbor->dev_nbor,
                         &this->_nbor_data->begin(),
                         &this->dev_short_nbor, &this->_off2_disp, &ainum,
                         &nbor_pitch, &this->_threads_per_atom);

  k_dispersion.set_size(GX,BX);
  k_dispersion.run(&this->atom->x, &this->atom->extra,
                         &coeff_amtype, &coeff_amclass, &sp_nonpolar,
                         &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                         &this->dev_short_nbor,
                         &this->ans->force, &this->ans->engv,
                         &eflag, &vflag, &ainum, &_nall, &nbor_pitch,
                         &this->_threads_per_atom,  &this->_aewald,
                         &this->_off2_disp);
  this->time_pair.stop();

  return GX;
}

// ---------------------------------------------------------------------------
// Reneighbor on GPU if necessary, and then compute multipole real-space
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int** HippoT::compute_multipole_real(const int ago, const int inum_full,
                                     const int nall, double **host_x,
                                     int *host_type, int *host_amtype,
                                     int *host_amgroup, double **host_rpole,
                                     double* host_pval, double *sublo,
                                     double *subhi, tagint *tag,
                                     int **nspecial, tagint **special,
                                     int *nspecial15, tagint **special15,
                                     const bool eflag_in, const bool vflag_in,
                                     const bool eatom, const bool vatom,
                                     int &host_start, int **ilist, int **jnum,
                                     const double cpu_time, bool &success,
                                     const double aewald, const double felec,
                                     const double off2_mpole, double *host_q,
                                     double *boxlo, double *prd, void **tep_ptr) {
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

  // reallocate per-atom arrays, transfer data from the host
  //   and build the neighbor lists if needed
  // NOTE: 
  //   For now we invoke precompute() again here,
  //     to be able to turn on/off the udirect2b kernel (which comes before this)
  //   Once all the kernels are ready, precompute() is needed only once
  //     in the first kernel in a time step.
  //   We only need to cast uind and uinp from host to device here
  //     if the neighbor lists are rebuilt and other per-atom arrays
  //     (x, type, amtype, amgroup, rpole) are ready on the device.

  int** firstneigh = nullptr;
  firstneigh = precompute(ago, inum_full, nall, host_x, host_type,
                          host_amtype, host_amgroup, host_rpole,
                          nullptr, nullptr, host_pval, sublo, subhi, tag,
                          nspecial, special, nspecial15, special15,
                          eflag_in, vflag_in, eatom, vatom,
                          host_start, ilist, jnum, cpu_time,
                          success, host_q, boxlo, prd);

  // ------------------- Resize _tep array ------------------------

  if (inum_full>this->_max_tep_size) {
    this->_max_tep_size=static_cast<int>(static_cast<double>(inum_full)*1.10);
    this->_tep.resize(this->_max_tep_size*4);
  }
  *tep_ptr=this->_tep.host.begin();

  this->_off2_mpole = off2_mpole;
  this->_felec = felec;
  this->_aewald = aewald;
  const int red_blocks=multipole_real(eflag,vflag);

  // leave the answers (forces, energies and virial) on the device,
  //   only copy them back in the last kernel (this one, or polar_real once done)
  this->ans->copy_answers(eflag_in,vflag_in,eatom,vatom,red_blocks);
  this->device->add_ans_object(this->ans);

  this->hd_balancer.stop_timer();

  // copy tep from device to host

  this->_tep.update_host(this->_max_tep_size*4,false);
/*
  printf("GPU lib: tep size = %d: max tep size = %d\n", this->_tep.cols(), _max_tep_size);
  for (int i = 0; i < 10; i++) {
    numtyp4* p = (numtyp4*)(&this->_tep[4*i]);
    printf("i = %d; tep = %f %f %f\n", i, p->x, p->y, p->z);
  }
*/  
  return firstneigh; // nbor->host_jlist.begin()-host_start;
}

// ---------------------------------------------------------------------------
// Calculate the multipole real-space term, returning tep
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int HippoT::multipole_real(const int eflag, const int vflag) {
  int ainum=this->ans->inum();
  if (ainum == 0)
    return 0;

  int _nall=this->atom->nall();
  int nbor_pitch=this->nbor->nbor_pitch();

  // Compute the block size and grid size to keep all cores busy
  const int BX=this->block_size();
  int GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/
                               (BX/this->_threads_per_atom)));
  this->time_pair.start();

  // Build the short neighbor list for the cutoff off2_mpole,
  //   at this point mpole is the first kernel in a time step
  
  this->k_short_nbor.set_size(GX,BX);
  this->k_short_nbor.run(&this->atom->x, &this->nbor->dev_nbor,
                         &this->_nbor_data->begin(),
                         &this->dev_short_nbor, &this->_off2_mpole, &ainum,
                         &nbor_pitch, &this->_threads_per_atom);

  this->k_multipole.set_size(GX,BX);
  this->k_multipole.run(&this->atom->x, &this->atom->extra,
                        &coeff_amtype, &coeff_amclass, &sp_polar,
                        &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                        &this->dev_short_nbor,
                        &this->ans->force, &this->ans->engv, &this->_tep,
                        &eflag, &vflag, &ainum, &_nall, &nbor_pitch,
                        &this->_threads_per_atom,  &this->_aewald, &this->_felec,
                        &this->_off2_mpole, &_polar_dscale, &_polar_uscale);
  this->time_pair.stop();

  return GX;
}

// ---------------------------------------------------------------------------
// Calculate the real-space permanent field, returning field and fieldp
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int HippoT::udirect2b(const int eflag, const int vflag) {
  int ainum=this->ans->inum(); 
  if (ainum == 0)
    return 0;

  int _nall=this->atom->nall();
  int nbor_pitch=this->nbor->nbor_pitch();

  // Compute the block size and grid size to keep all cores busy
  const int BX=this->block_size();
  int GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/
                               (BX/this->_threads_per_atom)));
  this->time_pair.start();

  // Build the short neighbor list if not done yet
  if (!this->short_nbor_polar_avail) {
    this->k_short_nbor.set_size(GX,BX);
    this->k_short_nbor.run(&this->atom->x, &this->nbor->dev_nbor,
                           &this->_nbor_data->begin(),
                           &this->dev_short_nbor, &this->_off2_polar, &ainum,
                           &nbor_pitch, &this->_threads_per_atom);
    this->short_nbor_polar_avail = true;
  }
  
  this->k_udirect2b.set_size(GX,BX);
  this->k_udirect2b.run(&this->atom->x, &this->atom->extra, &coeff_amtype, &sp_polar,
                        &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                        &this->dev_short_nbor,
                        &this->_fieldp, &ainum, &_nall, &nbor_pitch,
                        &this->_threads_per_atom, &this->_aewald, &this->_off2_polar,
                        &_polar_dscale, &_polar_uscale);

  this->time_pair.stop();
  return GX;
}

// ---------------------------------------------------------------------------
// Calculate the real-space induced field, returning field and fieldp
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int HippoT::umutual2b(const int eflag, const int vflag) {
  int ainum=this->ans->inum();
  if (ainum == 0)
    return 0;

  int _nall=this->atom->nall();
  int nbor_pitch=this->nbor->nbor_pitch();

  // Compute the block size and grid size to keep all cores busy
  const int BX=this->block_size();
  int GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/
                               (BX/this->_threads_per_atom)));
  this->time_pair.start();

  // Build the short neighbor list if not done yet
  if (!this->short_nbor_polar_avail) {
    this->k_short_nbor.set_size(GX,BX);
    this->k_short_nbor.run(&this->atom->x, &this->nbor->dev_nbor,
                           &this->_nbor_data->begin(), &this->dev_short_nbor,
                           &this->_off2_polar, &ainum, &nbor_pitch,
                           &this->_threads_per_atom);
    this->short_nbor_polar_avail = true;
  }

  this->k_umutual2b.set_size(GX,BX);
  this->k_umutual2b.run(&this->atom->x, &this->atom->extra, &coeff_amtype, &sp_polar,
                        &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                        &this->dev_short_nbor, &this->_fieldp, &ainum, &_nall,
                        &nbor_pitch, &this->_threads_per_atom, &this->_aewald,
                        &this->_off2_polar, &_polar_dscale, &_polar_uscale);

  this->time_pair.stop();
  return GX;
}

// ---------------------------------------------------------------------------
// Calculate the polar real-space term, returning tep
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int HippoT::polar_real(const int eflag, const int vflag) {
  int ainum=this->ans->inum();
  if (ainum == 0)
    return 0;

  int _nall=this->atom->nall();
  int nbor_pitch=this->nbor->nbor_pitch();

  // Compute the block size and grid size to keep all cores busy
  const int BX=this->block_size();
  int GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/
                               (BX/this->_threads_per_atom)));
  this->time_pair.start();

  // Build the short neighbor list if not done yet
  if (!this->short_nbor_polar_avail) {
    this->k_short_nbor.set_size(GX,BX);
    this->k_short_nbor.run(&this->atom->x, &this->nbor->dev_nbor,
                          &this->_nbor_data->begin(),
                          &this->dev_short_nbor, &this->_off2_polar, &ainum,
                          &nbor_pitch, &this->_threads_per_atom);
    this->short_nbor_polar_avail = true;
  }

  this->k_polar.set_size(GX,BX);
  this->k_polar.run(&this->atom->x, &this->atom->extra, &coeff_amtype, &sp_polar,
                    &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                    &this->dev_short_nbor,
                    &this->ans->force, &this->ans->engv, &this->_tep,
                    &eflag, &vflag, &ainum, &_nall, &nbor_pitch,
                    &this->_threads_per_atom,  &this->_aewald, &this->_felec,
                    &this->_off2_polar, &_polar_dscale, &_polar_uscale);
  this->time_pair.stop();

  // Signal that short nbor list is not avail for the next time step
  //   do it here because polar_real() is the last kernel in a time step at this point

  this->short_nbor_polar_avail = false;

  return GX;
}

template class Hippo<PRECISION,ACC_PRECISION>;
}
