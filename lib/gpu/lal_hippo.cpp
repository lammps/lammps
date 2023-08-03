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
  k_repulsion.clear();
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
                 const double *host_special_repel, const double *host_special_disp,
                 const double *host_special_mpole,
                 const double *host_special_polar_wscale,
                 const double *host_special_polar_piscale,
                 const double *host_special_polar_pscale,
                 const double *host_sizpr, const double *host_dmppr, const double *host_elepr,
                 const double *host_csix, const double *host_adisp,
                 const double *host_pcore, const double *host_palpha,
                 const int nlocal, const int nall, const int max_nbors,
                 const int maxspecial, const int maxspecial15,
                 const double cell_size, const double gpu_split, FILE *_screen,
                 const double polar_dscale, const double polar_uscale) {
  int success;
  success=this->init_atomic(nlocal,nall,max_nbors,maxspecial,maxspecial15,
                            cell_size,gpu_split,_screen,hippo,
                            "k_hippo_multipole", "k_hippo_udirect2b",
                            "k_hippo_umutual2b", "k_hippo_polar",
                            "k_hippo_fphi_uind", "k_hippo_fphi_mpole",
                            "k_hippo_short_nbor", "k_hippo_special15");
  if (success!=0)
    return success;

  // specific to HIPPO
  k_repulsion.set_function(*(this->pair_program),"k_hippo_repulsion");
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

  for (int i = 0; i < max_amtype; i++) {
    host_write[i].x = host_sizpr[i];
    host_write[i].y = host_dmppr[i];
    host_write[i].z = host_elepr[i];
    host_write[i].w = (numtyp)0;
  }

  coeff_rep.alloc(max_amtype,*(this->ucl_device), UCL_READ_ONLY);
  ucl_copy(coeff_rep,host_write,false);

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
    dview[i].x=host_special_repel[i];
    dview[i].y=host_special_disp[i];
    dview[i].z=(numtyp)0;
    dview[i].w=(numtyp)0;
  }
  ucl_copy(sp_nonpolar,dview,5,false);

  _polar_dscale = polar_dscale;
  _polar_uscale = polar_uscale;

  _allocated=true;
  this->_max_bytes=coeff_amtype.row_bytes() + coeff_rep.row_bytes()
    + coeff_amclass.row_bytes() + sp_polar.row_bytes()
    + sp_nonpolar.row_bytes() + this->_tep.row_bytes()
    + this->_fieldp.row_bytes() + this->_thetai1.row_bytes()
    + this->_thetai2.row_bytes()  + this->_thetai3.row_bytes()
    + this->_igrid.row_bytes() + this->_cgrid_brick.row_bytes();
  return 0;
}

template <class numtyp, class acctyp>
void HippoT::clear() {
  if (!_allocated)
    return;
  _allocated=false;

  coeff_amtype.clear();
  coeff_rep.clear();
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
// Compute the repulsion term, returning tep
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void HippoT::compute_repulsion(const int /*ago*/, const int inum_full,
                               const int /*nall*/, double ** /*host_x*/,
                               int * /*host_type*/, int * /*host_amtype*/,
                               int * /*host_amgroup*/, double ** /*host_rpole*/,
                               double * /*sublo*/, double * /*subhi*/, tagint * /*tag*/,
                               int ** /*nspecial*/, tagint ** /*special*/,
                               int * /*nspecial15*/, tagint ** /*special15*/,
                               const bool eflag_in, const bool vflag_in,
                               const bool eatom, const bool vatom,
                               int & /*host_start*/, int ** /*ilist*/, int ** /*jnum*/,
                               const double /*cpu_time*/, bool & /*success*/,
                               const double /*aewald*/, const double off2_repulse,
                               double * /*host_q*/, double * /*boxlo*/, double * /*prd*/,
                               double cut2, double c0, double c1, double c2,
                               double c3, double c4, double c5, void **tep_ptr) {
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

  // ------------------- Resize _tep array ------------------------

  if (inum_full>this->_max_tep_size) {
    this->_max_tep_size=static_cast<int>(static_cast<double>(inum_full)*1.10);
    this->_tep.resize(this->_max_tep_size*3);
  }
  *tep_ptr=this->_tep.host.begin();

  this->_off2_repulse = off2_repulse;
  _cut2 = cut2;
  _c0 = c0;
  _c1 = c1;
  _c2 = c2;
  _c3 = c3;
  _c4 = c4;
  _c5 = c5;
  repulsion(this->_eflag,this->_vflag);

  // copy tep from device to host
  this->_tep.update_host(this->_max_tep_size*3,false);
}

// ---------------------------------------------------------------------------
// Launch the repulsion kernel
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int HippoT::repulsion(const int eflag, const int vflag) {
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
  //   at this point repuslion is the first kernel in a time step for HIPPO

  this->k_short_nbor.set_size(GX,BX);
  this->k_short_nbor.run(&this->atom->x, &this->nbor->dev_nbor,
                         &this->_nbor_data->begin(),
                         &this->dev_short_nbor, &this->_off2_repulse, &ainum,
                         &nbor_pitch, &this->_threads_per_atom);

  k_repulsion.set_size(GX,BX);
  k_repulsion.run(&this->atom->x, &this->atom->extra,
                  &coeff_rep, &sp_nonpolar,
                  &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                  &this->dev_short_nbor,
                  &this->ans->force, &this->ans->engv, &this->_tep,
                  &eflag, &vflag, &ainum, &_nall, &nbor_pitch,
                  &this->_threads_per_atom,  &this->_aewald,
                  &this->_off2_repulse, &_cut2,
                  &_c0, &_c1, &_c2, &_c3, &_c4, &_c5);
  this->time_pair.stop();

  return GX;
}

// ---------------------------------------------------------------------------
// Compute dispersion real-space
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void HippoT::compute_dispersion_real(int *host_amtype, int *host_amgroup,
                                      double **host_rpole, const double aewald,
                                      const double off2_disp) {

  // cast necessary data arrays from host to device

  this->cast_extra_data(host_amtype, host_amgroup, host_rpole,
                        nullptr, nullptr, nullptr);
  this->atom->add_extra_data();

  this->_off2_disp = off2_disp;
  this->_aewald = aewald;
  dispersion_real(this->_eflag,this->_vflag);

  // only copy them back if this is the last kernel
  //   otherwise, commenting out these two lines to leave the answers
  //   (forces, energies and virial) on the device until the last kernel
  //this->ans->copy_answers(eflag_in,vflag_in,eatom,vatom,red_blocks);
  //this->device->add_ans_object(this->ans);
}

// ---------------------------------------------------------------------------
// Launch the dispersion real-space kernel
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
  //   at this point dispersion is the first kernel in a time step

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
// Compute the multipole real-space term, returning tep
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void HippoT::compute_multipole_real(const int /*ago*/, const int inum_full,
                                    const int /*nall*/, double ** /*host_x*/,
                                    int * /*host_type*/, int * /*host_amtype*/,
                                    int * /*host_amgroup*/, double ** /*host_rpole*/,
                                    double* host_pval, double * /*sublo*/,
                                    double * /*subhi*/, tagint * /*tag*/,
                                    int ** /*nspecial*/, tagint ** /*special*/,
                                    int * /*nspecial15*/, tagint ** /*special15*/,
                                    const bool /*eflag_in*/, const bool /*vflag_in*/,
                                    const bool /*eatom*/, const bool /*vatom*/,
                                    int & /*host_start*/, int ** /*ilist*/, int ** /*jnum*/,
                                    const double /*cpu_time*/, bool & /*success*/,
                                     const double aewald, const double felec,
                                    const double off2_mpole, double * /*host_q*/,
                                    double * /*boxlo*/, double * /*prd*/, void **tep_ptr) {

  // cast necessary data arrays from host to device

  this->cast_extra_data(nullptr, nullptr, nullptr, nullptr, nullptr, host_pval);
  this->atom->add_extra_data();

  // ------------------- Resize _tep array ------------------------

  if (inum_full>this->_max_tep_size) {
    this->_max_tep_size=static_cast<int>(static_cast<double>(inum_full)*1.10);
    this->_tep.resize(this->_max_tep_size*3);
  }
  *tep_ptr=this->_tep.host.begin();

  this->_off2_mpole = off2_mpole;
  this->_felec = felec;
  this->_aewald = aewald;
  multipole_real(this->_eflag,this->_vflag);

  // copy tep from device to host
  this->_tep.update_host(this->_max_tep_size*3,false);
}

// ---------------------------------------------------------------------------
// Launch the multipole real-space kernel
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

  // Build the short neighbor list for the cutoff off2_mpole

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
// Compute the direct real space part of the permanent field
//   returning field and fieldp
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void HippoT::compute_udirect2b(int * /*host_amtype*/, int * /*host_amgroup*/, double ** /*host_rpole*/,
                                double **host_uind, double **host_uinp, double* host_pval,
                                const double aewald, const double off2_polar,
                                void** fieldp_ptr) {

  // all the necessary data arrays are already copied from host to device

  this->cast_extra_data(nullptr, nullptr, nullptr, host_uind, host_uinp, host_pval);
  this->atom->add_extra_data();

  if (this->_max_tep_size>this->_max_fieldp_size) {
    this->_max_fieldp_size = this->_max_tep_size;
    this->_fieldp.resize(this->_max_fieldp_size*6);
  }
  *fieldp_ptr=this->_fieldp.host.begin();

  this->_off2_polar = off2_polar;
  this->_aewald = aewald;
  udirect2b(this->_eflag,this->_vflag);

  // copy field and fieldp from device to host (_fieldp store both arrays, one after another)

  this->_fieldp.update_host(this->_max_fieldp_size*6,false);
}

// ---------------------------------------------------------------------------
// Launch the real-space permanent field kernel
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int HippoT::udirect2b(const int /*eflag*/, const int /*vflag*/) {
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

  // Build the short neighbor list for the cutoff _off2_polar, if not done yet
  //   this is the first kernel in a time step where _off2_polar is used

  if (!this->short_nbor_polar_avail) {
    this->k_short_nbor.set_size(GX,BX);
    this->k_short_nbor.run(&this->atom->x, &this->nbor->dev_nbor,
                           &this->_nbor_data->begin(),
                           &this->dev_short_nbor, &this->_off2_polar, &ainum,
                           &nbor_pitch, &this->_threads_per_atom);
    this->short_nbor_polar_avail = true;
  }

  this->k_udirect2b.set_size(GX,BX);
  this->k_udirect2b.run(&this->atom->x, &this->atom->extra,
                        &coeff_amtype, &coeff_amclass, &sp_polar,
                        &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                        &this->dev_short_nbor,
                        &this->_fieldp, &ainum, &_nall, &nbor_pitch,
                        &this->_threads_per_atom, &this->_aewald, &this->_off2_polar,
                        &_polar_dscale, &_polar_uscale);

  this->time_pair.stop();
  return GX;
}

// ---------------------------------------------------------------------------
// Compute the direct real space term of the induced field
//   returning field and fieldp
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void HippoT::compute_umutual2b(int * /*host_amtype*/, int * /*host_amgroup*/, double ** /*host_rpole*/,
                               double **host_uind, double **host_uinp, double * /*host_pval*/,
                               const double aewald, const double off2_polar, void ** /*fieldp_ptr*/) {

  // cast necessary data arrays from host to device

  this->cast_extra_data(nullptr, nullptr, nullptr, host_uind, host_uinp, nullptr);
  this->atom->add_extra_data();

  this->_off2_polar = off2_polar;
  this->_aewald = aewald;
  umutual2b(this->_eflag,this->_vflag);

  // copy field and fieldp from device to host (_fieldp store both arrays, one after another)
  // NOTE: move this step to update_fieldp() to delay device-host transfer
  // *fieldp_ptr=this->_fieldp.host.begin();
  // this->_fieldp.update_host(this->_max_fieldp_size*8,false);
}

// ---------------------------------------------------------------------------
// Launch the real-space induced field kernel
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int HippoT::umutual2b(const int /*eflag*/, const int /*vflag*/) {
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
  this->k_umutual2b.run(&this->atom->x, &this->atom->extra,
                        &coeff_amtype, &coeff_amclass, &sp_polar,
                        &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                        &this->dev_short_nbor, &this->_fieldp, &ainum, &_nall,
                        &nbor_pitch, &this->_threads_per_atom, &this->_aewald,
                        &this->_off2_polar, &_polar_dscale, &_polar_uscale);

  this->time_pair.stop();
  return GX;
}

// ---------------------------------------------------------------------------
// Reneighbor on GPU if necessary, and then compute polar real-space
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void HippoT::compute_polar_real(int * /*host_amtype*/, int * /*host_amgroup*/, double ** /*host_rpole*/,
                                double **host_uind, double **host_uinp, double * /*host_pval*/,
                                const bool eflag_in, const bool vflag_in,
                                const bool eatom, const bool vatom,
                                const double aewald, const double felec,
                                const double off2_polar, void **tep_ptr) {
  // cast necessary data arrays from host to device

  this->cast_extra_data(nullptr, nullptr, nullptr, host_uind, host_uinp, nullptr);
  this->atom->add_extra_data();

  *tep_ptr=this->_tep.host.begin();

  this->_off2_polar = off2_polar;
  this->_felec = felec;
  this->_aewald = aewald;
  const int red_blocks=polar_real(this->_eflag,this->_vflag);

  // only copy answers (forces, energies and virial) back from the device
  //   in the last kernel in a timestep (which is polar_real here)
  this->ans->copy_answers(eflag_in,vflag_in,eatom,vatom,red_blocks);
  this->device->add_ans_object(this->ans);

  // copy tep from device to host
  this->_tep.update_host(this->_max_tep_size*3,false);
}

// ---------------------------------------------------------------------------
// Launch the polar real-space kernel
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
  const int GX=static_cast<int>(ceil(static_cast<double>(ainum)/(BX/this->_threads_per_atom)));
  /*
  const int cus = this->device->gpu->cus();
  while (GX < cus && GX > 1) {
    BX /= 2;
    GX=static_cast<int>(ceil(static_cast<double>(ainum)/(BX/this->_threads_per_atom)));
  }
  */
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
  this->k_polar.run(&this->atom->x, &this->atom->extra,
                    &coeff_amtype, &coeff_amclass, &sp_polar,
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
