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
int AmoebaT::init(const int ntypes, const int max_amtype, const int max_amclass,
                  const double *host_pdamp, const double *host_thole,
                  const double *host_dirdamp, const int *host_amtype2class,
                  const double *host_special_hal,
                  const double * /*host_special_repel*/,
                  const double * /*host_special_disp*/,
                  const double *host_special_mpole,
                  const double * /*host_special_polar_wscale*/,
                  const double *host_special_polar_piscale,
                  const double *host_special_polar_pscale,
                  const double *host_csix, const double *host_adisp,
                  const int nlocal, const int nall, const int max_nbors,
                  const int maxspecial, const int maxspecial15,
                  const double cell_size, const double gpu_split, FILE *_screen,
                  const double polar_dscale, const double polar_uscale) {
  int success;
  success=this->init_atomic(nlocal,nall,max_nbors,maxspecial,maxspecial15,
                            cell_size,gpu_split,_screen,amoeba,
                            "k_amoeba_multipole", "k_amoeba_udirect2b",
                            "k_amoeba_umutual2b", "k_amoeba_polar",
                            "k_amoeba_fphi_uind", "k_amoeba_fphi_mpole",
                            "k_amoeba_short_nbor", "k_amoeba_special15");
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
    host_write[i].w = host_amtype2class[i];
  }

  coeff_amtype.alloc(max_amtype,*(this->ucl_device), UCL_READ_ONLY);
  ucl_copy(coeff_amtype,host_write,false);

  UCL_H_Vec<numtyp4> host_write2(max_amclass, *(this->ucl_device), UCL_WRITE_ONLY);
  for (int i = 0; i < max_amclass; i++) {
    host_write2[i].x = host_csix[i];
    host_write2[i].y = host_adisp[i];
    host_write2[i].z = (numtyp)0;
    host_write2[i].w = (numtyp)0;
  }

  coeff_amclass.alloc(max_amclass,*(this->ucl_device), UCL_READ_ONLY);
  ucl_copy(coeff_amclass,host_write2,false);

  UCL_H_Vec<numtyp4> dview(5, *(this->ucl_device), UCL_WRITE_ONLY);
  sp_amoeba.alloc(5,*(this->ucl_device),UCL_READ_ONLY);
  for (int i=0; i<5; i++) {
    dview[i].x=host_special_hal[i];
    dview[i].y=host_special_polar_piscale[i];
    dview[i].z=host_special_polar_pscale[i];
    dview[i].w=host_special_mpole[i];
  }
  ucl_copy(sp_amoeba,dview,5,false);

  _polar_dscale = polar_dscale;
  _polar_uscale = polar_uscale;

  _allocated=true;
  this->_max_bytes=coeff_amtype.row_bytes() + coeff_amclass.row_bytes()
    + sp_amoeba.row_bytes() + this->_tep.row_bytes()
    + this->_fieldp.row_bytes() + this->_thetai1.row_bytes()
    + this->_thetai2.row_bytes()  + this->_thetai3.row_bytes()
    + this->_igrid.row_bytes() + this->_cgrid_brick.row_bytes();
  return 0;
}

template <class numtyp, class acctyp>
void AmoebaT::clear() {
  if (!_allocated)
    return;
  _allocated=false;

  coeff_amtype.clear();
  coeff_amclass.clear();
  sp_amoeba.clear();

  this->clear_atomic();
}

template <class numtyp, class acctyp>
double AmoebaT::host_memory_usage() const {
  return this->host_memory_usage_atomic()+sizeof(Amoeba<numtyp,acctyp>);
}

// ---------------------------------------------------------------------------
// Calculate the multipole real-space term, returning tep
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int AmoebaT::multipole_real(const int eflag, const int vflag) {
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
  //   at this point mpole is the first kernel in a time step for AMOEBA

  this->k_short_nbor.set_size(GX,BX);
  this->k_short_nbor.run(&this->atom->x, &this->nbor->dev_nbor,
                         &this->_nbor_data->begin(),
                         &this->dev_short_nbor, &this->_off2_mpole, &ainum,
                         &nbor_pitch, &this->_threads_per_atom);

  this->k_multipole.set_size(GX,BX);
  this->k_multipole.run(&this->atom->x, &this->atom->extra,
                        &coeff_amtype, &sp_amoeba,
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
// Launch the real-space permanent field kernel
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int AmoebaT::udirect2b(const int /*eflag*/, const int /*vflag*/) {
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
  this->k_udirect2b.run(&this->atom->x, &this->atom->extra, &coeff_amtype, &sp_amoeba,
                        &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                        &this->dev_short_nbor,
                        &this->_fieldp, &ainum, &_nall, &nbor_pitch,
                        &this->_threads_per_atom, &this->_aewald, &this->_off2_polar,
                        &_polar_dscale, &_polar_uscale);

  this->time_pair.stop();
  return GX;
}

// ---------------------------------------------------------------------------
// Launch the real-space induced field kernel, returning field and fieldp
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int AmoebaT::umutual2b(const int /*eflag*/, const int /*vflag*/) {
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
  this->k_umutual2b.run(&this->atom->x, &this->atom->extra, &coeff_amtype, &sp_amoeba,
                        &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                        &this->dev_short_nbor, &this->_fieldp, &ainum, &_nall,
                        &nbor_pitch, &this->_threads_per_atom, &this->_aewald,
                        &this->_off2_polar, &_polar_dscale, &_polar_uscale);

  this->time_pair.stop();
  return GX;
}

// ---------------------------------------------------------------------------
// Launch the polar real-space kernel, returning tep
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int AmoebaT::polar_real(const int eflag, const int vflag) {
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
  this->k_polar.run(&this->atom->x, &this->atom->extra, &coeff_amtype, &sp_amoeba,
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

template class Amoeba<PRECISION,ACC_PRECISION>;
}
