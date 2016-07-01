/***************************************************************************
                                lal_table.cpp
                             -------------------
                           Trung Dac Nguyen (ORNL)

  Class for acceleration of the table pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : nguyentd@ornl.gov
 ***************************************************************************/

#if defined(USE_OPENCL)
#include "table_cl.h"
#elif defined(USE_CUDART)
const char *table=0;
#else
#include "table_cubin.h"
#endif

#include "lal_table.h"
#include <cassert>
using namespace LAMMPS_AL;
#define TableT Table<numtyp, acctyp>

#define LOOKUP 0
#define LINEAR 1
#define SPLINE 2
#define BITMAP 3

extern Device<PRECISION,ACC_PRECISION> device;

template <class numtyp, class acctyp>
TableT::Table() : BaseAtomic<numtyp,acctyp>(),
  _allocated(false), _compiled_styles(false) {
}

template <class numtyp, class acctyp>
TableT::~Table() {
  clear();
}

template <class numtyp, class acctyp>
int TableT::bytes_per_atom(const int max_nbors) const {
  return this->bytes_per_atom_atomic(max_nbors);
}

template <class numtyp, class acctyp>
int TableT::init(const int ntypes,
                double **host_cutsq, double ***host_table_coeffs,
                double **host_table_data,
                double *host_special_lj, const int nlocal,
                const int nall, const int max_nbors,
                const int maxspecial, const double cell_size,
                const double gpu_split, FILE *_screen,
                int tabstyle, int ntables, int tablength) {
  int success;
  success=this->init_atomic(nlocal,nall,max_nbors,maxspecial,cell_size,
                            gpu_split,_screen,table,"k_table");
  if (success!=0)
    return success;

  k_pair_linear.set_function(*(this->pair_program),"k_table_linear");
  k_pair_linear_fast.set_function(*(this->pair_program),"k_table_linear_fast");
  k_pair_spline.set_function(*(this->pair_program),"k_table_spline");
  k_pair_spline_fast.set_function(*(this->pair_program),"k_table_spline_fast");
  k_pair_bitmap.set_function(*(this->pair_program),"k_table_bitmap");
  k_pair_bitmap_fast.set_function(*(this->pair_program),"k_table_bitmap_fast");
  _compiled_styles = true;

  // If atom type constants fit in shared memory use fast kernel
  int lj_types=ntypes;
  shared_types=false;
  int max_shared_types=this->device->max_shared_types();
  if (lj_types<=max_shared_types && this->_block_size>=max_shared_types) {
    lj_types=max_shared_types;
    shared_types=true;
  }
  _lj_types=lj_types;

  _tabstyle = tabstyle;
  _ntables = ntables;
  if (tabstyle != BITMAP) _tablength = tablength;
  else _tablength = 1 << tablength;

  // Allocate a host write buffer for data initialization
  UCL_H_Vec<int> host_write_int(lj_types*lj_types,*(this->ucl_device),
                               UCL_WRITE_ONLY);

  for (int i=0; i<lj_types*lj_types; i++)
    host_write_int[i] = 0;

  tabindex.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  nshiftbits.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  nmask.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);

  for (int ix=1; ix<ntypes; ix++)
    for (int iy=1; iy<ntypes; iy++)
      host_write_int[ix*lj_types+iy] = (int)host_table_coeffs[ix][iy][0]; // tabindex
  ucl_copy(tabindex,host_write_int,false);

  for (int ix=1; ix<ntypes; ix++)
    for (int iy=1; iy<ntypes; iy++)
      host_write_int[ix*lj_types+iy] = (int)host_table_coeffs[ix][iy][1]; // nshiftbits
  ucl_copy(nshiftbits,host_write_int,false);

  for (int ix=1; ix<ntypes; ix++)
    for (int iy=1; iy<ntypes; iy++)
      host_write_int[ix*lj_types+iy] = (int)host_table_coeffs[ix][iy][2]; // nmask
  ucl_copy(nmask,host_write_int,false);

  UCL_H_Vec<numtyp4> host_write(lj_types*lj_types,*(this->ucl_device),
                               UCL_WRITE_ONLY);

  coeff2.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  for (int ix=1; ix<ntypes; ix++)
    for (int iy=1; iy<ntypes; iy++) {
      host_write[ix*lj_types+iy].x = host_table_coeffs[ix][iy][3]; // innersq
      host_write[ix*lj_types+iy].y = host_table_coeffs[ix][iy][4]; // invdelta
      host_write[ix*lj_types+iy].z = host_table_coeffs[ix][iy][5]; // deltasq6
      host_write[ix*lj_types+iy].w = (numtyp)0.0;
  }
  ucl_copy(coeff2,host_write,false);

  // Allocate tablength arrays
  UCL_H_Vec<numtyp4> host_write2(_ntables*_tablength,*(this->ucl_device),
                                 UCL_WRITE_ONLY);
  for (int i=0; i<_ntables*_tablength; i++) {
    host_write2[i].x = 0.0;
    host_write2[i].y = 0.0;
    host_write2[i].z = 0.0;
    host_write2[i].w = 0.0;
  }

  coeff3.alloc(_ntables*_tablength,*(this->ucl_device),UCL_READ_ONLY);
  for (int n=0; n<_ntables; n++) {
    if (tabstyle == LOOKUP) {
      for (int k=0; k<_tablength-1; k++) {
          host_write2[n*_tablength+k].x = (numtyp)0;
          host_write2[n*_tablength+k].y = host_table_data[n][6*k+1]; // e
          host_write2[n*_tablength+k].z = host_table_data[n][6*k+2]; // f
          host_write2[n*_tablength+k].w = (numtyp)0;
      }
    } else if (tabstyle == LINEAR || tabstyle == SPLINE || tabstyle == BITMAP) {
      for (int k=0; k<_tablength; k++) {
          host_write2[n*_tablength+k].x = host_table_data[n][6*k+0]; // rsq
          host_write2[n*_tablength+k].y = host_table_data[n][6*k+1]; // e
          host_write2[n*_tablength+k].z = host_table_data[n][6*k+2]; // f
          host_write2[n*_tablength+k].w = (numtyp)0;
      }
    }
  }
  ucl_copy(coeff3,host_write2,false);

  coeff4.alloc(_ntables*_tablength,*(this->ucl_device),UCL_READ_ONLY);
  for (int i=0; i<_ntables*_tablength; i++) {
    host_write2[i].x = 0.0;
    host_write2[i].y = 0.0;
    host_write2[i].z = 0.0;
    host_write2[i].w = 0.0;
  }

  for (int n=0; n<_ntables; n++) {
    if (tabstyle == LINEAR) {
      for (int k=0; k<_tablength-1; k++) {
        host_write2[n*_tablength+k].x = (numtyp)0;
        host_write2[n*_tablength+k].y = host_table_data[n][6*k+3]; // de
        host_write2[n*_tablength+k].z = host_table_data[n][6*k+4]; // df
        host_write2[n*_tablength+k].w = (numtyp)0;
      }
    } else if (tabstyle == SPLINE) {
      for (int k=0; k<_tablength; k++) {
        host_write2[n*_tablength+k].x = (numtyp)0;
        host_write2[n*_tablength+k].y = host_table_data[n][6*k+3]; // e2
        host_write2[n*_tablength+k].z = host_table_data[n][6*k+4]; // f2
        host_write2[n*_tablength+k].w = (numtyp)0;
      }
    } else if (tabstyle == BITMAP) {
      for (int k=0; k<_tablength; k++) {
        host_write2[n*_tablength+k].x = (numtyp)0;
        host_write2[n*_tablength+k].y = host_table_data[n][6*k+3]; // de
        host_write2[n*_tablength+k].z = host_table_data[n][6*k+4]; // df
        host_write2[n*_tablength+k].w = host_table_data[n][6*k+5]; // drsq
      }
    }
  }
  ucl_copy(coeff4,host_write2,false);

  UCL_H_Vec<numtyp> host_rsq(lj_types*lj_types,*(this->ucl_device),
                             UCL_WRITE_ONLY);
  cutsq.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack1(ntypes,lj_types,cutsq,host_rsq,host_cutsq);

  UCL_H_Vec<double> dview;
  sp_lj.alloc(4,*(this->ucl_device),UCL_READ_ONLY);
  dview.view(host_special_lj,4,*(this->ucl_device));
  ucl_copy(sp_lj,dview,false);

  _allocated=true;
  this->_max_bytes=tabindex.row_bytes()+nshiftbits.row_bytes()
    +nmask.row_bytes()+coeff2.row_bytes()
    +coeff3.row_bytes()+coeff4.row_bytes()+cutsq.row_bytes()
    +sp_lj.row_bytes();
  return 0;
}

template <class numtyp, class acctyp>
void TableT::clear() {
  if (!_allocated)
    return;
  _allocated=false;

  tabindex.clear();
  nshiftbits.clear();
  nmask.clear();
  coeff2.clear();
  coeff3.clear();
  coeff4.clear();
  sp_lj.clear();

  if (_compiled_styles) {
    k_pair_linear_fast.clear();
    k_pair_linear.clear();
    k_pair_spline_fast.clear();
    k_pair_spline.clear();
    k_pair_bitmap_fast.clear();
    k_pair_bitmap.clear();
    _compiled_styles=false;
  }

  this->clear_atomic();
}

template <class numtyp, class acctyp>
double TableT::host_memory_usage() const {
  return this->host_memory_usage_atomic()+sizeof(Table<numtyp,acctyp>);
}

// ---------------------------------------------------------------------------
// Calculate energies, forces, and torques
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void TableT::loop(const bool _eflag, const bool _vflag) {
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

  int GX=static_cast<int>(ceil(static_cast<double>(this->ans->inum())/
                               (BX/this->_threads_per_atom)));

  int ainum=this->ans->inum();
  int nbor_pitch=this->nbor->nbor_pitch();
  this->time_pair.start();
  if (shared_types) {
    if (_tabstyle == LOOKUP) {
      this->k_pair_fast.set_size(GX,BX);
      this->k_pair_fast.run(&this->atom->x, &tabindex, &coeff2, &coeff3,
                            &coeff4, &cutsq, &sp_lj, &this->nbor->dev_nbor,
                            &this->_nbor_data->begin(), &this->ans->force,
                            &this->ans->engv, &eflag, &vflag, &ainum,
                            &nbor_pitch, &this->_threads_per_atom, &_tablength);
    } else if (_tabstyle == LINEAR) {
      this->k_pair_linear_fast.set_size(GX,BX);
      this->k_pair_linear_fast.run(&this->atom->x, &tabindex, &coeff2,
                                   &coeff3, &coeff4, &cutsq, &sp_lj,
                                   &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                                   &this->ans->force, &this->ans->engv,
                                   &eflag, &vflag, &ainum, &nbor_pitch,
                                   &this->_threads_per_atom, &_tablength);
    } else if (_tabstyle == SPLINE) {
      this->k_pair_spline_fast.set_size(GX,BX);
      this->k_pair_spline_fast.run(&this->atom->x, &tabindex, &coeff2,
                                   &coeff3, &coeff4, &cutsq, &sp_lj,
                                   &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                                   &this->ans->force, &this->ans->engv,
                                   &eflag, &vflag, &ainum, &nbor_pitch,
                                   &this->_threads_per_atom, &_tablength);
    } else if (_tabstyle == BITMAP) {
      this->k_pair_bitmap_fast.set_size(GX,BX);
      this->k_pair_bitmap_fast.run(&this->atom->x, &tabindex, &nshiftbits,
                                   &nmask, &coeff2, &coeff3, &coeff4, &cutsq,
                                   &sp_lj, &this->nbor->dev_nbor,
                                   &this->_nbor_data->begin(), &this->ans->force,
                                   &this->ans->engv, &eflag, &vflag,
                                   &ainum, &nbor_pitch,
                                   &this->_threads_per_atom, &_tablength);
    }
  } else {
    if (_tabstyle == LOOKUP) {
      this->k_pair.set_size(GX,BX);
      this->k_pair.run(&this->atom->x, &tabindex, &coeff2, &coeff3,
                       &coeff4, &_lj_types, &cutsq, &sp_lj,
                       &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                       &this->ans->force, &this->ans->engv, &eflag,
                       &vflag, &ainum, &nbor_pitch, &this->_threads_per_atom,
                       &_tablength);
    } else if (_tabstyle == LINEAR) {
      this->k_pair_linear.set_size(GX,BX);
      this->k_pair_linear.run(&this->atom->x, &tabindex, &coeff2, &coeff3,
                              &coeff4, &_lj_types, &cutsq, &sp_lj,
                              &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                              &this->ans->force, &this->ans->engv, &eflag,
                              &vflag, &ainum, &nbor_pitch,
                              &this->_threads_per_atom, &_tablength);
    } else if (_tabstyle == SPLINE) {
      this->k_pair_spline.set_size(GX,BX);
      this->k_pair_spline.run(&this->atom->x, &tabindex, &coeff2, &coeff3,
                              &coeff4, &_lj_types, &cutsq, &sp_lj,
                              &this->nbor->dev_nbor, &this->_nbor_data->begin(),
                              &this->ans->force, &this->ans->engv, &eflag,
                              &vflag, &ainum, &nbor_pitch,
                              &this->_threads_per_atom, &_tablength);
    } else if (_tabstyle == BITMAP) {
      this->k_pair_bitmap.set_size(GX,BX);
      this->k_pair_bitmap.run(&this->atom->x, &tabindex, &nshiftbits,
                              &nmask, &coeff2, &coeff3, &coeff4, &_lj_types,
                              &cutsq, &sp_lj, &this->nbor->dev_nbor,
                              &this->_nbor_data->begin(), &this->ans->force,
                              &this->ans->engv, &eflag, &vflag, &ainum,
                              &nbor_pitch, &this->_threads_per_atom,
                              &_tablength);
    }
  }
  this->time_pair.stop();
}

template class Table<PRECISION,ACC_PRECISION>;
