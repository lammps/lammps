/***************************************************************************
                                lal_eam.cpp
                             -------------------
                      W. Michael Brown, Trung Dac Nguyen (ORNL)

  Class for acceleration of the eam pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : 
    email                : brownw@ornl.gov nguyentd@ornl.gov
 ***************************************************************************/
 
#ifdef USE_OPENCL
#include "eam_cl.h"
#else
#include "eam_ptx.h"
#endif

#include "lal_eam.h"
#include <cassert>
using namespace LAMMPS_AL;
#define EAMT EAM<numtyp, acctyp>

extern Device<PRECISION,ACC_PRECISION> device;

template <class numtyp, class acctyp>
EAMT::EAM() : BaseCharge<numtyp,acctyp>(), 
                      _allocated(false) {
}

template <class numtyp, class acctyp>
EAMT::~EAM() {
  clear();
}
 
template <class numtyp, class acctyp>
int EAMT::bytes_per_atom(const int max_nbors) const {
  return this->bytes_per_atom_atomic(max_nbors);
}

template <class numtyp, class acctyp>
int EAMT::init(const int ntypes, double host_cutforcesq,
              int **host_type2rhor, int **host_type2z2r,
              double ***host_rhor_spline, double ***host_z2r_spline,
              double rdr, int nrhor, int nz2r, int nr,
              const int nlocal, const int nall, const int max_nbors,
              const int maxspecial, const double cell_size,
              const double gpu_split, FILE *_screen) 
{
   
  int success;
  success=this->init_atomic(nlocal,nall,max_nbors,maxspecial,cell_size,gpu_split,
                            _screen,eam);
  
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
  
  _ntypes=lj_types;
  _cutforcesq=host_cutforcesq;
  _rdr=rdr;
  _nrhor=nrhor;
  _nz2r=nz2r;
  _nr=nr;
 
  UCL_H_Vec<numtyp2> dview_type(lj_types*lj_types,*(this->ucl_device),
                               UCL_WRITE_OPTIMIZED);
  
  for (int i=0; i<lj_types*lj_types; i++) {
    dview_type[i].x=(numtyp)0.0; 
    dview_type[i].y=(numtyp)0.0; 
  }
                                
  // pack type2rhor and type2z2r
  type2rhor_z2r.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  
  for (int ix=0; ix<ntypes; ix++)
    for (int iy=0; iy<ntypes; iy++) {
      dview_type[ix*lj_types+iy].x = host_type2rhor[ix][iy];
      dview_type[ix*lj_types+iy].y = host_type2z2r[ix][iy];
  }
  
  ucl_copy(type2rhor_z2r,dview_type,false);
  


  // pack rhor_spline
  UCL_H_Vec<numtyp> dview_rhor_spline(nrhor*(nr+1)*7,*(this->ucl_device),
                               UCL_WRITE_OPTIMIZED);
                               
  for (int ix=0; ix<nrhor; ix++)
    for (int iy=0; iy<nr+1; iy++)
      for (int iz=0; iz<7; iz++) 
    dview_rhor_spline[ix*(nr+1)*7+iy*7+iz]=host_rhor_spline[ix][iy][iz];
  
  rhor_spline.alloc(nrhor*(nr+1)*7,*(this->ucl_device),UCL_READ_ONLY);
  ucl_copy(rhor_spline,dview_rhor_spline,false);
  
  // pack z2r_spline
  UCL_H_Vec<numtyp> dview_z2r_spline(nz2r*(nr+1)*7,*(this->ucl_device),
                               UCL_WRITE_OPTIMIZED);
                               
  for (int ix=0; ix<nz2r; ix++)
    for (int iy=0; iy<nr+1; iy++)
      for (int iz=0; iz<7; iz++) 
    dview_z2r_spline[ix*(nr+1)*7+iy*7+iz]=host_z2r_spline[ix][iy][iz];
  
  z2r_spline.alloc(nz2r*(nr+1)*7,*(this->ucl_device),UCL_READ_ONLY);
  ucl_copy(z2r_spline,dview_z2r_spline,false);

  _allocated=true;
  this->_max_bytes=type2rhor_z2r.row_bytes()+
        rhor_spline.row_bytes()+z2r_spline.row_bytes();
  return 0;
}

template <class numtyp, class acctyp>
void EAMT::clear() {
  if (!_allocated)
    return;
  _allocated=false;
  
  type2rhor_z2r.clear();
  rhor_spline.clear();
  z2r_spline.clear();
  this->clear_atomic();
}

template <class numtyp, class acctyp>
double EAMT::host_memory_usage() const {
  return this->host_memory_usage_atomic()+sizeof(EAM<numtyp,acctyp>);
}

// ---------------------------------------------------------------------------
// Calculate energies, forces, and torques
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void EAMT::loop(const bool _eflag, const bool _vflag) {
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
    this->k_pair_fast.set_size(GX,BX);
    this->k_pair_fast.run(&this->atom->dev_x.begin(), &this->atom->dev_q.begin(), 
                   &type2rhor_z2r.begin(),
                   &rhor_spline.begin(), &z2r_spline.begin(),
                   &this->nbor->dev_nbor.begin(),
                   &this->_nbor_data->begin(), &this->ans->dev_ans.begin(),
                   &this->ans->dev_engv.begin(), &eflag, &vflag, &ainum,
                   &nbor_pitch, &_cutforcesq, &_rdr,
                   &_nrhor, &_nz2r, &_nr,
                   &this->_threads_per_atom);
  } else {
    this->k_pair.set_size(GX,BX);
    this->k_pair.run(&this->atom->dev_x.begin(), &this->atom->dev_q.begin(), 
                   &type2rhor_z2r.begin(),
                   &rhor_spline.begin(), &z2r_spline.begin(),
                   &this->nbor->dev_nbor.begin(),
                   &this->_nbor_data->begin(), &this->ans->dev_ans.begin(),
                   &this->ans->dev_engv.begin(), &eflag, &vflag, &ainum,
                   &nbor_pitch, &_ntypes, &_cutforcesq, &_rdr,
                   &_nrhor, &_nz2r, &_nr,
                   &this->_threads_per_atom);
  }

  this->time_pair.stop();
}

template class EAM<PRECISION,ACC_PRECISION>;
