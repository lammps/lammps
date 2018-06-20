/***************************************************************************
                            born_coul_long_cs.cpp
                             -------------------
                           Trung Dac Nguyen (Northwestern)

  Class for acceleration of the born/coul/long/cs pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : ndactrung@gmail.com
 ***************************************************************************/

#ifdef USE_OPENCL
#include "born_coul_long_cs_cl.h"
#elif defined(USE_CUDART)
const char *born_coul_long_cs=0;
#else
#include "born_coul_long_cs_cubin.h"
#endif

#include "lal_born_coul_long_cs.h"
#include <cassert>
using namespace LAMMPS_AL;
#define BornCoulLongCST BornCoulLongCS<numtyp, acctyp>

extern Device<PRECISION,ACC_PRECISION> device;

template <class numtyp, class acctyp>
int BornCoulLongCST::init(const int ntypes, double **host_cutsq, double **host_rhoinv,
                       double **host_born1, double **host_born2, double **host_born3,
                       double **host_a, double **host_c, double **host_d,
                       double **host_sigma, double **host_offset,
                       double *host_special_lj, const int nlocal,
                       const int nall, const int max_nbors,
                       const int maxspecial, const double cell_size,
                       const double gpu_split, FILE *_screen,
                       double **host_cut_ljsq, const double host_cut_coulsq,
                       double *host_special_coul, const double qqrd2e,
                       const double g_ewald) {
  int success;
  success=this->init_atomic(nlocal,nall,max_nbors,maxspecial,cell_size,gpu_split,
                            _screen,born_coul_long_cs,"k_born_coul_long_cs");
  if (success!=0)
    return success;

  // If atom type constants fit in shared memory use fast kernel
  int lj_types=ntypes;
  this->shared_types=false;
  int max_shared_types=this->device->max_shared_types();
  if (lj_types<=max_shared_types && this->_block_size>=max_shared_types) {
    lj_types=max_shared_types;
    this->shared_types=true;
  }
  this->_lj_types=lj_types;

  // Allocate a host write buffer for data initialization
  UCL_H_Vec<numtyp> host_write(lj_types*lj_types*32,*(this->ucl_device),
                               UCL_WRITE_ONLY);

  for (int i=0; i<lj_types*lj_types; i++)
    host_write[i]=0.0;

  this->coeff1.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack4(ntypes,lj_types,this->coeff1,host_write,host_rhoinv,
                         host_born1,host_born2,host_born3);

  this->coeff2.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack4(ntypes,lj_types,this->coeff2,host_write,host_a,host_c,
                         host_d,host_offset);

  this->cutsq_sigma.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack4(ntypes,lj_types,this->cutsq_sigma,host_write,host_cutsq,
             host_cut_ljsq,host_sigma);

  this->sp_lj.alloc(8,*(this->ucl_device),UCL_READ_ONLY);
  for (int i=0; i<4; i++) {
    host_write[i]=host_special_lj[i];
    host_write[i+4]=host_special_coul[i];
  }
  ucl_copy(this->sp_lj,host_write,8,false);

  this->_cut_coulsq=host_cut_coulsq;
  this->_qqrd2e=qqrd2e;
  this->_g_ewald=g_ewald;

  this->_allocated=true;
  this->_max_bytes=this->coeff1.row_bytes()+this->coeff2.row_bytes()
      +this->cutsq_sigma.row_bytes()+this->sp_lj.row_bytes();
  return 0;
}

template class BornCoulLongCS<PRECISION,ACC_PRECISION>;
