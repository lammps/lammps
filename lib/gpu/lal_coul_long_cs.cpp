/***************************************************************************
                              coul_long_cs.cpp
                             -------------------
                           Trung Nguyen (Northwestern)

  Class for acceleration of the coul/long pair style.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : June 2018
    email                : ndactrung@gmail.com
 ***************************************************************************/

#if defined(USE_OPENCL)
#include "coul_long_cs_cl.h"
#elif defined(USE_CUDART)
const char *coul_long_cs=0;
#else
#include "coul_long_cs_cubin.h"
#endif

#include "lal_coul_long_cs.h"
#include <cassert>
namespace LAMMPS_AL {
#define CoulLongCST CoulLongCS<numtyp, acctyp>

extern Device<PRECISION,ACC_PRECISION> pair_gpu_device;

template <class numtyp, class acctyp>
int CoulLongCST::init(const int ntypes, double **host_scale,
                    const int nlocal, const int nall, const int max_nbors,
                    const int maxspecial, const double cell_size,
                    const double gpu_split, FILE *_screen,
                    const double host_cut_coulsq, double *host_special_coul,
                    const double qqrd2e, const double g_ewald) {
  int success;
  success=this->init_atomic(nlocal,nall,max_nbors,maxspecial,cell_size,
                            gpu_split,_screen,coul_long_cs,"k_coul_long_cs");
  if (success!=0)
    return success;

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

  this->scale.alloc(lj_types*lj_types,*(this->ucl_device),UCL_READ_ONLY);
  this->atom->type_pack1(ntypes,lj_types,this->scale,host_write,host_scale);

  this->sp_cl.alloc(4,*(this->ucl_device),UCL_READ_ONLY);
  for (int i=0; i<4; i++) {
    host_write[i]=host_special_coul[i];
  }
  ucl_copy(this->sp_cl,host_write,4,false);

  this->_cut_coulsq=host_cut_coulsq;
  this->_qqrd2e=qqrd2e;
  this->_g_ewald=g_ewald;

  this->_allocated=true;
  this->_max_bytes=this->scale.row_bytes()+this->sp_cl.row_bytes();
  return 0;
}

template class CoulLongCS<PRECISION,ACC_PRECISION>;
}
