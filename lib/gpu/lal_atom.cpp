/***************************************************************************
                                  atom.cpp
                             -------------------
                            W. Michael Brown (ORNL)

  Class for particle data management

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : 
    email                : brownw@ornl.gov
 ***************************************************************************/

#include "lal_atom.h"

using namespace LAMMPS_AL;
#define AtomT Atom<numtyp,acctyp>

template <class numtyp, class acctyp>
AtomT::Atom() : _compiled(false),_allocated(false),
                              _max_gpu_bytes(0) {
  #ifdef USE_CUDPP
  sort_config.op = CUDPP_ADD;
  sort_config.datatype = CUDPP_UINT;
  sort_config.algorithm = CUDPP_SORT_RADIX;
  sort_config.options = CUDPP_OPTION_KEY_VALUE_PAIRS;
  #endif
}

template <class numtyp, class acctyp>
int AtomT::bytes_per_atom() const { 
  int id_space=0;
  if (_gpu_nbor==1)
    id_space=2;
  else if (_gpu_nbor==2)
    id_space=4;
  int bytes=4*sizeof(numtyp)+id_space*sizeof(int);
  if (_rot)
    bytes+=4*sizeof(numtyp);
  if (_charge)
    bytes+=sizeof(numtyp);
  return bytes;
}

template <class numtyp, class acctyp>
bool AtomT::alloc(const int nall) {
  _max_atoms=static_cast<int>(static_cast<double>(nall)*1.10);

  bool success=true;
  
  // Ignore host/device transfers?
  bool cpuview=false;
  if (dev->device_type()==UCL_CPU)
    cpuview=true;
    
  // Allocate storage for CUDPP sort
  #ifdef USE_CUDPP
  if (_gpu_nbor==1) {
    CUDPPResult result = cudppPlan(&sort_plan, sort_config, _max_atoms, 1, 0);  
    if (CUDPP_SUCCESS != result)
      return false;
  }
  #endif

  // --------------------------   Host allocations
  // Get a host write only buffer
  #ifdef GPU_CAST
  success=success && (host_x_cast.alloc(_max_atoms*3,*dev,
                                        UCL_WRITE_OPTIMIZED)==UCL_SUCCESS);
  success=success && (host_type_cast.alloc(_max_atoms,*dev,
                                           UCL_WRITE_OPTIMIZED)==UCL_SUCCESS);
  #else
  success=success && (host_x.alloc(_max_atoms*4,*dev,
                      UCL_WRITE_OPTIMIZED)==UCL_SUCCESS);
  #endif                      
  // Buffer for casting only if different precisions
  if (_charge)
    success=success && (host_q.alloc(_max_atoms,*dev,
                                     UCL_WRITE_OPTIMIZED)==UCL_SUCCESS);
  // Buffer for casting only if different precisions
  if (_rot)
    success=success && (host_quat.alloc(_max_atoms*4,*dev,
                                        UCL_WRITE_OPTIMIZED)==UCL_SUCCESS);

    
  // ---------------------------  Device allocations
  int gpu_bytes=0;
  if (cpuview) {
    #ifdef GPU_CAST
    assert(0==1);
    #else
    dev_x.view(host_x);
    #endif
    if (_rot)
      dev_quat.view(host_quat);
    if (_charge)
      dev_q.view(host_q);
  } else {
    #ifdef GPU_CAST
    success=success && (UCL_SUCCESS==dev_x.alloc(_max_atoms*4,*dev));
    success=success && (UCL_SUCCESS==
                        dev_x_cast.alloc(_max_atoms*3,*dev,UCL_READ_ONLY));
    success=success && (UCL_SUCCESS==
                        dev_type_cast.alloc(_max_atoms,*dev,UCL_READ_ONLY));
    gpu_bytes+=dev_x_cast.row_bytes()+dev_type_cast.row_bytes();
    #else
    success=success && (UCL_SUCCESS==
                        dev_x.alloc(_max_atoms*4,*dev,UCL_READ_ONLY));
    #endif
    if (_charge) {
      success=success && (dev_q.alloc(_max_atoms,*dev,
                                      UCL_READ_ONLY)==UCL_SUCCESS);
      gpu_bytes+=dev_q.row_bytes();
    }
    if (_rot) {
      success=success && (dev_quat.alloc(_max_atoms*4,*dev,
                                      UCL_READ_ONLY)==UCL_SUCCESS);
      gpu_bytes+=dev_quat.row_bytes();
    }
  }
  if (_gpu_nbor>0) {
    success=success && (dev_particle_id.alloc(_max_atoms,*dev)==UCL_SUCCESS);
    gpu_bytes+=dev_particle_id.row_bytes();
    if (_bonds) {
      success=success && (dev_tag.alloc(_max_atoms,*dev)==UCL_SUCCESS);
      gpu_bytes+=dev_tag.row_bytes();
    }
    if (_gpu_nbor==1) {
      success=success && (dev_cell_id.alloc(_max_atoms,*dev)==UCL_SUCCESS);
      gpu_bytes+=dev_cell_id.row_bytes();
    } else {
      success=success && (host_particle_id.alloc(_max_atoms,*dev)==UCL_SUCCESS);
      success=success && 
             (host_cell_id.alloc(_max_atoms,*dev,UCL_NOT_PINNED)==UCL_SUCCESS);
    }             
  }

  gpu_bytes+=dev_x.row_bytes();
  if (gpu_bytes>_max_gpu_bytes)
    _max_gpu_bytes=gpu_bytes;
  
  _allocated=true;  
  return success;
}

template <class numtyp, class acctyp>
bool AtomT::add_fields(const bool charge, const bool rot,
                       const int gpu_nbor, const bool bonds) {
  bool realloc=false;
  if (charge && _charge==false) {
    _charge=true;
    realloc=true;
  }
  if (rot && _rot==false) {
    _rot=true;
    realloc=true;
  }
  if (gpu_nbor>0 && _gpu_nbor==0) {
    _gpu_nbor=gpu_nbor;
    realloc=true;
  }
  if (bonds && _bonds==false) {
    _bonds=true;
    realloc=true;
  }
  if (realloc) {
    _other=_charge || _rot;
    int max_atoms=_max_atoms;
    clear_resize();
    return alloc(max_atoms);
  }
  return true;
}

template <class numtyp, class acctyp>
bool AtomT::init(const int nall, const bool charge, const bool rot,
                 UCL_Device &devi, const int gpu_nbor, const bool bonds) {
  clear();

  bool success=true;
  _x_avail=false;
  _q_avail=false;
  _quat_avail=false;
  _resized=false;
  _gpu_nbor=gpu_nbor;
  _bonds=bonds;
  _charge=charge;
  _rot=rot;
  _other=_charge || _rot;
  dev=&devi;
  _time_transfer=0;

  // Initialize atom and nbor data
  int ef_nall=nall;
  if (ef_nall==0)
    ef_nall=2000;
  
  // Initialize timers for the selected device
  time_pos.init(*dev);
  time_q.init(*dev);
  time_quat.init(*dev);
  time_pos.zero();
  time_q.zero();
  time_quat.zero();
  _time_cast=0.0;
  
  #ifdef GPU_CAST
  compile_kernels(*dev);
  #endif
  
  return success && alloc(ef_nall);
}
  
template <class numtyp, class acctyp>
void AtomT::clear_resize() {
  if (!_allocated)
    return;
  _allocated=false;

  dev_x.clear();
  if (_charge) { 
    dev_q.clear();
    host_q.clear();
  }
  if (_rot) {
    dev_quat.clear();
    host_quat.clear();
  }
  #ifndef GPU_CAST
  host_x.clear();
  #else
  host_x_cast.clear();
  host_type_cast.clear();
  #endif
  dev_cell_id.clear();
  dev_particle_id.clear();
  dev_tag.clear();
  #ifdef GPU_CAST
  dev_x_cast.clear();
  dev_type_cast.clear();
  #endif

  #ifdef USE_CUDPP
  if (_gpu_nbor==1) cudppDestroyPlan(sort_plan);
  #endif
  
  if (_gpu_nbor==2) {
    host_particle_id.clear();
    host_cell_id.clear();
  }
}

template <class numtyp, class acctyp>
void AtomT::clear() {
  _max_gpu_bytes=0;
  if (!_allocated)
    return;

  time_pos.clear();
  time_q.clear();
  time_quat.clear();
  clear_resize();

  #ifdef GPU_CAST
  if (_compiled) {
    k_cast_x.clear();
    delete atom_program;
    _compiled=false;
  }
  #endif
}

template <class numtyp, class acctyp>
double AtomT::host_memory_usage() const {
  int atom_bytes=4;
  if (_charge) 
    atom_bytes+=1;
  if (_rot) 
    atom_bytes+=4;
  return _max_atoms*atom_bytes*sizeof(numtyp)+
         sizeof(Atom<numtyp,acctyp>);
}
  
// Sort arrays for neighbor list calculation
template <class numtyp, class acctyp>
void AtomT::sort_neighbor(const int num_atoms) {
  #ifdef USE_CUDPP
  CUDPPResult result = cudppSort(sort_plan, (unsigned *)dev_cell_id.begin(), 
                                 (int *)dev_particle_id.begin(), 
                                 8*sizeof(unsigned), num_atoms);
  if (CUDPP_SUCCESS != result) {
    printf("Error in cudppSort\n");
    NVD_GERYON_EXIT;
  }
  #endif
}

#ifdef GPU_CAST
#ifdef USE_OPENCL
#include "atom_cl.h"
#else
#include "atom_ptx.h"
#endif

template <class numtyp, class acctyp>
void AtomT::compile_kernels(UCL_Device &dev) {
  std::string flags = "-D"+std::string(OCL_VENDOR);
  atom_program=new UCL_Program(dev);
  atom_program->load_string(atom,flags);
  k_cast_x.set_function(*atom_program,"kernel_cast_x");
  _compiled=true;
}

#endif

template class Atom<PRECISION,ACC_PRECISION>;
