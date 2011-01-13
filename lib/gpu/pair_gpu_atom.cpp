/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
 
/* ----------------------------------------------------------------------
   Contributing authors: Mike Brown (ORNL), brownw@ornl.gov
------------------------------------------------------------------------- */

#include "pair_gpu_atom.h"

#define PairGPUAtomT PairGPUAtom<numtyp,acctyp>

#ifdef WINDLL
#include <windows.h>
typedef bool (*__win_sort_alloc)(const int max_atoms);
typedef void (*__win_sort)(const int max_atoms, unsigned *cell_begin,
                           int *particle_begin);
__win_sort_alloc _win_sort_alloc;
__win_sort _win_sort;
#endif

template <class numtyp, class acctyp>
PairGPUAtomT::PairGPUAtom() : _compiled(false),_allocated(false),_eflag(false),
                              _vflag(false),_inum(0),_ilist(NULL), 
                              _newton(false) {
  #ifndef USE_OPENCL
  sort_config.op = CUDPP_ADD;
  sort_config.datatype = CUDPP_UINT;
  sort_config.algorithm = CUDPP_SORT_RADIX;
  sort_config.options = CUDPP_OPTION_KEY_VALUE_PAIRS;

  #ifdef WINDLL
  HINSTANCE hinstLib = LoadLibrary(TEXT("gpu.dll"));
  if (hinstLib == NULL) {
    printf("\nUnable to load gpu.dll\n");
    exit(1);
  }
  _win_sort_alloc=(__win_sort_alloc)GetProcAddress(hinstLib,"_win_sort_alloc");
  _win_sort=(__win_sort)GetProcAddress(hinstLib,"_win_sort");
  #endif

  #endif
}

template <class numtyp, class acctyp>
int PairGPUAtomT::bytes_per_atom() const { 
  int id_space=0;
  if (_gpu_nbor)
    id_space=2;
  int bytes=4*sizeof(numtyp)+11*sizeof(acctyp)+id_space;
  if (_rot)
    bytes+=4*sizeof(numtyp)+4*sizeof(acctyp);
  if (_charge)
    bytes+=sizeof(numtyp);
  return bytes;
}

template <class numtyp, class acctyp>
bool PairGPUAtomT::alloc(const int inum, const int nall) {
  _max_atoms=static_cast<int>(static_cast<double>(nall)*1.10);
  if (_newton)
    _max_local=_max_atoms;
  else
    _max_local=static_cast<int>(static_cast<double>(inum)*1.10);

  bool success=true;
  
  int ans_elements=4;
  if (_rot)
    ans_elements+=4;
  
  // Ignore host/device transfers?
  bool cpuview=false;
  if (dev->device_type()==UCL_CPU)
    cpuview=true;
    
  // Allocate storage for CUDPP sort
  #ifndef USE_OPENCL
  #ifdef WINDLL
  _win_sort_alloc(_max_atoms);
  #else
  if (_gpu_nbor) {
    CUDPPResult result = cudppPlan(&sort_plan, sort_config, _max_atoms, 1, 0);  
    if (CUDPP_SUCCESS != result)
      return false;
  }
  #endif
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
  success=success &&(host_ans.alloc(ans_elements*_max_local,*dev)==UCL_SUCCESS);
  success=success &&(host_engv.alloc(_ev_fields*_max_local,*dev)==UCL_SUCCESS);
  // Buffer for casting only if different precisions
  if (_charge)
    success=success && (host_q.alloc(_max_atoms,*dev,
                                     UCL_WRITE_OPTIMIZED)==UCL_SUCCESS);
  // Buffer for casting only if different precisions
  if (_rot)
    success=success && (host_quat.alloc(_max_atoms*4,*dev,
                                        UCL_WRITE_OPTIMIZED)==UCL_SUCCESS);

    
  // ---------------------------  Device allocations
  _gpu_bytes=0;
  if (cpuview) {
    #ifdef GPU_CAST
    assert(0==1);
    #else
    dev_x.view(host_x);
    #endif
    dev_engv.view(host_engv);
    dev_ans.view(host_ans);
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
    _gpu_bytes+=dev_x_cast.row_bytes()+dev_type_cast.row_bytes();
    #else
    success=success && (UCL_SUCCESS==
                        dev_x.alloc(_max_atoms*4,*dev,UCL_READ_ONLY));
    #endif
    success=success && (dev_engv.alloc(_ev_fields*_max_local,*dev,
                                       UCL_WRITE_ONLY)==UCL_SUCCESS);
    success=success && (dev_ans.alloc(ans_elements*_max_local,
                                      *dev,UCL_WRITE_ONLY)==UCL_SUCCESS);
    if (_charge) {
      success=success && (dev_q.alloc(_max_atoms,*dev,
                                      UCL_READ_ONLY)==UCL_SUCCESS);
      _gpu_bytes+=dev_q.row_bytes();
    }
    if (_rot) {
      success=success && (dev_quat.alloc(_max_atoms*4,*dev,
                                      UCL_READ_ONLY)==UCL_SUCCESS);
      _gpu_bytes+=dev_quat.row_bytes();
    }
  }
  if (_gpu_nbor) {
    success=success && (dev_cell_id.alloc(_max_atoms,*dev)==UCL_SUCCESS);
    success=success && (dev_particle_id.alloc(_max_atoms,*dev)==UCL_SUCCESS);
    _gpu_bytes+=dev_cell_id.row_bytes()+dev_particle_id.row_bytes();
    if (_bonds) {
      success=success && (dev_tag.alloc(_max_atoms,*dev)==UCL_SUCCESS);
      _gpu_bytes+=dev_tag.row_bytes();
    }
  }

  _gpu_bytes+=dev_x.row_bytes()+dev_engv.row_bytes()+dev_ans.row_bytes();
  
  _allocated=true;  
  return success;
}

template <class numtyp, class acctyp>
bool PairGPUAtomT::init(const int inum, const int nall, const bool charge,
                        const bool rot, UCL_Device &devi, const bool gpu_nbor,
                        const bool bonds) {
  clear();

  bool success=true;
  _gpu_nbor=gpu_nbor;
  _bonds=bonds;
  _charge=charge;
  _rot=rot;
  _other=_charge || _rot;
  dev=&devi;

  _e_fields=1;
  if (_charge)
    _e_fields++;
  _ev_fields=6+_e_fields;
    
  // Initialize atom and nbor data
  int ef_inum=inum;
  if (ef_inum==0)
    ef_inum=1000;
  int ef_nall=nall;
  if (ef_nall<=ef_inum)
    ef_nall=ef_inum*2;
  
  // Initialize timers for the selected device
  time_pos.init(*dev);
  time_other.init(*dev);
  time_answer.init(*dev);
  time_pos.zero();
  time_other.zero();
  time_answer.zero();
  _time_cast=0.0;
  
  #ifdef GPU_CAST
  compile_kernels(*dev);
  #endif
  
  return success && alloc(ef_inum,ef_nall);
}
  
template <class numtyp, class acctyp>
void PairGPUAtomT::clear_resize() {
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
  dev_ans.clear();
  dev_engv.clear();
  #ifndef GPU_CAST
  host_x.clear();
  #else
  host_x_cast.clear();
  host_type_cast.clear();
  #endif
  host_ans.clear();
  host_engv.clear();
  dev_cell_id.clear();
  dev_particle_id.clear();
  dev_tag.clear();
  #ifdef GPU_CAST
  dev_x_cast.clear();
  dev_type_cast.clear();
  #endif

  #ifndef USE_OPENCL
  #ifndef WINDLL
  if (_gpu_nbor) cudppDestroyPlan(sort_plan);
  #endif
  #endif
}

template <class numtyp, class acctyp>
void PairGPUAtomT::clear() {
  _gpu_bytes=0;
  if (!_allocated)
    return;

  time_pos.clear();
  time_other.clear();
  time_answer.clear();
  clear_resize();
  _inum=0;
  _eflag=false;
  _vflag=false;

  #ifdef GPU_CAST
  if (_compiled) {
    k_cast_x.clear();
    delete atom_program;
    _compiled=false;
  }
  #endif
}

template <class numtyp, class acctyp>
double PairGPUAtomT::host_memory_usage() const {
  int atom_bytes=4;
  if (_charge) 
    atom_bytes+=1;
  if (_rot) 
    atom_bytes+=4;
  int ans_bytes=atom_bytes+_ev_fields;
  return _max_atoms*atom_bytes*sizeof(numtyp)+
         ans_bytes*(_max_local)*sizeof(acctyp)+
         sizeof(PairGPUAtom<numtyp,acctyp>);
}
  
template <class numtyp, class acctyp>
void PairGPUAtomT::copy_answers(const bool eflag, const bool vflag,
                                const bool ef_atom, const bool vf_atom) {
  time_answer.start();
  _eflag=eflag;
  _vflag=vflag;
  _ef_atom=ef_atom;
  _vf_atom=vf_atom;
    
  int csize=_ev_fields;    
  if (!eflag)
    csize-=_e_fields;
  if (!vflag)
    csize-=6;
      
  if (csize>0)
    ucl_copy(host_engv,dev_engv,_inum*csize,true);
  if (_rot)
    ucl_copy(host_ans,dev_ans,_inum*4*2,true);
  else
    ucl_copy(host_ans,dev_ans,_inum*4,true);
  time_answer.stop();
}

template <class numtyp, class acctyp>
void PairGPUAtomT::copy_answers(const bool eflag, const bool vflag,
                                const bool ef_atom, const bool vf_atom,
                                int *ilist) {
  _ilist=ilist;
  copy_answers(eflag,vflag,ef_atom,vf_atom);
}

template <class numtyp, class acctyp>
double PairGPUAtomT::energy_virial(double *eatom, double **vatom,
                                   double *virial) {
  if (_eflag==false && _vflag==false)
    return 0.0;

  double evdwl=0.0;
  if (_gpu_nbor) {
    for (int i=0; i<_inum; i++) {
      acctyp *ap=host_engv.begin()+i;
      if (_eflag) {
        if (_ef_atom) {
          evdwl+=*ap;
          eatom[i]+=*ap*0.5;
          ap+=_inum;
        } else {
          evdwl+=*ap;
          ap+=_inum;
        }
      }
      if (_vflag) {
        if (_vf_atom) {
          for (int j=0; j<6; j++) {
            vatom[i][j]+=*ap*0.5;
            virial[j]+=*ap;
            ap+=_inum;
          }
        } else {
          for (int j=0; j<6; j++) {
            virial[j]+=*ap;
            ap+=_inum;
          }
        }
      }
    }
    for (int j=0; j<6; j++)
      virial[j]*=0.5;
  } else {
    for (int i=0; i<_inum; i++) {
      acctyp *ap=host_engv.begin()+i;
      int ii=_ilist[i];
      if (_eflag) {
        if (_ef_atom) {
          evdwl+=*ap;
          eatom[ii]+=*ap*0.5;
          ap+=_inum;
        } else {
          evdwl+=*ap;
          ap+=_inum;
        }
      }
      if (_vflag) {
        if (_vf_atom) {
          for (int j=0; j<6; j++) {
            vatom[ii][j]+=*ap*0.5;
            virial[j]+=*ap;
            ap+=_inum;
          }
        } else {
          for (int j=0; j<6; j++) {
            virial[j]+=*ap;
            ap+=_inum;
          }
        }
      }
    }
    for (int j=0; j<6; j++)
      virial[j]*=0.5;
  }
  
  evdwl*=0.5;
  return evdwl;
}

template <class numtyp, class acctyp>
double PairGPUAtomT::energy_virial(double *eatom, double **vatom,
                                   double *virial, double &ecoul) {
  if (_eflag==false && _vflag==false) {
    ecoul=0.0;
    return 0.0;
  }

  if (_charge==false)
    return energy_virial(eatom,vatom,virial);

  double evdwl=0.0;
  double _ecoul=0.0;
  if (_gpu_nbor) {
    for (int i=0; i<_inum; i++) {
      acctyp *ap=host_engv.begin()+i;
      if (_eflag) {
        if (_ef_atom) {
          evdwl+=*ap;
          eatom[i]+=*ap*0.5;
          ap+=_inum;
          _ecoul+=*ap;
          eatom[i]+=*ap*0.5;
          ap+=_inum;
        } else {
          evdwl+=*ap;
          ap+=_inum;
          _ecoul+=*ap;
          ap+=_inum;
        }
      }
      if (_vflag) {
        if (_vf_atom) {
          for (int j=0; j<6; j++) {
            vatom[i][j]+=*ap*0.5;
            virial[j]+=*ap;
            ap+=_inum;
          }
        } else {
          for (int j=0; j<6; j++) {
            virial[j]+=*ap;
            ap+=_inum;
          }
        }
      }
    }
    for (int j=0; j<6; j++)
      virial[j]*=0.5;
  } else {
    for (int i=0; i<_inum; i++) {
      acctyp *ap=host_engv.begin()+i;
      int ii=_ilist[i];
      if (_eflag) {
        if (_ef_atom) {
          evdwl+=*ap;
          eatom[ii]+=*ap*0.5;
          ap+=_inum;
          _ecoul+=*ap;
          eatom[ii]+=*ap*0.5;
          ap+=_inum;
        } else {
          evdwl+=*ap;
          ap+=_inum;
          _ecoul+=*ap;
          ap+=_inum;
        }
      }
      if (_vflag) {
        if (_vf_atom) {
          for (int j=0; j<6; j++) {
            vatom[ii][j]+=*ap*0.5;
            virial[j]+=*ap;
            ap+=_inum;
          }
        } else {
          for (int j=0; j<6; j++) {
            virial[j]+=*ap;
            ap+=_inum;
          }
        }
      }
    }
    for (int j=0; j<6; j++)
      virial[j]*=0.5;
  }
  
  evdwl*=0.5;
  ecoul+=_ecoul*0.5;
  return evdwl;
}

template <class numtyp, class acctyp>
void PairGPUAtomT::get_answers(double **f, double **tor) {
  acctyp *ap=host_ans.begin();
  if (_gpu_nbor) {
    for (int i=0; i<_inum; i++) {
      f[i][0]+=*ap;
      ap++;
      f[i][1]+=*ap;
      ap++;
      f[i][2]+=*ap;
      ap+=2;
    }
    if (_rot) {
      for (int i=0; i<_inum; i++) {
        tor[i][0]+=*ap;
        ap++;
        tor[i][1]+=*ap;
        ap++;
        tor[i][2]+=*ap;
        ap+=2;
      }
    }
  } else {
    for (int i=0; i<_inum; i++) {
      int ii=_ilist[i];
      f[ii][0]+=*ap;
      ap++;
      f[ii][1]+=*ap;
      ap++;
      f[ii][2]+=*ap;
      ap+=2;
    }
    if (_rot) {
      for (int i=0; i<_inum; i++) {
        int ii=_ilist[i];
        tor[ii][0]+=*ap;
        ap++;
        tor[ii][1]+=*ap;
        ap++;
        tor[ii][2]+=*ap;
        ap+=2;
      }
    }
  }
}

// Sort arrays for neighbor list calculation
template <class numtyp, class acctyp>
void PairGPUAtomT::sort_neighbor(const int num_atoms) {
  #ifndef USE_OPENCL
  #ifdef WINDLL
  _win_sort(num_atoms,(unsigned *)dev_cell_id.begin(),
            (int *)dev_particle_id.begin());
  #else
  CUDPPResult result = cudppSort(sort_plan, (unsigned *)dev_cell_id.begin(), 
                                 (int *)dev_particle_id.begin(), 
                                 8*sizeof(unsigned), num_atoms);
  if (CUDPP_SUCCESS != result) {
    printf("Error in cudppSort\n");
    NVD_GERYON_EXIT;
  }
  #endif
  #endif
}

#ifdef GPU_CAST
#ifdef USE_OPENCL
#include "pair_gpu_atom_cl.h"
#else
#include "pair_gpu_atom_ptx.h"
#endif

template <class numtyp, class acctyp>
void PairGPUAtomT::compile_kernels(UCL_Device &dev) {
  atom_program=new UCL_Program(dev);
  atom_program->load_string(pair_gpu_atom_kernel,"");
  k_cast_x.set_function(*atom_program,"kernel_cast_x");
  _compiled=true;
}

#endif

template class PairGPUAtom<PRECISION,ACC_PRECISION>;
