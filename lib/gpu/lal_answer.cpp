/***************************************************************************
                                 answer.cpp
                             -------------------
                            W. Michael Brown (ORNL)

  Class for data management of forces, torques, energies, and virials

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : brownw@ornl.gov
 ***************************************************************************/

#include "lal_answer.h"
#if (LAL_USE_OMP == 1)
#include <omp.h>
#endif

namespace LAMMPS_AL {
#define AnswerT Answer<numtyp,acctyp>

template <class numtyp, class acctyp>
AnswerT::Answer() : _allocated(false),_eflag(false),_vflag(false),
                    _inum(0),_ilist(nullptr),_newton(false) {
}

template <class numtyp, class acctyp>
int AnswerT::bytes_per_atom() const {
  int bytes=11*sizeof(acctyp);
  if (_rot)
    bytes+=4*sizeof(acctyp);
  if (_charge)
    bytes+=sizeof(acctyp);
  return bytes;
}

template <class numtyp, class acctyp>
bool AnswerT::alloc(const int inum) {
  _max_local=static_cast<int>(static_cast<double>(inum)*1.10);

  bool success=true;

  _ans_fields=4;
  if (_rot)
    _ans_fields+=4;

  // ---------------------------  Device allocations
  success=success && (engv.alloc(_ev_fields*_max_local,*dev,UCL_READ_ONLY,
                                 UCL_READ_WRITE)==UCL_SUCCESS);
  success=success && (force.alloc(_ans_fields*_max_local,*dev,UCL_READ_ONLY,
                                UCL_READ_WRITE)==UCL_SUCCESS);
  _gpu_bytes=engv.device.row_bytes()+force.device.row_bytes();

  _allocated=true;
  return success;
}

template <class numtyp, class acctyp>
bool AnswerT::init(const int inum, const bool charge, const bool rot,
                   UCL_Device &devi) {
  clear();

  bool success=true;
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

  // Initialize timers for the selected device
  time_answer.init(*dev);
  time_answer.zero();
  _time_cast=0.0;
  _time_cpu_idle=0.0;

  success=success && (error_flag.alloc(1,*dev,UCL_READ_WRITE,
                                        UCL_WRITE_ONLY)==UCL_SUCCESS);
  if (success) error_flag.zero();

  return success && alloc(ef_inum);
}

template <class numtyp, class acctyp>
bool AnswerT::add_fields(const bool charge, const bool rot) {
  bool realloc=false;
  if (charge && _charge==false) {
    _charge=true;
    _e_fields++;
    _ev_fields++;
    realloc=true;
  }
  if (rot && _rot==false) {
    _rot=true;
    realloc=true;
  }
  if (realloc) {
    _other=_charge || _rot;
    int inum=_max_local;
    force.clear();
    engv.clear();
    _allocated=false;
    return alloc(inum);
  }
  return true;
}

template <class numtyp, class acctyp>
void AnswerT::clear() {
  _gpu_bytes=0;
  error_flag.clear();
  if (!_allocated)
    return;
  _allocated=false;

  force.clear();
  engv.clear();
  time_answer.clear();
  _inum=0;
  _ilist=nullptr;
  _eflag=false;
  _vflag=false;
}

template <class numtyp, class acctyp>
double AnswerT::host_memory_usage() const {
  int atom_bytes=4;
  if (_charge)
    atom_bytes+=1;
  if (_rot)
    atom_bytes+=4;
  int ans_bytes=atom_bytes+_ev_fields;
  return ans_bytes*(_max_local)*sizeof(acctyp)+
         sizeof(Answer<numtyp,acctyp>);
}

template <class numtyp, class acctyp>
void AnswerT::copy_answers(const bool eflag, const bool vflag,
                           const bool ef_atom, const bool vf_atom,
                           const int red_blocks) {
  time_answer.start();
  _eflag=eflag;
  _vflag=vflag;
  _ef_atom=ef_atom;
  _vf_atom=vf_atom;
  #ifdef LAL_NO_BLOCK_REDUCE
  _ev_stride=_inum;
  #else
  if (ef_atom || vf_atom)
    _ev_stride=_inum;
  else
    _ev_stride=red_blocks;
  #endif

  int csize=_ev_fields;
  if (!eflag)
    csize-=_e_fields;
  if (!vflag)
    csize-=6;

  if (csize>0)
    engv.update_host(_ev_stride*csize,true);
  if (_rot)
    force.update_host(_inum*4*2,true);
  else
    force.update_host(_inum*4,true);
  time_answer.stop();

  #ifndef GERYON_OCL_FLUSH
  force.flush();
  #endif
}

template <class numtyp, class acctyp>
void AnswerT::copy_answers(const bool eflag, const bool vflag,
                           const bool ef_atom, const bool vf_atom,
                           int *ilist, const int red_blocks) {
  _ilist=ilist;
  copy_answers(eflag,vflag,ef_atom,vf_atom,red_blocks);
}

template <class numtyp, class acctyp>
double AnswerT::energy_virial(double *eatom, double **vatom,
                                  double *virial) {
  if (_eflag==false && _vflag==false)
    return 0.0;

  double evdwl=0.0;
  int vstart=0;
  if (_eflag) {
    #if (LAL_USE_OMP_SIMD == 1)
    #pragma omp simd reduction(+:evdwl)
    #endif
    for (int i=0; i<_ev_stride; i++)
      evdwl+=engv[i];
    if (_ef_atom) {
      if (_ilist==nullptr) {
        for (int i=0; i<_ev_stride; i++)
          eatom[i]+=engv[i];
      } else {
        for (int i=0; i<_ev_stride; i++)
          eatom[_ilist[i]]+=engv[i];
      }
    }
    vstart=_ev_stride;
  }
  if (_vflag) {
    int iend=vstart+_ev_stride;
    for (int j=0; j<6; j++) {
      for (int i=vstart; i<iend; i++)
        virial[j]+=engv[i];
      if (_vf_atom){
        if (_ilist==nullptr) {
          int ii=0;
          for (int i=vstart; i<iend; i++)
            vatom[ii++][j]+=engv[i];
        } else {
          int ii=0;
          for (int i=vstart; i<iend; i++)
            vatom[_ilist[ii++]][j]+=engv[i];
        }
      }
      vstart+=_ev_stride;
      iend+=_ev_stride;
    }
  }

  return evdwl;
}

template <class numtyp, class acctyp>
double AnswerT::energy_virial(double *eatom, double **vatom,
                              double *virial, double &ecoul) {
  if (_eflag==false && _vflag==false)
    return 0.0;

  if (_charge==false)
    return energy_virial(eatom,vatom,virial);

  double evdwl=0.0;
  int ii, vstart=0, iend=_ev_stride;
  if (_eflag) {
    iend=_ev_stride*2;
    #if (LAL_USE_OMP_SIMD == 1)
    #pragma omp simd reduction(+:evdwl)
    #endif
    for (int i=0; i<_ev_stride; i++)
      evdwl+=engv[i];
    double ecv=0.0;
    #if (LAL_USE_OMP_SIMD == 1)
    #pragma omp simd reduction(+:ecv)
    #endif
    for (int i=_ev_stride; i<iend; i++)
      ecv+=engv[i];
    ecoul+=ecv;
    if (_ef_atom) {
      if (_ilist==nullptr) {
        for (int i=0; i<_ev_stride; i++)
          eatom[i]+=engv[i];
        for (int i=_ev_stride; i<iend; i++)
          eatom[i]+=engv[i];
      } else {
        for (int i=0, ii=0; i<_ev_stride; i++)
          eatom[_ilist[ii++]]+=engv[i];
        for (int i=_ev_stride, ii=0; i<iend; i++)
          eatom[_ilist[ii++]]+=engv[i];
      }
    }
    vstart=iend;
    iend+=_ev_stride;
  }
  if (_vflag) {
    for (int j=0; j<6; j++) {
      for (int i=vstart; i<iend; i++)
        virial[j]+=engv[i];
      if (_vf_atom) {
        if (_ilist==nullptr) {
          for (int i=vstart, ii=0; i<iend; i++)
            vatom[ii++][j]+=engv[i];
        } else {
          for (int i=vstart, ii=0; i<iend; i++)
            vatom[_ilist[ii++]][j]+=engv[i];
        }
      }
      vstart+=_ev_stride;
      iend+=_ev_stride;
    }
  }

  return evdwl;
}

template <class numtyp, class acctyp>
void AnswerT::get_answers(double **f, double **tor) {
  if (_ilist==nullptr) {
    typedef struct { double x,y,z; } vec3d;
    typedef struct { acctyp x,y,z,w; } vec4d_t;
    vec3d *fp=reinterpret_cast<vec3d*>(&(f[0][0]));
    vec4d_t *forcep=reinterpret_cast<vec4d_t*>(&(force[0]));

    #if (LAL_USE_OMP == 1)
    #pragma omp parallel
    #endif
    {
      #if (LAL_USE_OMP == 1)
      const int nthreads = omp_get_num_threads();
      const int tid = omp_get_thread_num();
      const int idelta = _inum / nthreads + 1;
      const int ifrom = tid * idelta;
      const int ito = std::min(ifrom + idelta, _inum);
      #else
      const int tid = 0;
      const int ifrom = 0;
      const int ito = _inum;
      #endif

      for (int i=ifrom; i<ito; i++) {
        fp[i].x+=forcep[i].x;
        fp[i].y+=forcep[i].y;
        fp[i].z+=forcep[i].z;
      }
      if (_rot) {
        vec3d *torp=reinterpret_cast<vec3d*>(&(tor[0][0]));
        vec4d_t *torquep=reinterpret_cast<vec4d_t*>(&(force[_inum*4]));
        for (int i=ifrom; i<ito; i++) {
          torp[i].x+=torquep[i].x;
          torp[i].y+=torquep[i].y;
          torp[i].z+=torquep[i].z;
        }
      }
    }
  } else {
    #if (LAL_USE_OMP == 1)
    #pragma omp parallel
    #endif
    {
      #if (LAL_USE_OMP == 1)
      const int nthreads = omp_get_num_threads();
      const int tid = omp_get_thread_num();
      const int idelta = _inum / nthreads + 1;
      const int ifrom = tid * idelta;
      const int ito = std::min(ifrom + idelta, _inum);
      int fl=ifrom*4;
      #else
      const int tid = 0;
      const int ifrom = 0;
      const int ito = _inum;
      int fl=0;
      #endif

      for (int i=ifrom; i<ito; i++) {
        int ii=_ilist[i];
        f[ii][0]+=force[fl];
        f[ii][1]+=force[fl+1];
        f[ii][2]+=force[fl+2];
        fl+=4;
      }
      if (_rot) {
        fl=_inum*4 + ifrom*4;
        for (int i=ifrom; i<ito; i++) {
          int ii=_ilist[i];
          tor[ii][0]+=force[fl];
          tor[ii][1]+=force[fl+1];
          tor[ii][2]+=force[fl+2];
          fl+=4;
        }
      }
    }
  }
}

template <class numtyp, class acctyp>
void AnswerT::cq(const int cq_index) {
  engv.cq(dev->cq(cq_index));
  force.cq(dev->cq(cq_index));
  time_answer.clear();
  time_answer.init(*dev,dev->cq(cq_index));
  time_answer.zero();
}

template class Answer<PRECISION,ACC_PRECISION>;
}
