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

using namespace LAMMPS_AL;
#define AnswerT Answer<numtyp,acctyp>

template <class numtyp, class acctyp>
AnswerT::Answer() : _allocated(false),_eflag(false),_vflag(false),
                            _inum(0),_ilist(NULL),_newton(false) {
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
  if (!_allocated)
    return;
  _allocated=false;

  force.clear();
  engv.clear();
  time_answer.clear();
  _inum=0;
  _ilist=NULL;
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
    engv.update_host(_inum*csize,true);
  if (_rot)
    force.update_host(_inum*4*2,true);
  else
    force.update_host(_inum*4,true);
  time_answer.stop();
}

template <class numtyp, class acctyp>
void AnswerT::copy_answers(const bool eflag, const bool vflag,
                               const bool ef_atom, const bool vf_atom,
                               int *ilist) {
  _ilist=ilist;
  copy_answers(eflag,vflag,ef_atom,vf_atom);
}

template <class numtyp, class acctyp>
double AnswerT::energy_virial(double *eatom, double **vatom,
                                  double *virial) {
  if (_eflag==false && _vflag==false)
    return 0.0;

  double evdwl=0.0;
  int vstart=0;
  if (_eflag) {
    for (int i=0; i<_inum; i++)
      evdwl+=engv[i];
    if (_ef_atom)
      if (_ilist==NULL) 
        for (int i=0; i<_inum; i++)
          eatom[i]+=engv[i];
      else
        for (int i=0; i<_inum; i++)
          eatom[_ilist[i]]+=engv[i];
    vstart=_inum;
  }
  if (_vflag) {
    int iend=vstart+_inum;
    for (int j=0; j<6; j++) {
      for (int i=vstart; i<iend; i++)
        virial[j]+=engv[i];
      if (_vf_atom)
        if (_ilist==NULL)
          for (int i=vstart; i<iend; i++)
            vatom[i][j]+=engv[i];
        else
          for (int i=vstart; i<iend; i++)
            vatom[_ilist[i]][j]+=engv[i];
      vstart+=_inum;
      iend+=_inum;
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
  int vstart=0, iend=_inum*2;
  if (_eflag) {
    for (int i=0; i<_inum; i++)
      evdwl+=engv[i];
    for (int i=_inum; i<iend; i++)
      ecoul+=engv[i];
    if (_ef_atom)
      if (_ilist==NULL) {
        for (int i=0; i<_inum; i++)
          eatom[i]+=engv[i];
        for (int i=_inum; i<iend; i++)
          eatom[i]+=engv[i];
      } else {
        for (int i=0; i<_inum; i++)
          eatom[_ilist[i]]+=engv[i];
        for (int i=_inum; i<iend; i++)
          eatom[_ilist[i]]+=engv[i];
      }
    vstart=iend;
    iend+=_inum;
  }
  if (_vflag) {
    for (int j=0; j<6; j++) {
      for (int i=vstart; i<iend; i++)
        virial[j]+=engv[i];
      if (_vf_atom)
        if (_ilist==NULL)
          for (int i=vstart; i<iend; i++)
              vatom[i][j]+=engv[i];
        else
          for (int i=vstart; i<iend; i++)
              vatom[_ilist[i]][j]+=engv[i];
      vstart+=_inum;
      iend+=_inum;
    }
  }
  
  return evdwl;
}

template <class numtyp, class acctyp>
void AnswerT::get_answers(double **f, double **tor) {
  int fl=0;
  if (_ilist==NULL) {
    for (int i=0; i<_inum; i++) {
      f[i][0]+=force[fl];
      f[i][1]+=force[fl+1];
      f[i][2]+=force[fl+2];
      fl+=4;
    }
    if (_rot) {
      for (int i=0; i<_inum; i++) {
        tor[i][0]+=force[fl];
        tor[i][1]+=force[fl+1];
        tor[i][2]+=force[fl+2];
        fl+=4;
      }
    }
  } else {
    for (int i=0; i<_inum; i++) {
      int ii=_ilist[i];
      f[ii][0]+=force[fl];
      f[ii][1]+=force[fl+1];
      f[ii][2]+=force[fl+2];
      fl+=4;
    }
    if (_rot) {
      for (int i=0; i<_inum; i++) {
        int ii=_ilist[i];
        tor[ii][0]+=force[fl];
        tor[ii][1]+=force[fl+1];
        tor[ii][2]+=force[fl+2];
        fl+=4;
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

