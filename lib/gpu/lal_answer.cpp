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
  
  int ans_elements=4;
  if (_rot)
    ans_elements+=4;
  
  // Ignore host/device transfers?
  bool cpuview=false;
  if (dev->device_type()==UCL_CPU)
    cpuview=true;
    
  // --------------------------   Host allocations
  success=success &&(host_ans.alloc(ans_elements*_max_local,*dev)==UCL_SUCCESS);
  success=success &&(host_engv.alloc(_ev_fields*_max_local,*dev)==UCL_SUCCESS);
    
  // ---------------------------  Device allocations
  if (cpuview) {
    dev_engv.view(host_engv);
    dev_ans.view(host_ans);
  } else {
    success=success && (dev_engv.alloc(_ev_fields*_max_local,*dev,
                                       UCL_WRITE_ONLY)==UCL_SUCCESS);
    success=success && (dev_ans.alloc(ans_elements*_max_local,
                                      *dev,UCL_WRITE_ONLY)==UCL_SUCCESS);
  }
  _gpu_bytes=dev_engv.row_bytes()+dev_ans.row_bytes();
  
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
    clear_resize();
    return alloc(inum);
  }
  return true;
}

template <class numtyp, class acctyp>
void AnswerT::clear_resize() {
  if (!_allocated)
    return;
  _allocated=false;

  dev_ans.clear();
  dev_engv.clear();
  host_ans.clear();
  host_engv.clear();
}

template <class numtyp, class acctyp>
void AnswerT::clear() {
  _gpu_bytes=0;
  if (!_allocated)
    return;

  time_answer.clear();
  clear_resize();
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
    ucl_copy(host_engv,dev_engv,_inum*csize,true);
  if (_rot)
    ucl_copy(host_ans,dev_ans,_inum*4*2,true);
  else
    ucl_copy(host_ans,dev_ans,_inum*4,true);
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
  double virial_acc[6];
  for (int i=0; i<6; i++) virial_acc[i]=0.0;
  if (_ilist==NULL) {
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
            virial_acc[j]+=*ap;
            ap+=_inum;
          }
        } else {
          for (int j=0; j<6; j++) {
            virial_acc[j]+=*ap;
            ap+=_inum;
          }
        }
      }
    }
    for (int j=0; j<6; j++)
      virial[j]+=virial_acc[j]*0.5;
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
            virial_acc[j]+=*ap;
            ap+=_inum;
          }
        } else {
          for (int j=0; j<6; j++) {
            virial_acc[j]+=*ap;
            ap+=_inum;
          }
        }
      }
    }
    for (int j=0; j<6; j++)
      virial[j]+=virial_acc[j]*0.5;
  }
  
  evdwl*=0.5;
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
  double _ecoul=0.0;
  double virial_acc[6];
  for (int i=0; i<6; i++) virial_acc[i]=0.0;
  if (_ilist==NULL) {
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
            virial_acc[j]+=*ap;
            ap+=_inum;
          }
        } else {
          for (int j=0; j<6; j++) {
            virial_acc[j]+=*ap;
            ap+=_inum;
          }
        }
      }
    }
    for (int j=0; j<6; j++)
      virial[j]+=virial_acc[j]*0.5;
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
            virial_acc[j]+=*ap;
            ap+=_inum;
          }
        } else {
          for (int j=0; j<6; j++) {
            virial_acc[j]+=*ap;
            ap+=_inum;
          }
        }
      }
    }
    for (int j=0; j<6; j++)
      virial[j]+=virial_acc[j]*0.5;
  }
  
  evdwl*=0.5;
  ecoul+=_ecoul*0.5;
  return evdwl;
}

template <class numtyp, class acctyp>
void AnswerT::get_answers(double **f, double **tor) {
  acctyp *ap=host_ans.begin();
  if (_ilist==NULL) {
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

template class Answer<PRECISION,ACC_PRECISION>;
