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
   Contributing authors: Mike Brown (SNL), wmbrown@sandia.gov
                         Peng Wang (Nvidia), penwang@nvidia.com
                         Paul Crozier (SNL), pscrozi@sandia.gov
------------------------------------------------------------------------- */

#include "pair_gpu_texture.h"
#include "pair_gpu_atom.h"

#define PairGPUAtomT PairGPUAtom<numtyp,acctyp>

template <class numtyp, class acctyp>
int PairGPUAtomT::bytes_per_atom() const { 
  return atom_fields()*sizeof(numtyp)+ans_fields()*sizeof(acctyp); 
}

template <class numtyp, class acctyp>
bool PairGPUAtomT::init(const int max_atoms) {
  bool success=true;
  
  if (allocated)
    clear();
    
  _max_atoms=max_atoms;

  // Initialize timers for the selected GPU
  time_atom.init();
  time_answer.init();

  // Device matrices for atom and force data
  success=success && dev_x.alloc(max_atoms*sizeof(vec4));
  success=success && dev_q.alloc(max_atoms*sizeof(vec4));
  success=success && ans.alloc(ans_fields()*max_atoms);
  // Get a host read/write buffer
  success=success && host_read.alloc_rw(max_atoms*ans_fields());

  // Get a host write only buffer
  success=success && host_write.alloc_w(max_atoms*atom_fields());
    
  allocated=true;
  
  return success;
}
  
template <class numtyp, class acctyp>
void PairGPUAtomT::resize(const int max_atoms, bool &success) {
  ans.clear();
  dev_x.clear();
  dev_q.clear();
  host_write.clear();
  host_read.clear();

  _max_atoms=max_atoms;

  success = success && dev_x.alloc(_max_atoms*sizeof(vec4));
  success = success && dev_q.alloc(_max_atoms*sizeof(vec4));
  success = success && ans.alloc(ans_fields()*_max_atoms);
  success = success && host_read.alloc_rw(_max_atoms*ans_fields());
  success = success && host_write.alloc_w(_max_atoms*atom_fields());
}  
 
template <class numtyp, class acctyp>
void PairGPUAtomT::clear() {
  if (!allocated)
      return;
  allocated=false;
      
  ans.clear();
  dev_x.clear();
  dev_q.clear();
  host_write.clear();
  host_read.clear();
}  
 
template <class numtyp, class acctyp>
double PairGPUAtomT::host_memory_usage(const int max_atoms) const {
  return max_atoms*atom_fields()*sizeof(numtyp)+
         ans_fields()*(max_atoms)*sizeof(acctyp)+
         sizeof(PairGPUAtom<numtyp,acctyp>);
}
  
template <class numtyp, class acctyp>
void PairGPUAtomT::copy_answers(const bool eflag, const bool vflag, 
                                cudaStream_t &s) {
  _eflag=eflag;
  _vflag=vflag;
    
  int csize=ans_fields();    
  if (!eflag)
    csize--;
  if (!vflag)
    csize-=6;
      
  host_read.copy_from_device(ans.begin(),_inum*csize,s);
}
  
template <class numtyp, class acctyp>
double PairGPUAtomT::energy_virial(const int *ilist, const bool eflag_atom,
                                   const bool vflag_atom, double *eatom, 
                                   double **vatom, double *virial,
                                   double **f, double **tor, const int n) {
  double evdwl=0.0;

  acctyp *ap=host_read.begin();
  for (int i=0; i<_inum; i++) {
    int ii=ilist[i];
    if (_eflag) {
      if (eflag_atom) {
        evdwl+=*ap;
        eatom[ii]+=*ap*0.5;
        ap++;
      } else {
        evdwl+=*ap;
        ap++;
      }
    }
    if (_vflag) {
      if (vflag_atom) {
        for (int j=0; j<6; j++) {
          vatom[ii][j]+=*ap*0.5;
          virial[j]+=*ap;
          ap++;
        }
      } else {
        for (int j=0; j<6; j++) {
          virial[j]+=*ap;
          ap++;
        }
      }
    }
    f[ii][0]+=*ap;
    ap++;
    f[ii][1]+=*ap;
    ap++;
    f[ii][2]+=*ap;
    ap++;
    if (i<n) {
      tor[ii][0]+=*ap;
      ap++;
      tor[ii][1]+=*ap;
      ap++;
      tor[ii][2]+=*ap;
      ap++;
    } else {
      ap+=3;
    }
  }
  for (int j=0; j<6; j++)
    virial[j]*=0.5;
  evdwl*=0.5;
  return evdwl;
}

template <class numtyp, class acctyp>
void PairGPUAtomT::copy_asphere(const int *ilist, double **f, double **tor, 
                                const int n) {
  acctyp *ap=host_read.begin();
  for (int i=0; i<_inum; i++) {
    int ii=ilist[i];
    f[ii][0]+=*ap;
    ap++;
    f[ii][1]+=*ap;
    ap++;
    f[ii][2]+=*ap;
    ap++;
    if (i<n) {
      tor[ii][0]+=*ap;
      ap++;
      tor[ii][1]+=*ap;
      ap++;
      tor[ii][2]+=*ap;
      ap++;
    } else {
      ap+=3;
    }
  }
}

template class PairGPUAtom<PRECISION,ACC_PRECISION>;
