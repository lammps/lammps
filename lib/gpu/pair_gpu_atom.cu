/***************************************************************************
                               pair_gpu_atom.cu
                             -------------------
                               W. Michael Brown

  Memory routines for moving atom and force data between host and gpu

 __________________________________________________________________________
    This file is part of the LAMMPS GPU Library
 __________________________________________________________________________

    begin                : Tue Aug 4 2009
    copyright            : (C) 2009 by W. Michael Brown
    email                : wmbrown@sandia.gov
 ***************************************************************************/

/* -----------------------------------------------------------------------
   Copyright (2009) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.
   ----------------------------------------------------------------------- */

#include "pair_gpu_texture.h"
#include "pair_gpu_atom.h"

#define PairGPUAtomT PairGPUAtom<numtyp,acctyp>

template <class numtyp, class acctyp>
int PairGPUAtomT::bytes_per_atom() const { 
  return atom_fields()*sizeof(numtyp)+ans_fields()*sizeof(acctyp); 
}

template <class numtyp, class acctyp>
void PairGPUAtomT::init(const int max_atoms) {
  if (allocated)
    clear();
    
  _max_atoms=max_atoms;

  // Initialize timers for the selected GPU
  time_atom.init();
  time_answer.init();

  // Device matrices for atom and force data
  dev_x.safe_alloc(atom_fields(),max_atoms,x_get_texture<numtyp>());
  ans.safe_alloc(ans_fields(),max_atoms);

  // Get a host write only buffer
  host_write.safe_alloc_w(max_atoms*4);
  // Get a host read/write buffer
  host_read.safe_alloc_rw(ans.row_size()*ans_fields());
    
  allocated=true;
}
  
template <class numtyp, class acctyp>
void PairGPUAtomT::clear() {
  if (!allocated)
      return;
  allocated=false;
      
  dev_x.unbind();
  ans.clear();                               
  host_write.clear();
  host_read.clear();
  dev_x.clear();
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
      
  host_read.copy_from_device(ans.begin(),ans.row_size()*csize,s);
}
  
template <class numtyp, class acctyp>
double PairGPUAtomT::energy_virial(const int *ilist, const bool eflag_atom,
                                   const bool vflag_atom, double *eatom, 
                                   double **vatom, double *virial) {
  double evdwl=0.0;
  int gap=ans.row_size()-_inum;

  acctyp *ap=host_read.begin();
  if (_eflag) {
    if (eflag_atom) {
      for (int i=0; i<_inum; i++) {
        evdwl+=*ap;
        eatom[ilist[i]]+=*ap*0.5;
        ap++;
      }
    } else
      for (int i=0; i<_inum; i++) {
        evdwl+=*ap;
        ap++;
      }
    ap+=gap;
    evdwl*=0.5;
  }
  _read_loc=ap;
  gap=ans.row_size();
  if (_vflag) {
    if (vflag_atom) {
      for (int ii=0; ii<_inum; ii++) {
        int i=ilist[ii];
        ap=_read_loc+ii;
        for (int j=0; j<6; j++) {
          vatom[i][j]+=*ap*0.5;
          virial[j]+=*ap;
          ap+=gap;
        }
      }
    } else {
      for (int ii=0; ii<_inum; ii++) {
        ap=_read_loc+ii;
        for (int j=0; j<6; j++) {
          virial[j]+=*ap;
          ap+=gap;
        }
      }
    }
    for (int j=0; j<6; j++)
      virial[j]*=0.5;
    _read_loc+=gap*6;
  }
  
  return evdwl;
}

template <class numtyp, class acctyp>
void PairGPUAtomT::add_forces(const int *ilist, double **f) {
  int gap=ans.row_size();
  for (int ii=0; ii<_inum; ii++) {
    acctyp *ap=_read_loc+ii;
    int i=ilist[ii];
    f[i][0]+=*ap;
    ap+=gap;
    f[i][1]+=*ap;
    ap+=gap;
    f[i][2]+=*ap;
  }
}

template <class numtyp, class acctyp>
void PairGPUAtomT::add_torques(const int *ilist, double **tor, const int n) {
  int gap=ans.row_size();
  _read_loc+=gap*3;
  for (int ii=0; ii<n; ii++) {
    acctyp *ap=_read_loc+ii;
    int i=ilist[ii];
    tor[i][0]+=*ap;
    ap+=gap;
    tor[i][1]+=*ap;
    ap+=gap;
    tor[i][2]+=*ap;
  }
}

template class PairGPUAtom<PRECISION,ACC_PRECISION>;
