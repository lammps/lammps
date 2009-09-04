/***************************************************************************
                              lj_gpu_kernel.cu
                             -------------------
                               W. Michael Brown

  Routines that actually perform the force computation

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

#ifndef LJ_GPU_KERNEL
#define LJ_GPU_KERNEL

template<class numtyp, class acctyp>
__global__ void kernel_lj(const numtyp *special_lj, const int *dev_nbor, 
                          const int *dev_ij, const int nbor_pitch, acctyp *ans, 
                          size_t ans_pitch, const bool eflag, 
                          const bool vflag, const int inum, const int nall) {
  __shared__ numtyp sp_lj[4];
                            
  // ii indexes the two interacting particles in gi
  int ii=threadIdx.x;
  if (ii<4)
    sp_lj[ii]=special_lj[ii];    
  ii+=INT_MUL(blockIdx.x,blockDim.x);
  __syncthreads();

  if (ii<inum) {
  
    acctyp energy=(numtyp)0;
    acctyp fx=(numtyp)0;
    acctyp fy=(numtyp)0;
    acctyp fz=(numtyp)0;
    acctyp virial[6];
    for (int i=0; i<6; i++)
      virial[i]=(numtyp)0;
  
    const int *nbor=dev_nbor+ii;
    int i=*nbor;
    nbor+=nbor_pitch;
    int numj=*nbor;
    nbor+=nbor_pitch;
    const int *list=dev_ij+*nbor;
    const int *list_end=list+numj;
  
    numtyp ix=_x_<numtyp>(i,0);
    numtyp iy=_x_<numtyp>(i,1);
    numtyp iz=_x_<numtyp>(i,2);
    int itype=_x_<numtyp>(i,3);

    numtyp factor_lj;
    for ( ; list<list_end; list++) {
  
      int j=*list;
      if (j < nall) 
        factor_lj = 1.0;
      else {
        factor_lj = sp_lj[j/nall];
        j %= nall;
      }
      int jtype=_x_<numtyp>(j,3);

      // Compute r12
      numtyp delx = ix-_x_<numtyp>(j,0);
      numtyp dely = iy-_x_<numtyp>(j,1);
      numtyp delz = iz-_x_<numtyp>(j,2);
      numtyp r2inv = delx*delx+dely*dely+delz*delz;
        
      if (r2inv<_cutsq_<numtyp>(itype,jtype)) {
        r2inv=(numtyp)1.0/r2inv;
        numtyp r6inv =r2inv*r2inv*r2inv;
        numtyp force =factor_lj*r2inv*r6inv*(_lj1_<numtyp>(itype,jtype).x*r6inv-
                                             _lj1_<numtyp>(itype,jtype).y);
      
        fx+=delx*force;
        fy+=dely*force;
        fz+=delz*force;

        if (eflag) {
          numtyp e=r6inv*(_lj3_<numtyp>(itype,jtype).x*r6inv-
                          _lj3_<numtyp>(itype,jtype).y);
          energy+=factor_lj*(e-_offset_<numtyp>(1,1)); 
        }
        if (vflag) {
          virial[0] += delx*delx*force;
          virial[1] += dely*dely*force;
          virial[2] += delz*delz*force;
          virial[3] += delx*dely*force;
          virial[4] += delx*delz*force;
          virial[5] += dely*delz*force;
        }
      }

    } // for nbor

    // Store answers
    acctyp *ap1=ans+ii;
    if (eflag) {
      *ap1=energy;
      ap1+=ans_pitch;
    }
    if (vflag) {
      for (int i=0; i<6; i++) {
        *ap1=virial[i];
        ap1+=ans_pitch;
      }
    }
    *ap1=fx;
    ap1+=ans_pitch;
    *ap1=fy;
    ap1+=ans_pitch;
    *ap1=fz;

  } // if ii
}

template<class numtyp, class acctyp>
__global__ void kernel_lj_fast(const numtyp *special_lj, const int *dev_nbor, 
                               const int *dev_ij, const int nbor_pitch, 
                               acctyp *ans, size_t ans_pitch,const bool eflag, 
                               const bool vflag, const int inum, 
                               const int nall) {
                                
  // ii indexes the two interacting particles in gi
  int ii=threadIdx.x;
  __shared__ numtyp sp_lj[4];
  __shared__ numtyp cutsq[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __shared__ numtyp lj1[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __shared__ numtyp lj2[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __shared__ numtyp lj3[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __shared__ numtyp lj4[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __shared__ numtyp offset[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  if (ii<4)
    sp_lj[ii]=special_lj[ii];    
  if (ii<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
    int itype=ii/MAX_SHARED_TYPES;
    int jtype=ii%MAX_SHARED_TYPES;
    cutsq[ii]=_cutsq_<numtyp>(itype,jtype);
    lj1[ii]=_lj1_<numtyp>(itype,jtype).x;
    lj2[ii]=_lj1_<numtyp>(itype,jtype).y;
    if (eflag) {
      lj3[ii]=_lj3_<numtyp>(itype,jtype).x;
      lj4[ii]=_lj3_<numtyp>(itype,jtype).y;
      offset[ii]=_offset_<numtyp>(itype,jtype);
    }
  }
  ii+=INT_MUL(blockIdx.x,blockDim.x);
  __syncthreads();
  
  if (ii<inum) {
  
    acctyp energy=(numtyp)0;
    acctyp fx=(numtyp)0;
    acctyp fy=(numtyp)0;
    acctyp fz=(numtyp)0;
    acctyp virial[6];
    for (int i=0; i<6; i++)
      virial[i]=(numtyp)0;
  
    const int *nbor=dev_nbor+ii;
    int i=*nbor;
    nbor+=nbor_pitch;
    int numj=*nbor;
    nbor+=nbor_pitch;
    const int *list=dev_ij+*nbor;
    const int *list_end=list+numj;
  
    numtyp ix=_x_<numtyp>(i,0);
    numtyp iy=_x_<numtyp>(i,1);
    numtyp iz=_x_<numtyp>(i,2);
    int itype=INT_MUL(MAX_SHARED_TYPES,_x_<numtyp>(i,3));

    numtyp factor_lj;
    for ( ; list<list_end; list++) {
  
      int j=*list;
      if (j < nall) 
        factor_lj = 1.0;
      else {
        factor_lj = sp_lj[j/nall];
        j %= nall;
      }
      int mtype=itype+_x_<numtyp>(j,3);

      // Compute r12
      numtyp delx = ix-_x_<numtyp>(j,0);
      numtyp dely = iy-_x_<numtyp>(j,1);
      numtyp delz = iz-_x_<numtyp>(j,2);
      numtyp r2inv = delx*delx+dely*dely+delz*delz;

      if (r2inv<cutsq[mtype]) {
        r2inv=(numtyp)1.0/r2inv;
        numtyp r6inv = r2inv*r2inv*r2inv;
        numtyp force = factor_lj*r2inv*r6inv*(lj1[mtype]*r6inv-lj2[mtype]);
      
        fx+=delx*force;
        fy+=dely*force;
        fz+=delz*force;

        if (eflag) {
          numtyp e=r6inv*(lj3[mtype]*r6inv-lj4[mtype]);
          energy+=factor_lj*(e-offset[mtype]); 
        }
        if (vflag) {
          virial[0] += delx*delx*force;
          virial[1] += dely*dely*force;
          virial[2] += delz*delz*force;
          virial[3] += delx*dely*force;
          virial[4] += delx*delz*force;
          virial[5] += dely*delz*force;
        }
      }

    } // for nbor

    // Store answers
    acctyp *ap1=ans+ii;
    if (eflag) {
      *ap1=energy;
      ap1+=ans_pitch;
    }
    if (vflag) {
      for (int i=0; i<6; i++) {
        *ap1=virial[i];
        ap1+=ans_pitch;
      }
    }
    *ap1=fx;
    ap1+=ans_pitch;
    *ap1=fy;
    ap1+=ans_pitch;
    *ap1=fz;

  } // if ii
}

#endif
