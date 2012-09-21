// **************************************************************************
//                                   lal_table.cu
//                             -------------------
//                           Trung Dac Nguyen (ORNL)
//
//  Device code for acceleration of the table pair style
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                : 
//    email                : nguyentd@ornl.gov
// ***************************************************************************/

#ifdef NV_KERNEL
#include "lal_aux_fun1.h"
#ifndef _DOUBLE_DOUBLE
texture<float4> pos_tex;
#else
texture<int4,1> pos_tex;
#endif
#else
#define pos_tex x_
#endif

#define LOOKUP 0
#define LINEAR 1
#define SPLINE 2
#define BITMAP 3

#ifndef __UNION_INT_FLOAT
#define __UNION_INT_FLOAT
typedef union {
  int i;
  float f;
} union_int_float;
#endif

/// ---------------- LOOKUP -------------------------------------------------

__kernel void k_table(const __global numtyp4 *restrict x_, 
                      const __global int *restrict tabindex,
                      const __global numtyp4 *restrict coeff2, 
                      const __global numtyp4 *restrict coeff3,
                      const __global numtyp4 *restrict coeff4,
                      const int lj_types,
                      const __global numtyp *restrict cutsq,
                      const __global numtyp *restrict sp_lj_in, 
                      const __global int *dev_nbor, 
                      const __global int *dev_packed, 
                      __global acctyp4 *restrict ans, 
                      __global acctyp *restrict engv, 
                      const int eflag, const int vflag, const int inum, 
                      const int nbor_pitch, const int t_per_atom, 
                      int tablength) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);
  
  __local numtyp sp_lj[4];
  sp_lj[0]=sp_lj_in[0];
  sp_lj[1]=sp_lj_in[1];
  sp_lj[2]=sp_lj_in[2];
  sp_lj[3]=sp_lj_in[3];
  
  acctyp energy=(acctyp)0;
  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(acctyp)0;
  
  int tlm1 = tablength - 1;
  
  if (ii<inum) {
    const __global int *nbor, *list_end;
    int i, numj, n_stride;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,list_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int itype=ix.w;
    
    numtyp factor_lj;
    for ( ; nbor<list_end; nbor+=n_stride) {
  
      int j=*nbor;
      factor_lj = sp_lj[sbmask(j)];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int mtype=itype*lj_types+jx.w;
      int tbindex = tabindex[mtype];
      
      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;
          
      if (rsq<cutsq[mtype]) {
        int itable=0,idx;
        numtyp force = (numtyp)0;
        itable = (rsq - coeff2[mtype].x) * coeff2[mtype].y;
        if (itable < tlm1) {
          idx = itable + tbindex*tablength;
          force = factor_lj * coeff3[idx].z;
        } else force = (numtyp)0.0;
                       
        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (eflag>0) {
          numtyp e = (numtyp)0.0;
          if (itable < tlm1) 
            e = coeff3[idx].y;
          energy+=factor_lj*e;
        }
        if (vflag>0) {
          virial[0] += delx*delx*force;
          virial[1] += dely*dely*force;
          virial[2] += delz*delz*force;
          virial[3] += delx*dely*force;
          virial[4] += delx*delz*force;
          virial[5] += dely*delz*force;
        }
      }

    } // for nbor
    store_answers(f,energy,virial,ii,inum,tid,t_per_atom,offset,eflag,vflag,
                  ans,engv);
  } // if ii
}

__kernel void k_table_fast(const __global numtyp4 *restrict x_,
                           const __global int *restrict tabindex,
                           const __global numtyp4 *restrict coeff2, 
                           const __global numtyp4 *restrict coeff3,
                           const __global numtyp4 *restrict coeff4,
                           const __global numtyp *restrict cutsq_in,
                           const __global numtyp *restrict sp_lj_in, 
                           const __global int *dev_nbor, 
                           const __global int *dev_packed, 
                           __global acctyp4 *restrict ans, 
                           __global acctyp *restrict engv, 
                           const int eflag, const int vflag, const int inum, 
                           const int nbor_pitch, const int t_per_atom, 
                           int tablength) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);
  
  __local numtyp cutsq[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp sp_lj[4];
  if (tid<4)
    sp_lj[tid]=sp_lj_in[tid];
  if (tid<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
    cutsq[tid]=cutsq_in[tid];
  }
  
  acctyp energy=(acctyp)0;
  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(acctyp)0;
  
  __syncthreads();
 
  int tlm1 = tablength - 1;
  
  if (ii<inum) {
    const __global int *nbor, *list_end;
    int i, numj, n_stride;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,list_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int iw=ix.w;
    int itype=fast_mul((int)MAX_SHARED_TYPES,iw);
    
    numtyp factor_lj;
    for ( ; nbor<list_end; nbor+=n_stride) {
  
      int j=*nbor;
      factor_lj = sp_lj[sbmask(j)];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int mtype=itype+jx.w;
      int tbindex = tabindex[mtype];
      
      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;
          
      if (rsq<cutsq[mtype]) {
        int itable=0,idx;
        numtyp force = (numtyp)0;
        itable = (rsq - coeff2[mtype].x) * coeff2[mtype].y;
        if (itable < tlm1) {
          idx = itable + tbindex*tablength;
          force = factor_lj * coeff3[idx].z;
        } else force = (numtyp)0.0;
                       
        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (eflag>0) {
          numtyp e = (numtyp)0.0;
          if (itable < tlm1) 
            e = coeff3[idx].y;
          energy+=factor_lj*e;
        }
        if (vflag>0) {
          virial[0] += delx*delx*force;
          virial[1] += dely*dely*force;
          virial[2] += delz*delz*force;
          virial[3] += delx*dely*force;
          virial[4] += delx*delz*force;
          virial[5] += dely*delz*force;
        }
      }

    } // for nbor
    store_answers(f,energy,virial,ii,inum,tid,t_per_atom,offset,eflag,vflag,
                  ans,engv);
  } // if ii
}

/// ---------------- LINEAR -------------------------------------------------

__kernel void k_table_linear(const __global numtyp4 *restrict x_, 
                             const __global int *restrict tabindex,
                             const __global numtyp4 *restrict coeff2, 
                             const __global numtyp4 *restrict coeff3,
                             const __global numtyp4 *restrict coeff4,
                             const int lj_types,
                             const __global numtyp *restrict cutsq,
                             const __global numtyp *restrict sp_lj_in, 
                             const __global int *dev_nbor, 
                             const __global int *dev_packed, 
                             __global acctyp4 *restrict ans, 
                             __global acctyp *restrict engv, 
                             const int eflag, const int vflag, const int inum, 
                             const int nbor_pitch, const int t_per_atom, 
                             int tablength) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);
  
  __local numtyp sp_lj[4];
  sp_lj[0]=sp_lj_in[0];
  sp_lj[1]=sp_lj_in[1];
  sp_lj[2]=sp_lj_in[2];
  sp_lj[3]=sp_lj_in[3];

  acctyp energy=(acctyp)0;
  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(acctyp)0;
  
  int tlm1 = tablength - 1;
  
  if (ii<inum) {
    const __global int *nbor, *list_end;
    int i, numj, n_stride;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,list_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int itype=ix.w;
    
    numtyp factor_lj;
    for ( ; nbor<list_end; nbor+=n_stride) {
  
      int j=*nbor;
      factor_lj = sp_lj[sbmask(j)];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int mtype=itype*lj_types+jx.w;
      int tbindex = tabindex[mtype];
      
      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;
          
      if (rsq<cutsq[mtype]) {
        int itable=0,idx;
        numtyp fraction=(numtyp)0;
        numtyp value = (numtyp)0;
        numtyp force = (numtyp)0;
        itable = (rsq - coeff2[mtype].x) * coeff2[mtype].y;
        if (itable < tlm1) {
          idx = itable + tbindex*tablength;
          fraction = (rsq - coeff3[idx].x) * coeff2[mtype].y;
          value = coeff3[idx].z + fraction*coeff4[idx].z;
          force = factor_lj * value;
        } else force = (numtyp)0.0;
             
        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (eflag>0) {
          numtyp e = (numtyp)0.0;
          if (itable < tlm1) 
            e = coeff3[idx].y + fraction*coeff4[idx].y;
          energy+=factor_lj*e;
        }
        if (vflag>0) {
          virial[0] += delx*delx*force;
          virial[1] += dely*dely*force;
          virial[2] += delz*delz*force;
          virial[3] += delx*dely*force;
          virial[4] += delx*delz*force;
          virial[5] += dely*delz*force;
        }
      }

    } // for nbor
    store_answers(f,energy,virial,ii,inum,tid,t_per_atom,offset,eflag,vflag,
                  ans,engv);
  } // if ii
}

__kernel void k_table_linear_fast(const __global numtyp4 *restrict x_, 
                                  const __global int *restrict tabindex,
                                  const __global numtyp4 *restrict coeff2, 
                                  const __global numtyp4 *restrict coeff3,
                                  const __global numtyp4 *restrict coeff4,
                                  const __global numtyp *restrict cutsq_in,
                                  const __global numtyp *restrict sp_lj_in, 
                                  const __global int *dev_nbor, 
                                  const __global int *dev_packed, 
                                  __global acctyp4 *restrict ans, 
                                  __global acctyp *restrict engv, 
                                  const int eflag, const int vflag, 
                                  const int inum, const int nbor_pitch, 
                                  const int t_per_atom, int tablength) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);
  
  __local numtyp cutsq[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp sp_lj[4];
  if (tid<4)
    sp_lj[tid]=sp_lj_in[tid];
  if (tid<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
    cutsq[tid]=cutsq_in[tid];
  }
  
  acctyp energy=(acctyp)0;
  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(acctyp)0;

  __syncthreads();

  int tlm1 = tablength - 1;
  
  if (ii<inum) {
    const __global int *nbor, *list_end;
    int i, numj, n_stride;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,list_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int iw=ix.w;
    int itype=fast_mul((int)MAX_SHARED_TYPES,iw);
    
    numtyp factor_lj;
    for ( ; nbor<list_end; nbor+=n_stride) {
  
      int j=*nbor;
      factor_lj = sp_lj[sbmask(j)];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int mtype=itype+jx.w;
      int tbindex = tabindex[mtype];
      
      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;
          
      if (rsq<cutsq[mtype]) {
        int itable=0,idx;
        numtyp fraction=(numtyp)0;
        numtyp value = (numtyp)0;
        numtyp force = (numtyp)0;
        itable = (rsq - coeff2[mtype].x) * coeff2[mtype].y;
        if (itable < tlm1) {
          idx = itable + tbindex*tablength;
          fraction = (rsq - coeff3[idx].x) * coeff2[mtype].y;
          value = coeff3[idx].z + fraction*coeff4[idx].z;
          force = factor_lj * value;
        } else force = (numtyp)0.0;
             
        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (eflag>0) {
          numtyp e = (numtyp)0.0;
          if (itable < tlm1) 
            e = coeff3[idx].y + fraction*coeff4[idx].y;
          energy+=factor_lj*e;
        }
        if (vflag>0) {
          virial[0] += delx*delx*force;
          virial[1] += dely*dely*force;
          virial[2] += delz*delz*force;
          virial[3] += delx*dely*force;
          virial[4] += delx*delz*force;
          virial[5] += dely*delz*force;
        }
      }

    } // for nbor
    store_answers(f,energy,virial,ii,inum,tid,t_per_atom,offset,eflag,vflag,
                  ans,engv);
  } // if ii
}

/// ---------------- SPLINE -------------------------------------------------

__kernel void k_table_spline(const __global numtyp4 *restrict x_, 
                             const __global int *restrict tabindex,
                             const __global numtyp4 *restrict coeff2, 
                             const __global numtyp4 *restrict coeff3,
                             const __global numtyp4 *restrict coeff4,
                             const int lj_types,
                             const __global numtyp *restrict cutsq,
                             const __global numtyp *restrict sp_lj_in, 
                             const __global int *dev_nbor, 
                             const __global int *dev_packed, 
                             __global acctyp4 *restrict ans, 
                             __global acctyp *restrict engv, 
                             const int eflag, const int vflag, const int inum, 
                             const int nbor_pitch, const int t_per_atom, 
                             int tablength) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);
  
  __local numtyp sp_lj[4];
  sp_lj[0]=sp_lj_in[0];
  sp_lj[1]=sp_lj_in[1];
  sp_lj[2]=sp_lj_in[2];
  sp_lj[3]=sp_lj_in[3];
  
  acctyp energy=(acctyp)0;
  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(acctyp)0;
    
  int tlm1 = tablength - 1;
  
  if (ii<inum) {
    const __global int *nbor, *list_end;
    int i, numj, n_stride;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,list_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int itype=ix.w;
    
    numtyp factor_lj;
    for ( ; nbor<list_end; nbor+=n_stride) {
  
      int j=*nbor;
      factor_lj = sp_lj[sbmask(j)];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int mtype=itype*lj_types+jx.w;
      int tbindex = tabindex[mtype];
      
      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;
          
      if (rsq<cutsq[mtype]) {
        int itable=0,idx;
        numtyp a = (numtyp)0;
        numtyp b = (numtyp)0;
        numtyp value = (numtyp)0;
        numtyp force = (numtyp)0;
        itable = (rsq - coeff2[mtype].x) * coeff2[mtype].y;
        if (itable < tlm1) {
          idx = itable + tbindex*tablength;
          b = (rsq - coeff3[idx].x) * coeff2[mtype].y;
          a = (numtyp)1.0 - b;
          value = a * coeff3[idx].z + b * coeff3[idx+1].z + 
            ((a*a*a-a)*coeff4[idx].z + (b*b*b-b)*coeff4[idx+1].z) * 
                  coeff2[mtype].z;
          force = factor_lj * value;
        } else force = (numtyp)0.0;
              
        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (eflag>0) {
          numtyp e = (numtyp)0.0;
          if (itable < tlm1) {
            e = a * coeff3[idx].y + b * coeff3[idx+1].y + 
                ((a*a*a-a)*coeff4[idx].y + (b*b*b-b)*coeff4[idx+1].y) * 
                  coeff2[mtype].z;
          }  
          energy+=factor_lj*e;
        }
        if (vflag>0) {
          virial[0] += delx*delx*force;
          virial[1] += dely*dely*force;
          virial[2] += delz*delz*force;
          virial[3] += delx*dely*force;
          virial[4] += delx*delz*force;
          virial[5] += dely*delz*force;
        }
      }

    } // for nbor
    store_answers(f,energy,virial,ii,inum,tid,t_per_atom,offset,eflag,vflag,
                  ans,engv);
  } // if ii
}

__kernel void k_table_spline_fast(const __global numtyp4 *x_, 
                                  const __global int *tabindex,
                                  const __global numtyp4* coeff2, 
                                  const __global numtyp4 *coeff3,
                                  const __global numtyp4 *coeff4,
                                  const __global numtyp *cutsq_in,
                                  const __global numtyp* sp_lj_in, 
                                  const __global int *dev_nbor, 
                                  const __global int *dev_packed, 
                                  __global acctyp4 *ans, 
                                  __global acctyp *engv, 
                                  const int eflag, const int vflag, 
                                  const int inum, const int nbor_pitch, 
                                  const int t_per_atom, int tablength) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);
  
  __local numtyp cutsq[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp sp_lj[4];
  if (tid<4)
    sp_lj[tid]=sp_lj_in[tid];
  if (tid<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
    cutsq[tid]=cutsq_in[tid];
  }
  
  acctyp energy=(acctyp)0;
  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(acctyp)0;

  __syncthreads();
  
  int tlm1 = tablength - 1;
  
  if (ii<inum) {
    const __global int *nbor, *list_end;
    int i, numj, n_stride;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,list_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int iw=ix.w;
    int itype=fast_mul((int)MAX_SHARED_TYPES,iw);
    
    numtyp factor_lj;
    for ( ; nbor<list_end; nbor+=n_stride) {
  
      int j=*nbor;
      factor_lj = sp_lj[sbmask(j)];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int mtype=itype+jx.w;
      int tbindex = tabindex[mtype];
      
      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;
          
      if (rsq<cutsq[mtype]) {
        int itable=0,idx;
        numtyp a = (numtyp)0;
        numtyp b = (numtyp)0;
        numtyp value = (numtyp)0;
        numtyp force = (numtyp)0;
        itable = (rsq - coeff2[mtype].x) * coeff2[mtype].y;
        if (itable < tlm1) {
          idx = itable + tbindex*tablength;
          b = (rsq - coeff3[idx].x) * coeff2[mtype].y;
          a = (numtyp)1.0 - b;
          value = a * coeff3[idx].z + b * coeff3[idx+1].z + 
            ((a*a*a-a)*coeff4[idx].z + (b*b*b-b)*coeff4[idx+1].z) * 
                  coeff2[mtype].z;
          force = factor_lj * value;
        } else force = (numtyp)0.0;
              
        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (eflag>0) {
          numtyp e = (numtyp)0.0;
          if (itable < tlm1) {
            e = a * coeff3[idx].y + b * coeff3[idx+1].y + 
                ((a*a*a-a)*coeff4[idx].y + (b*b*b-b)*coeff4[idx+1].y) * 
                  coeff2[mtype].z;
          }  
          energy+=factor_lj*e;
        }
        if (vflag>0) {
          virial[0] += delx*delx*force;
          virial[1] += dely*dely*force;
          virial[2] += delz*delz*force;
          virial[3] += delx*dely*force;
          virial[4] += delx*delz*force;
          virial[5] += dely*delz*force;
        }
      }

    } // for nbor
    store_answers(f,energy,virial,ii,inum,tid,t_per_atom,offset,eflag,vflag,
                  ans,engv);
  } // if ii
}

/// ---------------- BITMAP -------------------------------------------------

__kernel void k_table_bitmap(const __global numtyp4 *x_, 
                             const __global int *tabindex,
                             const __global int *nshiftbits, 
                             const __global int *nmask,
                             const __global numtyp4* coeff2, 
                             const __global numtyp4 *coeff3,
                             const __global numtyp4 *coeff4,
                             const int lj_types,
                             const __global numtyp *cutsq,
                             const __global numtyp* sp_lj_in, 
                             const __global int *dev_nbor, 
                             const __global int *dev_packed, 
                             __global acctyp4 *ans, 
                             __global acctyp *engv, 
                             const int eflag, const int vflag, const int inum, 
                             const int nbor_pitch, const int t_per_atom, 
                             int tablength) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);
  
  __local numtyp sp_lj[4];
  sp_lj[0]=sp_lj_in[0];
  sp_lj[1]=sp_lj_in[1];
  sp_lj[2]=sp_lj_in[2];
  sp_lj[3]=sp_lj_in[3];
  
  acctyp energy=(acctyp)0;
  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(acctyp)0;
  
  int tlm1 = tablength - 1;
  
  if (ii<inum) {
    const __global int *nbor, *list_end;
    int i, numj, n_stride;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,list_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int itype=ix.w;
    
    numtyp factor_lj;
    for ( ; nbor<list_end; nbor+=n_stride) {
  
      int j=*nbor;
      factor_lj = sp_lj[sbmask(j)];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int mtype=itype*lj_types+jx.w;
      int tbindex = tabindex[mtype];
      
      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;
          
      if (rsq<cutsq[mtype]) {
        int itable=0,idx;
        numtyp fraction=(numtyp)0;
        numtyp value = (numtyp)0;
        numtyp force = (numtyp)0;
        union_int_float rsq_lookup; 
        rsq_lookup.f = rsq;
        itable = rsq_lookup.i & nmask[mtype];
        itable >>= nshiftbits[mtype];
        if (itable <= tlm1) {
          idx = itable + tbindex*tablength;
          fraction = (rsq_lookup.f - coeff3[idx].x) * coeff4[idx].w;
          value = coeff3[idx].z + fraction*coeff4[idx].z;
          force = factor_lj * value;
        } else force = (numtyp)0.0;
          
        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (eflag>0) {
          numtyp e = (numtyp)0.0;
          if (itable <= tlm1) 
            e = coeff3[idx].y + fraction*coeff4[idx].y;
          energy+=factor_lj*e;
        }
        if (vflag>0) {
          virial[0] += delx*delx*force;
          virial[1] += dely*dely*force;
          virial[2] += delz*delz*force;
          virial[3] += delx*dely*force;
          virial[4] += delx*delz*force;
          virial[5] += dely*delz*force;
        }
      }

    } // for nbor
    store_answers(f,energy,virial,ii,inum,tid,t_per_atom,offset,eflag,vflag,
                  ans,engv);
  } // if ii
}

__kernel void k_table_bitmap_fast(const __global numtyp4 *x_, 
                                  const __global int *tabindex,
                                  const __global int *nshiftbits, 
                                  const __global int *nmask,
                                  const __global numtyp4* coeff2, 
                                  const __global numtyp4 *coeff3,
                                  const __global numtyp4 *coeff4,
                                  const __global numtyp *cutsq_in,
                                  const __global numtyp* sp_lj_in, 
                                  const __global int *dev_nbor, 
                                  const __global int *dev_packed, 
                                  __global acctyp4 *ans, 
                                  __global acctyp *engv, 
                                  const int eflag, const int vflag, 
                                  const int inum, const int nbor_pitch, 
                                  const int t_per_atom, int tablength) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);
  
  __local numtyp cutsq[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp sp_lj[4];
  if (tid<4)
    sp_lj[tid]=sp_lj_in[tid];
  if (tid<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
    cutsq[tid]=cutsq_in[tid];
  }
  
  acctyp energy=(acctyp)0;
  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(acctyp)0;
  
  __syncthreads();
  
  int tlm1 = tablength - 1;
  
  if (ii<inum) {
    const __global int *nbor, *list_end;
    int i, numj, n_stride;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,list_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int iw=ix.w;
    int itype=fast_mul((int)MAX_SHARED_TYPES,iw);
    
    numtyp factor_lj;
    for ( ; nbor<list_end; nbor+=n_stride) {
  
      int j=*nbor;
      factor_lj = sp_lj[sbmask(j)];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int mtype=itype+jx.w;
      int tbindex = tabindex[mtype];
      
      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;
          
      if (rsq<cutsq[mtype]) {
        int itable=0,idx;
        numtyp fraction=(numtyp)0;
        numtyp value = (numtyp)0;
        numtyp force = (numtyp)0;
        union_int_float rsq_lookup; 
        rsq_lookup.f = rsq;
        itable = rsq_lookup.i & nmask[mtype];
        itable >>= nshiftbits[mtype];
        if (itable <= tlm1) {
          idx = itable + tbindex*tablength;
          fraction = (rsq_lookup.f - coeff3[idx].x) * coeff4[idx].w;
          value = coeff3[idx].z + fraction*coeff4[idx].z;
          force = factor_lj * value;
        } else force = (numtyp)0.0;
          
        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (eflag>0) {
          numtyp e = (numtyp)0.0;
          if (itable <= tlm1) 
            e = coeff3[idx].y + fraction*coeff4[idx].y;
          energy+=factor_lj*e;
        }
        if (vflag>0) {
          virial[0] += delx*delx*force;
          virial[1] += dely*dely*force;
          virial[2] += delz*delz*force;
          virial[3] += delx*dely*force;
          virial[4] += delx*delz*force;
          virial[5] += dely*delz*force;
        }
      }

    } // for nbor
    store_answers(f,energy,virial,ii,inum,tid,t_per_atom,offset,eflag,vflag,
                  ans,engv);
  } // if ii
}
