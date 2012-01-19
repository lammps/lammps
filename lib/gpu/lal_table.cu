// **************************************************************************
//                                   lal_table.cu
//                             -------------------
//                      Trung Dac Nguyen, W. Michael Brown (ORNL)
//
//  Device code for acceleration of the table pair style
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                : 
//    email                : brownw@ornl.gov nguyentd@ornl.gov
// ***************************************************************************/

#ifdef NV_KERNEL
#include "lal_aux_fun1.h"
texture<float4> pos_tex;
#ifndef _DOUBLE_DOUBLE
ucl_inline float4 fetch_pos(const int& i, const float4 *pos) 
  { return tex1Dfetch(pos_tex, i); }
#endif
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

/// ------ -------inline functions per table style ----------------------------

ucl_inline void lookup(__global numtyp4 *x_, __global int *tabindex,
                     __global numtyp4* coeff2, 
                     __global numtyp4 *coeff3,
                     __global numtyp4 *coeff4,
                     const int lj_types,
                     __global numtyp *cutsq_in,
                     __global numtyp* sp_lj_in, 
                     __global int *dev_nbor, __global int *dev_packed, 
                     __global acctyp4 *ans, __global acctyp *engv, 
                     const int eflag, const int vflag, const int inum, 
                     const int nbor_pitch, const int t_per_atom, 
                     int tablength, bool shared_types) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);
  
  acctyp energy=(acctyp)0;
  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(acctyp)0;
  
  __local numtyp cutsq[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp sp_lj[4];
  __local numtyp *_cutsq, *_sp_lj; 
   
  int type_pitch = lj_types;  
  if (shared_types) {
    if (tid<4)
      sp_lj[tid]=sp_lj_in[tid];
    if (tid<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
      cutsq[tid]=cutsq_in[tid];
    }
    type_pitch =  MAX_SHARED_TYPES;
    _cutsq = cutsq;
    _sp_lj = sp_lj;
    __syncthreads();
  } else {
    sp_lj[0]=sp_lj_in[0];
    sp_lj[1]=sp_lj_in[1];
    sp_lj[2]=sp_lj_in[2];
    sp_lj[3]=sp_lj_in[3];
  
    _cutsq = cutsq_in;
    _sp_lj = sp_lj_in;
  } 
  
  int tlm1 = tablength - 1;
  
  if (ii<inum) {
    __global int *nbor, *list_end;
    int i, numj, n_stride;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,list_end,nbor);

    numtyp4 ix=fetch_pos(i,x_); //x_[i];
    int iw=ix.w;
    int itype=fast_mul((int)type_pitch,iw);
    
    numtyp factor_lj;
    for ( ; nbor<list_end; nbor+=n_stride) {
  
      int j=*nbor;
      factor_lj = _sp_lj[sbmask(j)];
      j &= NEIGHMASK;

      numtyp4 jx=fetch_pos(j,x_); //x_[j];
      int mtype=itype+jx.w;
      int tbindex = tabindex[mtype];
      
      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;
          
      if (rsq<_cutsq[mtype]) {
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

ucl_inline void linear(__global numtyp4 *x_, __global int *tabindex,
                               __global numtyp4* coeff2, 
                               __global numtyp4 *coeff3,
                               __global numtyp4 *coeff4,
                               const int lj_types,
                               __global numtyp *cutsq_in,
                               __global numtyp* sp_lj_in, 
                               __global int *dev_nbor, __global int *dev_packed, 
                               __global acctyp4 *ans, __global acctyp *engv, 
                               const int eflag, const int vflag, const int inum, 
                               const int nbor_pitch, const int t_per_atom, 
                               int tablength, bool shared_types) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);
  
  acctyp energy=(acctyp)0;
  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(acctyp)0;
  
  __local numtyp cutsq[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp sp_lj[4];
  __local numtyp *_cutsq, *_sp_lj;
        
  int type_pitch = lj_types;  
  if (shared_types) {
    if (tid<4)
      sp_lj[tid]=sp_lj_in[tid];
    if (tid<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
      cutsq[tid]=cutsq_in[tid];
    }
    type_pitch =  MAX_SHARED_TYPES;
    _cutsq = cutsq;
    _sp_lj = sp_lj;
    __syncthreads();
  } else {
    sp_lj[0]=sp_lj_in[0];
    sp_lj[1]=sp_lj_in[1];
    sp_lj[2]=sp_lj_in[2];
    sp_lj[3]=sp_lj_in[3];
  
    _cutsq = cutsq_in;
    _sp_lj = sp_lj_in;
  } 

  int tlm1 = tablength - 1;
  
  if (ii<inum) {
    __global int *nbor, *list_end;
    int i, numj, n_stride;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,list_end,nbor);

    numtyp4 ix=fetch_pos(i,x_); //x_[i];
    int iw=ix.w;
    int itype=fast_mul(type_pitch,iw);
    
    numtyp factor_lj;
    for ( ; nbor<list_end; nbor+=n_stride) {
  
      int j=*nbor;
      factor_lj = _sp_lj[sbmask(j)];
      j &= NEIGHMASK;

      numtyp4 jx=fetch_pos(j,x_); //x_[j];
      int mtype=itype+jx.w;
      int tbindex = tabindex[mtype];
      
      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;
          
      if (rsq<_cutsq[mtype]) {
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

ucl_inline void spline(__global numtyp4 *x_, __global int *tabindex,
                     __global numtyp4* coeff2, 
                     __global numtyp4 *coeff3,
                     __global numtyp4 *coeff4,
                     const int lj_types,
                     __global numtyp *cutsq_in,
                     __global numtyp* sp_lj_in, 
                     __global int *dev_nbor, __global int *dev_packed, 
                     __global acctyp4 *ans, __global acctyp *engv, 
                     const int eflag, const int vflag, const int inum, 
                     const int nbor_pitch, const int t_per_atom, 
                     int tablength, bool shared_types) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);
  
  acctyp energy=(acctyp)0;
  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(acctyp)0;
  
  __local numtyp cutsq[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp sp_lj[4];
  __local numtyp *_cutsq, *_sp_lj;
    
  int type_pitch = lj_types;  
  if (shared_types) {
    if (tid<4)
      sp_lj[tid]=sp_lj_in[tid];
    if (tid<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
      cutsq[tid]=cutsq_in[tid];
    }
    type_pitch =  MAX_SHARED_TYPES;
    _cutsq = cutsq;
    _sp_lj = sp_lj;
    __syncthreads();
  } else {
    sp_lj[0]=sp_lj_in[0];
    sp_lj[1]=sp_lj_in[1];
    sp_lj[2]=sp_lj_in[2];
    sp_lj[3]=sp_lj_in[3];
  
    _cutsq = cutsq_in;
    _sp_lj = sp_lj_in;
  } 
  
  int tlm1 = tablength - 1;
  
  if (ii<inum) {
    __global int *nbor, *list_end;
    int i, numj, n_stride;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,list_end,nbor);

    numtyp4 ix=fetch_pos(i,x_); //x_[i];
    int iw=ix.w;
    int itype=fast_mul(type_pitch,iw);
    
    numtyp factor_lj;
    for ( ; nbor<list_end; nbor+=n_stride) {
  
      int j=*nbor;
      factor_lj = _sp_lj[sbmask(j)];
      j &= NEIGHMASK;

      numtyp4 jx=fetch_pos(j,x_); //x_[j];
      int mtype=itype+jx.w;
      int tbindex = tabindex[mtype];
      
      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;
          
      if (rsq<_cutsq[mtype]) {
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

ucl_inline void bitmap(__global numtyp4 *x_, __global int *tabindex,
                       __global int *nshiftbits, __global int *nmask,
                       __global numtyp4* coeff2, 
                       __global numtyp4 *coeff3,
                       __global numtyp4 *coeff4,
                       const int lj_types,
                       __global numtyp *cutsq_in,
                       __global numtyp* sp_lj_in, 
                       __global int *dev_nbor, __global int *dev_packed, 
                       __global acctyp4 *ans, __global acctyp *engv, 
                       const int eflag, const int vflag, const int inum, 
                       const int nbor_pitch, const int t_per_atom, 
                       int tablength, bool shared_types) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);
  
  acctyp energy=(acctyp)0;
  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(acctyp)0;
  
  __local numtyp cutsq[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp sp_lj[4];
  __local numtyp *_cutsq, *_sp_lj;
  
  int type_pitch = lj_types;  
  if (shared_types) {
    if (tid<4)
      sp_lj[tid]=sp_lj_in[tid];
    if (tid<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
      cutsq[tid]=cutsq_in[tid];
    }
    type_pitch =  MAX_SHARED_TYPES;
    _cutsq = cutsq;
    _sp_lj = sp_lj;
    __syncthreads();
  } else {
    sp_lj[0]=sp_lj_in[0];
    sp_lj[1]=sp_lj_in[1];
    sp_lj[2]=sp_lj_in[2];
    sp_lj[3]=sp_lj_in[3];
   
    _cutsq = cutsq_in;
    _sp_lj = sp_lj_in;
  } 
  
  int tlm1 = tablength - 1;
  
  if (ii<inum) {
    __global int *nbor, *list_end;
    int i, numj, n_stride;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,list_end,nbor);

    numtyp4 ix=fetch_pos(i,x_); //x_[i];
    int iw=ix.w;
    int itype=fast_mul(type_pitch,iw);
    
    numtyp factor_lj;
    for ( ; nbor<list_end; nbor+=n_stride) {
  
      int j=*nbor;
      factor_lj = _sp_lj[sbmask(j)];
      j &= NEIGHMASK;

      numtyp4 jx=fetch_pos(j,x_); //x_[j];
      int mtype=itype+jx.w;
      int tbindex = tabindex[mtype];
      
      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;
          
      if (rsq<_cutsq[mtype]) {
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

/// -------------------------------------------------------------------------

__kernel void kernel_pair(__global numtyp4 *x_, __global int *tabindex,
                          __global int *nshiftbits, __global int *nmask,
                          __global numtyp4* coeff2, __global numtyp4* coeff3,
                          __global numtyp4* coeff4, __global numtyp *cutsq,
                          const int lj_types, 
                          __global numtyp *sp_lj_in, __global int *dev_nbor, 
                          __global int *dev_packed, __global acctyp4 *ans,
                          __global acctyp *engv, const int eflag, 
                          const int vflag, const int inum,
                          const int nbor_pitch, const int t_per_atom, 
                          int tabstyle, int tablength) {
  bool shared_types = false;
  if (tabstyle == LOOKUP)
    lookup(x_, tabindex, coeff2, coeff3, coeff4, lj_types,
           cutsq, sp_lj_in, dev_nbor, dev_packed, 
           ans, engv, eflag, vflag, inum, 
           nbor_pitch, t_per_atom, tablength, shared_types);
  else if (tabstyle == LINEAR)
    linear(x_, tabindex, coeff2, coeff3, coeff4, lj_types,
           cutsq, sp_lj_in, dev_nbor, dev_packed, 
           ans, engv, eflag, vflag, inum, 
           nbor_pitch, t_per_atom, tablength, shared_types);
  else if (tabstyle == SPLINE)
    spline(x_, tabindex, coeff2, coeff3, coeff4, lj_types,
           cutsq, sp_lj_in, dev_nbor, dev_packed, 
           ans, engv, eflag, vflag, inum, 
           nbor_pitch, t_per_atom, tablength, shared_types);
  else // tabstyle == BITMAP
    bitmap(x_, tabindex, nshiftbits, nmask, 
           coeff2, coeff3, coeff4, lj_types,
           cutsq, sp_lj_in, dev_nbor, dev_packed, 
           ans, engv, eflag, vflag, inum, 
           nbor_pitch, t_per_atom, tablength, shared_types);
}


__kernel void kernel_pair_fast(__global numtyp4 *x_, __global int *tabindex,
                                __global int *nshiftbits, __global int *nmask,
                               __global numtyp4* coeff2, 
                               __global numtyp4 *coeff3,
                               __global numtyp4 *coeff4,
                               __global numtyp *cutsq_in,
                               __global numtyp* sp_lj_in, 
                               __global int *dev_nbor, __global int *dev_packed, 
                               __global acctyp4 *ans, __global acctyp *engv, 
                               const int eflag, const int vflag, const int inum, 
                               const int nbor_pitch, const int t_per_atom, 
                               int tabstyle, int tablength) {
  bool shared_types = true;
  int lj_types = 1; // dummy 
  if (tabstyle == LOOKUP)
    lookup(x_, tabindex, coeff2, coeff3, coeff4, lj_types,
           cutsq_in, sp_lj_in, dev_nbor, dev_packed, 
           ans, engv, eflag, vflag, inum, 
           nbor_pitch, t_per_atom, tablength, shared_types);
  else if (tabstyle == LINEAR)
    linear(x_, tabindex, coeff2, coeff3, coeff4, lj_types,
           cutsq_in, sp_lj_in, dev_nbor, dev_packed, 
           ans, engv, eflag, vflag, inum, 
           nbor_pitch, t_per_atom, tablength, shared_types);
  else if (tabstyle == SPLINE)
    spline(x_, tabindex, coeff2, coeff3, coeff4, lj_types,
           cutsq_in, sp_lj_in, dev_nbor, dev_packed, 
           ans, engv, eflag, vflag, inum, 
           nbor_pitch, t_per_atom, tablength, shared_types);
  else // tabstyle == BITMAP
    bitmap(x_, tabindex, nshiftbits, nmask, 
           coeff2, coeff3, coeff4, lj_types,
           cutsq_in, sp_lj_in, dev_nbor, dev_packed, 
           ans, engv, eflag, vflag, inum, 
           nbor_pitch, t_per_atom, tablength, shared_types);
  
}
