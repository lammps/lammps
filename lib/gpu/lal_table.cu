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

__kernel void kernel_pair(__global numtyp4 *x_, __global numtyp4 *coeff1,
                          __global numtyp4* coeff2, __global numtyp4* coeff3,
                          __global numtyp4* coeff4, __global numtyp *cutsq,
                          const int lj_types, 
                          __global numtyp *sp_lj_in, __global int *dev_nbor, 
                          __global int *dev_packed, __global acctyp4 *ans,
                          __global acctyp *engv, const int eflag, 
                          const int vflag, const int inum,
                          const int nbor_pitch, const int t_per_atom, 
                          int tabstyle, int tablength) {
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
  
  int tlm1;
  if (tabstyle != BITMAP) tlm1 = tablength - 1;
  else tlm1 = (1 << tablength) - 1;
  
  if (ii<inum) {
    __global int *nbor, *list_end;
    int i, numj, n_stride;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,list_end,nbor);
  
    numtyp4 ix=fetch_pos(i,x_); //x_[i];
    int itype=ix.w;

    numtyp factor_lj;
    for ( ; nbor<list_end; nbor+=n_stride) {
  
      int j=*nbor;
      factor_lj = sp_lj[sbmask(j)];
      j &= NEIGHMASK;

      numtyp4 jx=fetch_pos(j,x_); //x_[j];
      int jtype=jx.w;

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;
        
      int mtype=itype*lj_types+jtype;
      int tabindex = coeff1[mtype].x;
      if (rsq<cutsq[mtype]) {
        int itable;
        numtyp fraction=(numtyp)0;
        numtyp a = (numtyp)0;
        numtyp b = (numtyp)0;
        numtyp value = (numtyp)0;
        numtyp force = (numtyp)0;
        union_int_float rsq_lookup;
        if (tabstyle == LOOKUP) { 
          itable = (rsq - coeff2[mtype].x) * coeff2[mtype].y;
          if (itable < tlm1) {
            itable += tabindex*tablength;
            force = factor_lj * coeff3[itable].z;
          } else force = (numtyp)0.0;
        } else if (tabstyle == LINEAR) {  
          itable = (rsq - coeff2[mtype].x) * coeff2[mtype].y;
          if (itable < tlm1) {
            itable += tabindex*tablength;
            fraction = (rsq - coeff3[itable].x) * coeff2[mtype].y;
            value = coeff3[itable].z + fraction*coeff4[itable].z;
            force = factor_lj * value;
          } else force = (numtyp)0.0;
        } else if (tabstyle == SPLINE) {
          itable = (rsq - coeff2[mtype].x) * coeff2[mtype].y;
          if (itable < tlm1) {
            itable += tabindex*tablength;
            b = (rsq - coeff3[itable].x) * coeff2[mtype].y;
            a = (numtyp)1.0 - b;
            value = a * coeff3[itable].z + b * coeff3[itable+1].z + 
              ((a*a*a-a)*coeff4[itable].z + (b*b*b-b)*coeff4[itable+1].z) * 
                    coeff2[mtype].z;
            force = factor_lj * value;
          } else force = (numtyp)0.0;
        } else { 
          rsq_lookup.f = rsq;
   	      itable = rsq_lookup.i & ((int)coeff1[mtype].w);
	        itable >>= (int)coeff1[mtype].z;
	        fraction = (rsq_lookup.f - coeff3[itable].x) * coeff4[itable].w;
	        value = coeff3[itable].z + fraction*coeff4[itable].z;
	        force = factor_lj * value;
        }
  
        force*=factor_lj;
      
        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (eflag>0) {
          numtyp e = (numtyp)0.0;
          if (itable < tlm1) {
            if (tabstyle == LOOKUP) {
              e = coeff3[itable].y;
            } else if (tabstyle == LINEAR || tabstyle == BITMAP) {
              e = coeff3[itable].y + fraction*coeff4[itable].y;
            }
            else {
              e = a * coeff3[itable].y + b * coeff3[itable+1].y + 
                ((a*a*a-a)*coeff4[itable].y + (b*b*b-b)*coeff4[itable+1].y) * 
                  coeff2[mtype].z;
            }
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

__kernel void kernel_pair_fast(__global numtyp4 *x_, __global numtyp4 *coeff1,
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
  
  int tlm1;
  if (tabstyle != BITMAP) tlm1 = tablength - 1;
  else tlm1 = (1 << tablength) - 1;
  
  if (ii<inum) {
    __global int *nbor, *list_end;
    int i, numj, n_stride;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,list_end,nbor);

    numtyp4 ix=fetch_pos(i,x_); //x_[i];
    int iw=ix.w;
    int itype=fast_mul((int)MAX_SHARED_TYPES,iw);

    numtyp factor_lj;
    for ( ; nbor<list_end; nbor+=n_stride) {
  
      int j=*nbor;
      factor_lj = sp_lj[sbmask(j)];
      j &= NEIGHMASK;

      numtyp4 jx=fetch_pos(j,x_); //x_[j];
      int mtype=itype+jx.w;
      int tabindex = coeff1[mtype].x;
      
      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;
          
      if (rsq<cutsq[mtype]) {
        int itable=0;
        numtyp fraction=(numtyp)0;
        numtyp a = (numtyp)0;
        numtyp b = (numtyp)0;
        numtyp value = (numtyp)0;
        numtyp force = (numtyp)0;
        union_int_float rsq_lookup; 
        if (tabstyle == LOOKUP) {
          itable = (rsq - coeff2[mtype].x) * coeff2[mtype].y;
          if (itable < tlm1) {
            itable += tabindex*tablength;
            force = factor_lj * coeff3[itable].z;
          } else force = (numtyp)0.0;
        } else if (tabstyle == LINEAR) {
          itable = (rsq - coeff2[mtype].x) * coeff2[mtype].y;
          if (itable < tlm1) {
            itable += tabindex*tablength;
            fraction = (rsq - coeff3[itable].x) * coeff2[mtype].y;
            value = coeff3[itable].z + fraction*coeff4[itable].z;
            force = factor_lj * value;
          } else force = (numtyp)0.0;
        } else if (tabstyle == SPLINE) {
          itable = (rsq - coeff2[mtype].x) * coeff2[mtype].y;
          if (itable < tlm1) {
            itable += tabindex*tablength;
            b = (rsq - coeff3[itable].x) * coeff2[mtype].y;
            a = (numtyp)1.0 - b;
            value = a * coeff3[itable].z + b * coeff3[itable+1].z + 
              ((a*a*a-a)*coeff4[itable].z + (b*b*b-b)*coeff4[itable+1].z) * 
                    coeff2[mtype].z;
            force = factor_lj * value;
          } else force = (numtyp)0.0;
        } else { 
          rsq_lookup.f = rsq;
   	      itable = rsq_lookup.i & ((int)coeff1[mtype].w);
	        itable >>= (int)coeff1[mtype].z;
	        fraction = (rsq_lookup.f - coeff3[itable].x) * coeff4[itable].w;
	        value = coeff3[itable].z + fraction*coeff4[itable].z;
	        force = factor_lj * value;
        }
              
        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (eflag>0) {
          numtyp e = (numtyp)0.0;
          if (itable < tlm1) {
            if (tabstyle == LOOKUP) {
              e = coeff3[itable].y;
            } else if (tabstyle == LINEAR || tabstyle == BITMAP) {
              e = coeff3[itable].y + fraction*coeff4[itable].y;
            }
            else {
              e = a * coeff3[itable].y + b * coeff3[itable+1].y + 
                ((a*a*a-a)*coeff4[itable].y + (b*b*b-b)*coeff4[itable+1].y) * 
                  coeff2[mtype].z;
            }
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

