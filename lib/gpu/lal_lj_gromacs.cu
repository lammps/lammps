// **************************************************************************
//                               lj_gromacs.cu
//                             -------------------
//                           Trung Dac Nguyen (ORNL)
//
//  Device code for acceleration of the gromacs/coul/long pair style
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

__kernel void k_lj_gromacs(const __global numtyp4 *restrict x_,
                           const __global numtyp4 *restrict lj1,
                           const __global numtyp4 *restrict lj3,
                           const __global numtyp4 *restrict ljsw,
                           const int lj_types,
                           const __global numtyp *restrict sp_lj_in,
                           const __global int *dev_nbor,
                           const __global int *dev_packed,
                           __global acctyp4 *restrict ans,
                           __global acctyp *restrict engv, 
                           const int eflag, const int vflag, const int inum, 
                           const int nbor_pitch, const int t_per_atom) {
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

  if (ii<inum) {
    const __global int *nbor, *list_end;
    int i, numj;
    __local int n_stride;
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
      int jtype=jx.w;

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      int mtype=itype*lj_types+jtype;
      if (rsq<lj1[mtype].z) {
        numtyp r2inv=ucl_recip(rsq);
        numtyp force_lj, force, r6inv, t;
        
        r6inv = r2inv*r2inv*r2inv;
        force_lj = r6inv*(lj1[mtype].x*r6inv-lj1[mtype].y);
        if (rsq > lj1[mtype].w) {
          numtyp r = ucl_sqrt(rsq);
          t = r - lj3[mtype].z;
          numtyp fswitch = r*t*t*(ljsw[mtype].x + ljsw[mtype].y*t);
          force_lj += fswitch;
        } 

        force = factor_lj*force_lj * r2inv;

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (eflag>0) {
          numtyp e=r6inv*(lj3[mtype].x*r6inv-lj3[mtype].y);
          e += lj3[mtype].w;
          if (rsq > lj1[mtype].w) {
            numtyp eswitch = t*t*t*(ljsw[mtype].z + ljsw[mtype].w*t);
            e += eswitch;
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
    store_answers(f,energy,virial,ii,inum,tid,t_per_atom,offset,eflag,
                  vflag,ans,engv);
  } // if ii
}

__kernel void k_lj_gromacs_fast(const __global numtyp4 *restrict x_,
                                const __global numtyp4 *restrict lj1_in,
                                const __global numtyp4 *restrict lj3_in,
                                const __global numtyp4 *restrict ljsw_in,
                                const __global numtyp *restrict sp_lj_in,
                                const __global int *dev_nbor,
                                const __global int *dev_packed,
                                __global acctyp4 *restrict ans,
                                __global acctyp *restrict engv,
                                const int eflag, const int vflag, const int inum,
                                const int nbor_pitch, const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  __local numtyp4 lj1[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp4 lj3[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp4 ljsw[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp sp_lj[4];
  if (tid<4)
    sp_lj[tid]=sp_lj_in[tid];
  if (tid<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
    lj1[tid]=lj1_in[tid];
    lj3[tid]=lj3_in[tid];
    ljsw[tid]=ljsw_in[tid];
  }
  
  acctyp energy=(acctyp)0;
  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(acctyp)0;
  
  __syncthreads();
  
  if (ii<inum) {
    const __global int *nbor, *list_end;
    int i, numj;
    __local int n_stride;
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

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;
      
      if (rsq<lj1[mtype].z) {
        numtyp r2inv=ucl_recip(rsq);
        numtyp force_lj, force, r6inv, t;
        
        r6inv = r2inv*r2inv*r2inv;
        force_lj = r6inv*(lj1[mtype].x*r6inv-lj1[mtype].y);
        if (rsq > lj1[mtype].w) {
          numtyp r = ucl_sqrt(rsq);
          t = r - lj3[mtype].z;
          numtyp fswitch = r*t*t*(ljsw[mtype].x + ljsw[mtype].y*t);
          force_lj += fswitch;
        } 

        force = factor_lj*force_lj * r2inv;

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (eflag>0) {
          numtyp e=r6inv*(lj3[mtype].x*r6inv-lj3[mtype].y);
          e += lj3[mtype].w;
          if (rsq > lj1[mtype].w) {
            numtyp eswitch = t*t*t*(ljsw[mtype].z + ljsw[mtype].w*t);
            e += eswitch;
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
    store_answers(f,energy,virial,ii,inum,tid,t_per_atom,offset,eflag,
                  vflag,ans,engv);
  } // if ii
}

