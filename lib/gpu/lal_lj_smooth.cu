// **************************************************************************
//                                   lj_smooth.cu
//                             -------------------
//                           Gurgen Melikyan (HSE University)
//
//  Device code for acceleration of the lj/smooth pair style
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                :
//    email                : gkmeliyan@edu.hse.ru
// ***************************************************************************

#if defined(NV_KERNEL) || defined(USE_HIP)
#include "lal_aux_fun1.h"
#ifndef _DOUBLE_DOUBLE
_texture( pos_tex,float4);
#else
_texture_2d( pos_tex,int4);
#endif
#else
#define pos_tex x_
#endif

__kernel void k_lj_smooth(const __global numtyp4 *restrict x_,
                   const __global numtyp4 *restrict lj1,
                   const __global numtyp4 *restrict lj3,
                   const __global numtyp4 *restrict ljsw,
                   const __global numtyp2 *restrict ljsw0,
                   const int lj_types,
                   const __global numtyp *restrict sp_lj,
                   const __global int * dev_nbor,
                   const __global int * dev_packed,
                   __global acctyp3 *restrict ans,
                   __global acctyp *restrict engv,
                   const int eflag, const int vflag, const int inum,
                   const int nbor_pitch, const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  int n_stride;
  local_allocate_store_pair();

  acctyp3 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp energy, virial[6];
  if (EVFLAG) {
    energy=(acctyp)0;
    for (int i=0; i<6; i++) virial[i]=(acctyp)0;
  }

  if (ii<inum) {
    int i, numj, nbor, nbor_end;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int itype=ix.w;

    numtyp force, r6inv, factor_lj, forcelj;
    numtyp r, t, tsq, fskin;

    for ( ; nbor<nbor_end; nbor+=n_stride) {
      ucl_prefetch(dev_packed+nbor+n_stride);

      int j=dev_packed[nbor];
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
        if (rsq < lj1[mtype].w) {
          r6inv = r2inv*r2inv*r2inv;
          forcelj = r6inv*(lj1[mtype].x*r6inv-lj1[mtype].y);
        }
        else {
          r = sqrt(rsq);
          t = r - ljsw0[mtype].y;
          tsq = t*t;
          fskin = ljsw[mtype].x + ljsw[mtype].y*t +
            ljsw[mtype].z*tsq + ljsw[mtype].w*tsq*t;
          forcelj = fskin*r;
        }
        force = factor_lj*r2inv*forcelj;

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (EVFLAG && eflag) {
          numtyp e;
          if (rsq < lj1[mtype].w)
            e = r6inv * (lj3[mtype].x*r6inv - lj3[mtype].y) - lj3[mtype].z;
          else
            e = ljsw0[mtype].x - ljsw[mtype].x*t -
              ljsw[mtype].y*tsq/2.0 - ljsw[mtype].z*tsq*t/3.0 -
              ljsw[mtype].w*tsq*tsq/4.0 - lj3[mtype].z;

          //numtyp e=r6inv*(lj3[mtype].x*r6inv-lj3[mtype].y);
          energy+=factor_lj*e;
        }
        if (EVFLAG && vflag) {
          virial[0] += delx*delx*force;
          virial[1] += dely*dely*force;
          virial[2] += delz*delz*force;
          virial[3] += delx*dely*force;
          virial[4] += delx*delz*force;
          virial[5] += dely*delz*force;
        }
      }

    } // for nbor
  } // if ii
  store_answers(f,energy,virial,ii,inum,tid,t_per_atom,offset,eflag,vflag,
                ans,engv);
}

__kernel void k_lj_smooth_fast(const __global numtyp4 *restrict x_,
                        const __global numtyp4 *restrict lj1_in,
                        const __global numtyp4 *restrict lj3_in,
                        const __global numtyp4 *restrict ljsw,
                        const __global numtyp2 *restrict ljsw0,
                        const __global numtyp *restrict sp_lj_in,
                        const __global int * dev_nbor,
                        const __global int * dev_packed,
                        __global acctyp3 *restrict ans,
                        __global acctyp *restrict engv,
                        const int eflag, const int vflag, const int inum,
                        const int nbor_pitch, const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  #ifndef ONETYPE
  __local numtyp4 lj1[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp4 lj3[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp sp_lj[4];
  if (tid<4)
    sp_lj[tid]=sp_lj_in[tid];
  if (tid<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
    lj1[tid]=lj1_in[tid];
    if (EVFLAG && eflag)
      lj3[tid]=lj3_in[tid];
  }
  __syncthreads();
  #else
  const numtyp lj1x=lj1_in[ONETYPE].x;
  const numtyp lj1y=lj1_in[ONETYPE].y;
  const numtyp cutsq=lj1_in[ONETYPE].z;
  numtyp lj3x, lj3y, lj3z;
  if (EVFLAG && eflag) {
    lj3x=lj3_in[ONETYPE].x;
    lj3y=lj3_in[ONETYPE].y;
    lj3z=lj3_in[ONETYPE].z;
  }
  #endif

  int n_stride;
  local_allocate_store_pair();

  acctyp3 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp energy, virial[6];
  if (EVFLAG) {
    energy=(acctyp)0;
    for (int i=0; i<6; i++) virial[i]=(acctyp)0;
  }

  if (ii<inum) {
    int i, numj, nbor, nbor_end;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    #ifndef ONETYPE
    int iw=ix.w;
    int itype=fast_mul((int)MAX_SHARED_TYPES,iw);

    numtyp force, r6inv, factor_lj, forcelj;
    numtyp r, t, tsq, fskin;

    #endif

    NOUNROLL
    for ( ; nbor<nbor_end; nbor+=n_stride) {
      ucl_prefetch(dev_packed+nbor+n_stride);
      int j=dev_packed[nbor];
      #ifndef ONETYPE
      factor_lj = sp_lj[sbmask(j)];
      j &= NEIGHMASK;
      #endif

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      #ifndef ONETYPE
      int mtype=itype+jx.w;
      #endif

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      if (rsq<lj1[mtype].z) {
        numtyp r2inv=ucl_recip(rsq);
        if (rsq < lj1[mtype].w) {
          r6inv = r2inv*r2inv*r2inv;
          forcelj = r6inv*(lj1[mtype].x*r6inv-lj1[mtype].y);
        }
        else {
          r = sqrt(rsq);
          t = r - ljsw0[mtype].y; //?
          tsq = t*t;
          fskin = ljsw[mtype].x + ljsw[mtype].y*t +
            ljsw[mtype].z*tsq + ljsw[mtype].w*tsq*t;
          forcelj = fskin*r;
        }
        force = factor_lj*r2inv*forcelj;

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (EVFLAG && eflag) {
          numtyp e;
          if (rsq < lj1[mtype].w)
            e = r6inv * (lj3[mtype].x*r6inv - lj3[mtype].y) - lj3[mtype].z;
          else
            e = ljsw0[mtype].x - ljsw[mtype].x*t -
              ljsw[mtype].y*tsq/2.0 - ljsw[mtype].z*tsq*t/3.0 -
              ljsw[mtype].w*tsq*tsq/4.0 - lj3[mtype].z; //???

          energy+=factor_lj*e;
        }
        if (EVFLAG && vflag) {
          virial[0] += delx*delx*force;
          virial[1] += dely*dely*force;
          virial[2] += delz*delz*force;
          virial[3] += delx*dely*force;
          virial[4] += delx*delz*force;
          virial[5] += dely*delz*force;
        }
      }

    } // for nbor
  } // if ii
  store_answers(f,energy,virial,ii,inum,tid,t_per_atom,offset,eflag,vflag,
                ans,engv);
}
