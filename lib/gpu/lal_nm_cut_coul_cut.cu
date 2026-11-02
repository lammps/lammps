// **************************************************************************
//                                 nm_cut_coul_cut.cu
//                             -------------------
//
//  Device code for acceleration of the nm/cut/coul/cut pair style
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                : October 2026
// ***************************************************************************

#if defined(NV_KERNEL) || defined(USE_HIP)

#include "lal_aux_fun1.h"
#ifndef _DOUBLE_DOUBLE
_texture( pos_tex,float4);
_texture( q_tex,float);
#else
_texture_2d( pos_tex,int4);
_texture( q_tex,int2);
#endif

#else
#define pos_tex x_
#define q_tex q_
#endif

__kernel void k_nm_cut_coul_cut(const __global numtyp4 *restrict x_,
                                const __global numtyp4 *restrict nm1,
                                const __global numtyp4 *restrict nm2,
                                const __global numtyp2 *restrict cutsq,
                                const int lj_types,
                                const __global numtyp *restrict sp_lj_in,
                                const __global int *dev_nbor,
                                const __global int *dev_packed,
                                __global acctyp3 *restrict ans,
                                __global acctyp *restrict engv,
                                const int eflag, const int vflag, const int inum,
                                const int nbor_pitch,
                                const __global numtyp *restrict q_,
                                const numtyp qqrd2e, const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  __local numtyp sp_lj[8];
  int n_stride;
  local_allocate_store_charge();

  sp_lj[0]=sp_lj_in[0];
  sp_lj[1]=sp_lj_in[1];
  sp_lj[2]=sp_lj_in[2];
  sp_lj[3]=sp_lj_in[3];
  sp_lj[4]=sp_lj_in[4];
  sp_lj[5]=sp_lj_in[5];
  sp_lj[6]=sp_lj_in[6];
  sp_lj[7]=sp_lj_in[7];

  acctyp3 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp energy, e_coul, virial[6];
  if (EVFLAG) {
    energy=(acctyp)0;
    e_coul=(acctyp)0;
    for (int i=0; i<6; i++) virial[i]=(acctyp)0;
  }

  if (ii<inum) {
    int nbor, nbor_end;
    int i, numj;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    numtyp qtmp; fetch(qtmp,i,q_tex);
    int itype=ix.w;

    for ( ; nbor<nbor_end; nbor+=n_stride) {
      ucl_prefetch(dev_packed+nbor+n_stride);

      int j=dev_packed[nbor];
      numtyp factor_lj, factor_coul;
      factor_lj = sp_lj[sbmask(j)];
      factor_coul = sp_lj[sbmask(j)+4];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int jtype=jx.w;

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      int mtype=itype*lj_types+jtype;
      if (rsq<cutsq[mtype].x) {
        numtyp r2inv = ucl_recip(rsq);
        numtyp forcecoul, force_nm, rninv, rminv;
        numtyp qj; fetch(qj,j,q_tex);
        if (rsq < cutsq[mtype].y)
          forcecoul = qqrd2e * qtmp*qj * ucl_rsqrt(rsq);
        else
          forcecoul = (numtyp)0.0;

        if (rsq < nm2[mtype].z) {
          rninv = ucl_powr(r2inv,(numtyp)0.5*nm2[mtype].x);
          rminv = ucl_powr(r2inv,(numtyp)0.5*nm2[mtype].y);
          force_nm = nm1[mtype].x*nm1[mtype].y*(nm1[mtype].z*rninv - nm1[mtype].w*rminv);
        } else
          force_nm = (numtyp)0.0;

        numtyp force = (factor_coul*forcecoul + factor_lj*force_nm) * r2inv;

        if (EVFLAG && eflag) {
          if (rsq < cutsq[mtype].y)
            e_coul += factor_coul*forcecoul;
          if (rsq < nm2[mtype].z) {
            numtyp e = nm1[mtype].x*(nm2[mtype].y*nm1[mtype].z*rninv -
                                     nm2[mtype].x*nm1[mtype].w*rminv);
            energy += factor_lj*(e-nm2[mtype].w);
          }
        }

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

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
  store_answers_q(f,energy,e_coul,virial,ii,inum,tid,t_per_atom,offset,eflag,
                  vflag,ans,engv);
}

__kernel void k_nm_cut_coul_cut_fast(const __global numtyp4 *restrict x_,
                                     const __global numtyp4 *restrict nm1_in,
                                     const __global numtyp4 *restrict nm2_in,
                                     const __global numtyp2 *restrict cutsq_in,
                                     const __global numtyp *restrict sp_lj_in,
                                     const __global int *dev_nbor,
                                     const __global int *dev_packed,
                                     __global acctyp3 *restrict ans,
                                     __global acctyp *restrict engv,
                                     const int eflag, const int vflag, const int inum,
                                     const int nbor_pitch,
                                     const __global numtyp *restrict q_,
                                     const numtyp qqrd2e, const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  __local numtyp4 nm1[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp4 nm2[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp2 cutsq[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp sp_lj[8];
  int n_stride;
  local_allocate_store_charge();

  if (tid<8)
    sp_lj[tid]=sp_lj_in[tid];
  if (tid<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
    nm1[tid]=nm1_in[tid];
    nm2[tid]=nm2_in[tid];
    cutsq[tid]=cutsq_in[tid];
  }

  acctyp3 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp energy, e_coul, virial[6];
  if (EVFLAG) {
    energy=(acctyp)0;
    e_coul=(acctyp)0;
    for (int i=0; i<6; i++) virial[i]=(acctyp)0;
  }

  __syncthreads();

  if (ii<inum) {
    int nbor, nbor_end;
    int i, numj;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    numtyp qtmp; fetch(qtmp,i,q_tex);
    int iw=ix.w;
    int itype=fast_mul((int)MAX_SHARED_TYPES,iw);

    for ( ; nbor<nbor_end; nbor+=n_stride) {
      ucl_prefetch(dev_packed+nbor+n_stride);

      int j=dev_packed[nbor];
      numtyp factor_lj, factor_coul;
      factor_lj = sp_lj[sbmask(j)];
      factor_coul = sp_lj[sbmask(j)+4];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int mtype=itype+jx.w;

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      if (rsq<cutsq[mtype].x) {
        numtyp r2inv = ucl_recip(rsq);
        numtyp forcecoul, force_nm, rninv, rminv;
        numtyp qj; fetch(qj,j,q_tex);
        if (rsq < cutsq[mtype].y)
          forcecoul = qqrd2e * qtmp*qj * ucl_rsqrt(rsq);
        else
          forcecoul = (numtyp)0.0;

        if (rsq < nm2[mtype].z) {
          rninv = ucl_powr(r2inv,(numtyp)0.5*nm2[mtype].x);
          rminv = ucl_powr(r2inv,(numtyp)0.5*nm2[mtype].y);
          force_nm = nm1[mtype].x*nm1[mtype].y*(nm1[mtype].z*rninv - nm1[mtype].w*rminv);
        } else
          force_nm = (numtyp)0.0;

        numtyp force = (factor_coul*forcecoul + factor_lj*force_nm) * r2inv;

        if (EVFLAG && eflag) {
          if (rsq < cutsq[mtype].y)
            e_coul += factor_coul*forcecoul;
          if (rsq < nm2[mtype].z) {
            numtyp e = nm1[mtype].x*(nm2[mtype].y*nm1[mtype].z*rninv -
                                     nm2[mtype].x*nm1[mtype].w*rminv);
            energy += factor_lj*(e-nm2[mtype].w);
          }
        }

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

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
  store_answers_q(f,energy,e_coul,virial,ii,inum,tid,t_per_atom,offset,eflag,
                  vflag,ans,engv);
}
