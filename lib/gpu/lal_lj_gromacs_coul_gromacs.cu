// **************************************************************************
//                                 lj_gromacs_coul_gromacs.cu
//                             -------------------
//
//  Device code for acceleration of the lj/gromacs/coul/gromacs pair style
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

__kernel void k_lj_gromacs_coul_gromacs(const __global numtyp4 *restrict x_,
                                        const __global numtyp4 *restrict lj1,
                                        const __global numtyp4 *restrict lj3,
                                        const __global numtyp *restrict ljsw5,
                                        const int lj_types,
                                        const __global numtyp *restrict sp_lj_in,
                                        const __global int *dev_nbor,
                                        const __global int *dev_packed,
                                        __global acctyp3 *restrict ans,
                                        __global acctyp *restrict engv,
                                        const int eflag, const int vflag, const int inum,
                                        const int nbor_pitch,
                                        const __global numtyp *restrict q_,
                                        const numtyp qqrd2e, const numtyp cut_bothsq,
                                        const numtyp cut_ljsq, const numtyp cut_lj_innersq,
                                        const numtyp cut_lj_inner, const numtyp cut_coulsq,
                                        const numtyp cut_coul_innersq,
                                        const numtyp cut_coul_inner, const numtyp coulsw1,
                                        const numtyp coulsw2, const numtyp coulsw3,
                                        const numtyp coulsw4, const numtyp coulsw5,
                                        const int t_per_atom) {
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
      if (rsq<cut_bothsq) {
        numtyp r2inv = ucl_recip(rsq);
        numtyp qj; fetch(qj,j,q_tex);
        numtyp forcecoul, force_lj, r6inv, tlj;

        // skip if qi or qj = 0.0 since this potential may be used as
        // coarse-grain model with many uncharged atoms

        if (rsq < cut_coulsq && qtmp != (numtyp)0.0 && qj != (numtyp)0.0) {
          forcecoul = qqrd2e * qtmp*qj * ucl_rsqrt(rsq);
          if (rsq > cut_coul_innersq) {
            numtyp r = ucl_sqrt(rsq);
            numtyp tc = r - cut_coul_inner;
            forcecoul += qqrd2e * qtmp*qj * r*tc*tc * (coulsw1 + coulsw2*tc);
          }
        } else
          forcecoul = (numtyp)0.0;

        if (rsq < cut_ljsq) {
          r6inv = r2inv*r2inv*r2inv;
          force_lj = r6inv * (lj1[mtype].x*r6inv - lj1[mtype].y);
          if (rsq > cut_lj_innersq) {
            numtyp r = ucl_sqrt(rsq);
            tlj = r - cut_lj_inner;
            force_lj += r*tlj*tlj * (lj1[mtype].z + lj1[mtype].w*tlj);
          }
        } else
          force_lj = (numtyp)0.0;

        numtyp force = (factor_coul*forcecoul + factor_lj*force_lj) * r2inv;

        if (EVFLAG && eflag) {
          if (rsq < cut_coulsq) {
            numtyp e = qqrd2e * qtmp*qj * (ucl_rsqrt(rsq) - coulsw5);
            if (rsq > cut_coul_innersq) {
              numtyp tc = ucl_sqrt(rsq) - cut_coul_inner;
              e += qqrd2e * qtmp*qj * tc*tc*tc * (coulsw3 + coulsw4*tc);
            }
            e_coul += factor_coul*e;
          }
          if (rsq < cut_ljsq) {
            numtyp e = r6inv*(lj3[mtype].x*r6inv - lj3[mtype].y) + ljsw5[mtype];
            if (rsq > cut_lj_innersq)
              e += tlj*tlj*tlj * (lj3[mtype].z + lj3[mtype].w*tlj);
            energy += factor_lj*e;
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

__kernel void k_lj_gromacs_coul_gromacs_fast(const __global numtyp4 *restrict x_,
                                             const __global numtyp4 *restrict lj1_in,
                                             const __global numtyp4 *restrict lj3_in,
                                             const __global numtyp *restrict ljsw5_in,
                                             const __global numtyp *restrict sp_lj_in,
                                             const __global int *dev_nbor,
                                             const __global int *dev_packed,
                                             __global acctyp3 *restrict ans,
                                             __global acctyp *restrict engv,
                                             const int eflag, const int vflag,
                                             const int inum, const int nbor_pitch,
                                             const __global numtyp *restrict q_,
                                             const numtyp qqrd2e, const numtyp cut_bothsq,
                                             const numtyp cut_ljsq,
                                             const numtyp cut_lj_innersq,
                                             const numtyp cut_lj_inner,
                                             const numtyp cut_coulsq,
                                             const numtyp cut_coul_innersq,
                                             const numtyp cut_coul_inner,
                                             const numtyp coulsw1, const numtyp coulsw2,
                                             const numtyp coulsw3, const numtyp coulsw4,
                                             const numtyp coulsw5, const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  __local numtyp4 lj1[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp4 lj3[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp ljsw5[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp sp_lj[8];
  int n_stride;
  local_allocate_store_charge();

  if (tid<8)
    sp_lj[tid]=sp_lj_in[tid];
  if (tid<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
    lj1[tid]=lj1_in[tid];
    if (EVFLAG && eflag) lj3[tid]=lj3_in[tid];
    if (EVFLAG && eflag) ljsw5[tid]=ljsw5_in[tid];
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

      if (rsq<cut_bothsq) {
        numtyp r2inv = ucl_recip(rsq);
        numtyp qj; fetch(qj,j,q_tex);
        numtyp forcecoul, force_lj, r6inv, tlj;

        // skip if qi or qj = 0.0 since this potential may be used as
        // coarse-grain model with many uncharged atoms

        if (rsq < cut_coulsq && qtmp != (numtyp)0.0 && qj != (numtyp)0.0) {
          forcecoul = qqrd2e * qtmp*qj * ucl_rsqrt(rsq);
          if (rsq > cut_coul_innersq) {
            numtyp r = ucl_sqrt(rsq);
            numtyp tc = r - cut_coul_inner;
            forcecoul += qqrd2e * qtmp*qj * r*tc*tc * (coulsw1 + coulsw2*tc);
          }
        } else
          forcecoul = (numtyp)0.0;

        if (rsq < cut_ljsq) {
          r6inv = r2inv*r2inv*r2inv;
          force_lj = r6inv * (lj1[mtype].x*r6inv - lj1[mtype].y);
          if (rsq > cut_lj_innersq) {
            numtyp r = ucl_sqrt(rsq);
            tlj = r - cut_lj_inner;
            force_lj += r*tlj*tlj * (lj1[mtype].z + lj1[mtype].w*tlj);
          }
        } else
          force_lj = (numtyp)0.0;

        numtyp force = (factor_coul*forcecoul + factor_lj*force_lj) * r2inv;

        if (EVFLAG && eflag) {
          if (rsq < cut_coulsq) {
            numtyp e = qqrd2e * qtmp*qj * (ucl_rsqrt(rsq) - coulsw5);
            if (rsq > cut_coul_innersq) {
              numtyp tc = ucl_sqrt(rsq) - cut_coul_inner;
              e += qqrd2e * qtmp*qj * tc*tc*tc * (coulsw3 + coulsw4*tc);
            }
            e_coul += factor_coul*e;
          }
          if (rsq < cut_ljsq) {
            numtyp e = r6inv*(lj3[mtype].x*r6inv - lj3[mtype].y) + ljsw5[mtype];
            if (rsq > cut_lj_innersq)
              e += tlj*tlj*tlj * (lj3[mtype].z + lj3[mtype].w*tlj);
            energy += factor_lj*e;
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
