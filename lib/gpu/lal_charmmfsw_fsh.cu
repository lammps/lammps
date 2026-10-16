// **************************************************************************
//                                 charmmfsw_fsh.cu
//                             -------------------
//      Claude Opus 5.5 (Anthropic), under the direction of Axel Kohlmeyer
//
//  Device code for acceleration of the lj/charmmfsw/coul/charmmfsh pair style
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

// Lennard-Jones energy with CHARMM force switching, as in pair style
// lj/charmmfsw/coul/charmmfsh; fsw holds cut_lj6, denom_lj12, cut_lj6inv, cut_lj3,
// denom_lj6, cut_lj3inv, cut_lj_inner6inv, and cut_lj_inner3inv

ucl_inline numtyp charmmfsw_lj_energy(const numtyp rsq, const numtyp r6inv,
                                      const numtyp lj3, const numtyp lj4,
                                      const numtyp cut_lj_innersq,
                                      const __global numtyp *restrict fsw) {
  numtyp evdwl12, evdwl6;
  if (rsq > cut_lj_innersq) {
    const numtyp rinv = ucl_rsqrt(rsq);
    const numtyp r3inv = rinv*rinv*rinv;
    evdwl12 = lj3*fsw[0]*fsw[1]*(r6inv-fsw[2])*(r6inv-fsw[2]);
    evdwl6 = -lj4*fsw[3]*fsw[4]*(r3inv-fsw[5])*(r3inv-fsw[5]);
  } else {
    evdwl12 = r6inv*lj3*r6inv - lj3*fsw[6]*fsw[2];
    evdwl6 = -lj4*r6inv + lj4*fsw[7]*fsw[5];
  }
  return evdwl12 + evdwl6;
}

__kernel void k_charmmfsw_fsh(const __global numtyp4 *restrict x_,
                              const __global numtyp4 *restrict lj1,
                              const int lj_types,
                              const __global numtyp *restrict sp_lj,
                              const __global numtyp *restrict fsw,
                              const __global int *dev_nbor,
                              const __global int *dev_packed,
                              __global acctyp3 *restrict ans,
                              __global acctyp *restrict engv,
                              const int eflag, const int vflag, const int inum,
                              const int nbor_pitch,
                              const __global numtyp *restrict q_,
                              const numtyp cut_coulsq, const numtyp qqrd2e,
                              const numtyp cut_coulinv, const numtyp denom_lj,
                              const numtyp cut_bothsq, const numtyp cut_ljsq,
                              const numtyp cut_lj_innersq, const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  int n_stride;
  local_allocate_store_bio();

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
        numtyp r2inv=ucl_recip(rsq);
        numtyp forcecoul, force_lj, force, r6inv, prefactor, r, rinv, switch1;

        if (rsq < cut_ljsq) {
          r6inv = r2inv*r2inv*r2inv;
          force_lj = factor_lj*r6inv*(lj1[mtype].x*r6inv-lj1[mtype].y);
          if (rsq > cut_lj_innersq) {
            switch1 = (cut_ljsq-rsq);
            switch1 *= switch1;
            switch1 *= (cut_ljsq+(numtyp)2.0*rsq-(numtyp)3.0*cut_lj_innersq)/
                       denom_lj;
            force_lj *= switch1;
          }
        } else
          force_lj = (numtyp)0.0;

        if (rsq < cut_coulsq) {
          rinv = ucl_rsqrt(rsq);
          r = ucl_recip(rinv);
          fetch(prefactor,j,q_tex);
          prefactor *= factor_coul * qqrd2e * qtmp;
          forcecoul = prefactor * (rinv - r*cut_coulinv*cut_coulinv);
        } else
          forcecoul = (numtyp)0.0;

        force = (force_lj + forcecoul) * r2inv;

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (EVFLAG && eflag) {
          if (rsq < cut_coulsq)
            e_coul += prefactor*(rinv + cut_coulinv*cut_coulinv*r - (numtyp)2.0*cut_coulinv);
          if (rsq < cut_ljsq)
            energy+=factor_lj*charmmfsw_lj_energy(rsq,r6inv,lj1[mtype].z,lj1[mtype].w,
                                                  cut_lj_innersq,fsw);
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
  store_answers_q(f,energy,e_coul,virial,ii,inum,tid,t_per_atom,offset,eflag,
                  vflag,ans,engv);
}

__kernel void k_charmmfsw_fsh_fast(const __global numtyp4 *restrict x_,
                                   const __global numtyp2 *restrict ljd_in,
                                   const __global numtyp *restrict sp_lj_in,
                                   const __global numtyp *restrict fsw,
                                   const __global int *dev_nbor,
                                   const __global int *dev_packed,
                                   __global acctyp3 *restrict ans,
                                   __global acctyp *restrict engv,
                                   const int eflag, const int vflag,
                                   const int inum, const int nbor_pitch,
                                   const __global numtyp *restrict q_,
                                   const numtyp cut_coulsq, const numtyp qqrd2e,
                                   const numtyp cut_coulinv, const numtyp denom_lj,
                                   const numtyp cut_bothsq, const numtyp cut_ljsq,
                                   const numtyp cut_lj_innersq,
                                   const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  __local numtyp2 ljd[MAX_BIO_SHARED_TYPES];
  __local numtyp sp_lj[8];
  int n_stride;
  local_allocate_store_bio();

  if (tid<8)
    sp_lj[tid]=sp_lj_in[tid];
  if (tid<MAX_BIO_SHARED_TYPES)
    ljd[tid]=ljd_in[tid];
  if (tid+BLOCK_BIO_PAIR<MAX_BIO_SHARED_TYPES)
    ljd[tid+BLOCK_BIO_PAIR]=ljd_in[tid+BLOCK_BIO_PAIR];

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

      if (rsq<cut_bothsq) {
        numtyp r2inv=ucl_recip(rsq);
        numtyp forcecoul, force_lj, force, prefactor, r, rinv, switch1;
        numtyp lj3, lj4, r6inv;

        if (rsq < cut_ljsq) {
          // mixed coefficients lj3 = 4 eps sigma^12 and lj4 = 4 eps sigma^6
          numtyp eps = ucl_sqrt(ljd[itype].x*ljd[jtype].x);
          numtyp sig6 = (numtyp)0.5 * (ljd[itype].y+ljd[jtype].y);
          sig6 = sig6*sig6;
          sig6 = sig6*sig6*sig6;
          lj4 = (numtyp)4.0*eps*sig6;
          lj3 = lj4*sig6;
          r6inv = r2inv*r2inv*r2inv;
          force_lj = factor_lj*r6inv*((numtyp)12.0*lj3*r6inv-(numtyp)6.0*lj4);
          if (rsq > cut_lj_innersq) {
            switch1 = (cut_ljsq-rsq);
            switch1 *= switch1;
            switch1 *= (cut_ljsq+(numtyp)2.0*rsq-(numtyp)3.0*cut_lj_innersq)/
                       denom_lj;
            force_lj *= switch1;
          }
        } else
          force_lj = (numtyp)0.0;

        if (rsq < cut_coulsq) {
          rinv = ucl_rsqrt(rsq);
          r = ucl_recip(rinv);
          fetch(prefactor,j,q_tex);
          prefactor *= factor_coul * qqrd2e * qtmp;
          forcecoul = prefactor * (rinv - r*cut_coulinv*cut_coulinv);
        } else
          forcecoul = (numtyp)0.0;

        force = (force_lj + forcecoul) * r2inv;

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (EVFLAG && eflag) {
          if (rsq < cut_coulsq)
            e_coul += prefactor*(rinv + cut_coulinv*cut_coulinv*r - (numtyp)2.0*cut_coulinv);
          if (rsq < cut_ljsq)
            energy+=factor_lj*charmmfsw_lj_energy(rsq,r6inv,lj3,lj4,cut_lj_innersq,fsw);
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
  store_answers_q(f,energy,e_coul,virial,ii,inum,tid,t_per_atom,offset,eflag,
                  vflag,ans,engv);
}

