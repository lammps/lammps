// **************************************************************************
//                                dipole_lj_sf.cu
//                             -------------------
//                           Trung Dac Nguyen (ORNL)
//
//  Device code for acceleration of the dipole/sf pair style
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                :
//    email                : nguyentd@ornl.gov
// ***************************************************************************

#if defined(NV_KERNEL) || defined(USE_HIP)
#include "lal_aux_fun1.h"

#ifndef _DOUBLE_DOUBLE
_texture( pos_tex,float4);
_texture( q_tex,float);
_texture( mu_tex,float4);
#else
_texture_2d( pos_tex,int4);
_texture( q_tex,int2);
_texture_2d( mu_tex,int4);
#endif

#else
#define pos_tex x_
#define q_tex q_
#define mu_tex mu_
#endif

#if (SHUFFLE_AVAIL == 0)

#define store_answers_tq(f, tor, energy, e_coul, virial, ii, inum, tid,     \
                         t_per_atom, offset, eflag, vflag, ans, engv)       \
  if (t_per_atom>1) {                                                       \
    simd_reduce_add6(t_per_atom, red_acc, offset, tid, f.x, f.y, f.z,       \
                     tor.x, tor.y, tor.z);                                  \
    if (EVFLAG && (vflag==2 || eflag==2)) {                                 \
      if (eflag) {                                                          \
        simdsync();                                                         \
        simd_reduce_add2(t_per_atom, red_acc, offset, tid, energy, e_coul); \
      }                                                                     \
      if (vflag) {                                                          \
        simdsync();                                                         \
        simd_reduce_arr(6, t_per_atom, red_acc, offset, tid, virial);       \
      }                                                                     \
    }                                                                       \
  }                                                                         \
  if (offset==0 && ii<inum) {                                               \
    ans[ii]=f;                                                              \
    ans[ii+inum]=tor;                                                       \
  }                                                                         \
  if (EVFLAG && (eflag || vflag)) {                                         \
    int ei=BLOCK_ID_X;                                                      \
    if (eflag!=2 && vflag!=2) {                                             \
      const int ev_stride=NUM_BLOCKS_X;                                     \
      if (eflag) {                                                          \
        simdsync();                                                         \
        block_reduce_add2(simd_size(), red_acc, tid, energy, e_coul);       \
        if (vflag) __syncthreads();                                         \
        if (tid==0) {                                                       \
          engv[ei]=energy*(acctyp)0.5;                                      \
          ei+=ev_stride;                                                    \
          engv[ei]=e_coul*(acctyp)0.5;                                      \
          ei+=ev_stride;                                                    \
        }                                                                   \
      }                                                                     \
      if (vflag) {                                                          \
        simdsync();                                                         \
        block_reduce_arr(6, simd_size(), red_acc, tid, virial);             \
        if (tid==0) {                                                       \
          for (int r=0; r<6; r++) {                                         \
            engv[ei]=virial[r]*(acctyp)0.5;                                 \
            ei+=ev_stride;                                                  \
          }                                                                 \
        }                                                                   \
      }                                                                     \
    } else if (offset==0 && ii<inum) {                                      \
      int ei=ii;                                                            \
      if (EVFLAG && eflag) {                                                \
        engv[ei]=energy*(acctyp)0.5;                                        \
        ei+=inum;                                                           \
        engv[ei]=e_coul*(acctyp)0.5;                                        \
        ei+=inum;                                                           \
      }                                                                     \
      if (EVFLAG && vflag) {                                                \
        for (int i=0; i<6; i++) {                                           \
          engv[ei]=virial[i]*(acctyp)0.5;                                   \
          ei+=inum;                                                         \
        }                                                                   \
      }                                                                     \
    }                                                                       \
  }

#else

#if (EVFLAG == 1)

#define store_answers_tq(f, tor, energy, e_coul, virial, ii, inum, tid,     \
                         t_per_atom, offset, eflag, vflag, ans, engv)       \
  if (t_per_atom>1) {                                                       \
    simd_reduce_add6(t_per_atom, f.x, f.y, f.z, tor.x, tor.y, tor.z);       \
    if (vflag==2 || eflag==2) {                                             \
      if (eflag)                                                            \
        simd_reduce_add2(t_per_atom,energy,e_coul);                         \
      if (vflag)                                                            \
        simd_reduce_arr(6, t_per_atom,virial);                              \
    }                                                                       \
  }                                                                         \
  if (offset==0 && ii<inum) {                                               \
    ans[ii]=f;                                                              \
    ans[ii+inum]=tor;                                                       \
  }                                                                         \
  if (eflag || vflag) {                                                     \
    if (eflag!=2 && vflag!=2) {                                             \
      const int vwidth = simd_size();                                       \
      const int voffset = tid & (simd_size() - 1);                          \
      const int bnum = tid/simd_size();                                     \
      int active_subgs = BLOCK_SIZE_X/simd_size();                          \
      for ( ; active_subgs > 1; active_subgs /= vwidth) {                   \
        if (active_subgs < BLOCK_SIZE_X/simd_size()) __syncthreads();       \
        if (bnum < active_subgs) {                                          \
          if (eflag) {                                                      \
            simd_reduce_add2(vwidth, energy, e_coul);                       \
            if (voffset==0) {                                               \
              red_acc[6][bnum] = energy;                                    \
              red_acc[7][bnum] = e_coul;                                    \
            }                                                               \
          }                                                                 \
          if (vflag) {                                                      \
            simd_reduce_arr(6, vwidth, virial);                             \
            if (voffset==0)                                                 \
              for (int r=0; r<6; r++) red_acc[r][bnum]=virial[r];           \
          }                                                                 \
        }                                                                   \
                                                                            \
        __syncthreads();                                                    \
        if (tid < active_subgs) {                                           \
          if (eflag) {                                                      \
            energy = red_acc[6][tid];                                       \
            e_coul = red_acc[7][tid];                                       \
          }                                                                 \
          if (vflag)                                                        \
            for (int r = 0; r < 6; r++) virial[r] = red_acc[r][tid];        \
        } else {                                                            \
          if (eflag) energy = e_coul = (acctyp)0;                           \
          if (vflag) for (int r = 0; r < 6; r++) virial[r] = (acctyp)0;     \
        }                                                                   \
      }                                                                     \
                                                                            \
      if (bnum == 0) {                                                      \
        int ei=BLOCK_ID_X;                                                  \
        const int ev_stride=NUM_BLOCKS_X;                                   \
        if (eflag) {                                                        \
          simd_reduce_add2(vwidth, energy, e_coul);                         \
          if (tid==0) {                                                     \
            engv[ei]=energy*(acctyp)0.5;                                    \
            ei+=ev_stride;                                                  \
            engv[ei]=e_coul*(acctyp)0.5;                                    \
            ei+=ev_stride;                                                  \
          }                                                                 \
        }                                                                   \
        if (vflag) {                                                        \
          simd_reduce_arr(6, vwidth, virial);                               \
          if (tid==0) {                                                     \
            for (int r=0; r<6; r++) {                                       \
              engv[ei]=virial[r]*(acctyp)0.5;                               \
              ei+=ev_stride;                                                \
            }                                                               \
          }                                                                 \
        }                                                                   \
      }                                                                     \
    } else if (offset==0 && ii<inum) {                                      \
      int ei=ii;                                                            \
      if (eflag) {                                                          \
        engv[ei]=energy*(acctyp)0.5;                                        \
        ei+=inum;                                                           \
        engv[ei]=e_coul*(acctyp)0.5;                                        \
        ei+=inum;                                                           \
      }                                                                     \
      if (vflag) {                                                          \
        for (int i=0; i<6; i++) {                                           \
          engv[ei]=virial[i]*(acctyp)0.5;                                   \
          ei+=inum;                                                         \
        }                                                                   \
      }                                                                     \
    }                                                                       \
  }

#else

#define store_answers_tq(f, tor, energy, e_coul, virial, ii, inum, tid,     \
                        t_per_atom, offset, eflag, vflag, ans, engv)        \
  if (t_per_atom>1)                                                         \
    simd_reduce_add6(t_per_atom, f.x, f.y, f.z, tor.x, tor.y, tor.z);       \
  if (offset==0 && ii<inum) {                                               \
    ans[ii]=f;                                                              \
    ans[ii+inum]=tor;                                                       \
  }

#endif
#endif

__kernel void k_dipole_lj_sf(const __global numtyp4 *restrict x_,
                             const __global numtyp4 *restrict lj1,
                             const __global numtyp4 *restrict lj3,
                             const int lj_types,
                             const __global numtyp *restrict sp_lj_in,
                             const __global int *dev_nbor,
                             const __global int *dev_packed,
                             __global acctyp4 *restrict ans,
                             __global acctyp *restrict engv,
                             const int eflag, const int vflag, const int inum,
                             const int nbor_pitch,
                             const __global numtyp *restrict q_ ,
                             const __global numtyp4 *restrict mu_,
                             const __global numtyp *restrict cutsq,
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

  acctyp4 f, tor;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  tor.x=(acctyp)0; tor.y=(acctyp)0; tor.z=(acctyp)0;
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
    numtyp4 mui; fetch4(mui,i,mu_tex); //mu_[i];
    int itype=ix.w;

    for ( ; nbor<nbor_end; nbor+=n_stride) {
      int j=dev_packed[nbor];

      numtyp factor_lj, factor_coul;
      factor_lj = sp_lj[sbmask(j)];
      factor_coul = sp_lj[sbmask(j)+4];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      numtyp qj; fetch(qj,j,q_tex);
      numtyp4 muj; fetch4(muj,j,mu_tex); //mu_[j];
      int jtype=jx.w;

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      int mtype=itype*lj_types+jtype;
      if (rsq<cutsq[mtype]) {
        numtyp r2inv=ucl_recip(rsq);
        numtyp force_lj, r6inv;
        numtyp rinv, r3inv, r5inv;
        numtyp pre1, pre2, pre4;
        numtyp pdotp, pidotr, pjdotr;
        numtyp presf,afac,bfac,pqfac,qpfac,rcutlj2inv,rcutlj6inv,rcutcoul2inv;
        numtyp4 aforcecoul, bforcecoul;

        acctyp4 forcecoul, ticoul;
        acctyp4 force;

        forcecoul.x = forcecoul.y = forcecoul.z = (acctyp)0;
        ticoul.x = ticoul.y = ticoul.z = (acctyp)0;

        if (rsq < lj1[mtype].z) {
          r6inv = r2inv*r2inv*r2inv;
          numtyp forceljcut = factor_lj*r6inv*(lj1[mtype].x*r6inv-lj1[mtype].y)*r2inv;

          rcutlj2inv = ucl_recip(lj1[mtype].z);
          rcutlj6inv = rcutlj2inv * rcutlj2inv * rcutlj2inv;
          numtyp forceljsf = rcutlj6inv*(lj1[mtype].x*rcutlj6inv-lj1[mtype].y)*rcutlj2inv;

          force_lj = factor_lj * (forceljcut - forceljsf);
        } else force_lj = (numtyp)0.0;

        if (rsq < lj1[mtype].w) {
          rinv = ucl_rsqrt(rsq);
          rcutcoul2inv = ucl_recip(lj1[mtype].w);

          // charge-charge
          if (qtmp != (numtyp)0.0 && qj != (numtyp)0.0) {
            r3inv = r2inv*rinv;
            pre1 = qtmp*qj*rinv*(r2inv-rcutcoul2inv);

            forcecoul.x += pre1*delx;
            forcecoul.y += pre1*dely;
            forcecoul.z += pre1*delz;
          }

          // dipole-dipole
          if (mui.w > (numtyp)0.0 && muj.w > (numtyp)0.0) {
            r3inv = r2inv*rinv;
            r5inv = r3inv*r2inv;

            pdotp  = mui.x*muj.x + mui.y*muj.y + mui.z*muj.z;
            pidotr = mui.x*delx + mui.y*dely + mui.z*delz;
            pjdotr = muj.x*delx + muj.y*dely + muj.z*delz;

            afac = (numtyp)1.0 - rsq*rsq * rcutcoul2inv*rcutcoul2inv;
            pre1 = afac * (pdotp - (numtyp)3.0*r2inv*pidotr*pjdotr);
            aforcecoul.x = pre1*delx;
            aforcecoul.y = pre1*dely;
            aforcecoul.z = pre1*delz;

            bfac = (numtyp)1.0-(numtyp)4.0*rsq*ucl_sqrt(rsq)*rcutcoul2inv*ucl_sqrt(rcutcoul2inv)+
              (numtyp)3.0*rsq*rsq*rcutcoul2inv*rcutcoul2inv;
            presf = (numtyp)2.0*r2inv*pidotr*pjdotr;
            bforcecoul.x = bfac * (pjdotr*mui.x+pidotr*muj.x-presf*delx);
            bforcecoul.y = bfac * (pjdotr*mui.y+pidotr*muj.y-presf*dely);
            bforcecoul.z = bfac * (pjdotr*mui.z+pidotr*muj.z-presf*delz);

            forcecoul.x += (numtyp)3.0*r5inv*(aforcecoul.x + bforcecoul.x);
            forcecoul.y += (numtyp)3.0*r5inv*(aforcecoul.y + bforcecoul.y);
            forcecoul.z += (numtyp)3.0*r5inv*(aforcecoul.z + bforcecoul.z);

            pre2 = (numtyp)3.0*bfac*r5inv*pjdotr;
            pre4 = -bfac*r3inv;

            numtyp crossx = pre4 * (mui.y*muj.z - mui.z*muj.y);
            numtyp crossy = pre4 * (mui.z*muj.x - mui.x*muj.z);
            numtyp crossz = pre4 * (mui.x*muj.y - mui.y*muj.x);

            ticoul.x += crossx + pre2 * (mui.y*delz - mui.z*dely);
            ticoul.y += crossy + pre2 * (mui.z*delx - mui.x*delz);
            ticoul.z += crossz + pre2 * (mui.x*dely - mui.y*delx);
          }

          // dipole-charge
          if (mui.w > (numtyp)0.0 && qj != (numtyp)0.0) {
            r3inv = r2inv*rinv;
            r5inv = r3inv*r2inv;
            pidotr = mui.x*delx + mui.y*dely + mui.z*delz;
            rcutcoul2inv=ucl_recip(lj1[mtype].w);
            pre1 = (numtyp)3.0*qj*r5inv * pidotr*((numtyp)1.0-rsq*rcutcoul2inv);
            pqfac = (numtyp)1.0 - (numtyp)3.0*rsq*rcutcoul2inv +
              (numtyp)2.0*rsq*ucl_sqrt(rsq)*rcutcoul2inv*ucl_sqrt(rcutcoul2inv);
            pre2 = qj*r3inv * pqfac;

            forcecoul.x += pre2*mui.x - pre1*delx;
            forcecoul.y += pre2*mui.y - pre1*dely;
            forcecoul.z += pre2*mui.z - pre1*delz;
            ticoul.x += pre2 * (mui.y*delz - mui.z*dely);
            ticoul.y += pre2 * (mui.z*delx - mui.x*delz);
            ticoul.z += pre2 * (mui.x*dely - mui.y*delx);
          }

          // charge-dipole
          if (muj.w > (numtyp)0.0 && qtmp != (numtyp)0.0) {
            r3inv = r2inv*rinv;
            r5inv = r3inv*r2inv;
            pjdotr = muj.x*delx + muj.y*dely + muj.z*delz;
            rcutcoul2inv=ucl_recip(lj1[mtype].w);
            pre1 = (numtyp)3.0*qtmp*r5inv * pjdotr*((numtyp)1.0-rsq*rcutcoul2inv);
            qpfac = (numtyp)1.0 - (numtyp)3.0*rsq*rcutcoul2inv +
              (numtyp)2.0*rsq*ucl_sqrt(rsq)*rcutcoul2inv*ucl_sqrt(rcutcoul2inv);
            pre2 = qtmp*r3inv * qpfac;

            forcecoul.x += pre1*delx - pre2*muj.x;
            forcecoul.y += pre1*dely - pre2*muj.y;
            forcecoul.z += pre1*delz - pre2*muj.z;
          }
        } else {
          forcecoul.x = forcecoul.y = forcecoul.z = (acctyp)0;
          ticoul.x = ticoul.y = ticoul.z = (acctyp)0;
        }

        numtyp fq = factor_coul*qqrd2e;
        force.x = fq*forcecoul.x + delx*force_lj;
        force.y = fq*forcecoul.y + dely*force_lj;
        force.z = fq*forcecoul.z + delz*force_lj;
        f.x+=force.x;
        f.y+=force.y;
        f.z+=force.z;
        tor.x+=fq*ticoul.x;
        tor.y+=fq*ticoul.y;
        tor.z+=fq*ticoul.z;

        if (EVFLAG && eflag) {
          acctyp e = (acctyp)0.0;
          if (rsq < lj1[mtype].w) {
            numtyp fac = (numtyp)1.0-ucl_sqrt(rsq*rcutcoul2inv);
            e = qtmp*qj*rinv*fac*fac;
            if (mui.w > (numtyp)0.0 && muj.w > (numtyp)0.0)
              e += bfac* (r3inv*pdotp - (numtyp)3.0*r5inv*pidotr*pjdotr);
            if (mui.w > (numtyp)0.0 && qj != (numtyp)0.0)
              e += -qj*r3inv*pidotr * pqfac;
            if (muj.w > (numtyp)0.0 && qtmp != (numtyp)0.0)
              e += qtmp*r3inv*pjdotr * qpfac;
              e *= fq;
          } else e = (acctyp)0.0;
          e_coul += e;

          if (rsq < lj1[mtype].z) {
            e=r6inv*(lj3[mtype].x*r6inv-lj3[mtype].y) +
              rcutlj6inv*((numtyp)6.0*lj3[mtype].x*rcutlj6inv -
              (numtyp)3.0*lj3[mtype].y)*rsq*rcutlj2inv +
              rcutlj6inv*((numtyp)(-7.0)*lj3[mtype].x*rcutlj6inv +
              (numtyp)4.0*lj3[mtype].y);
            energy+=factor_lj*e;
          }
        }
        if (EVFLAG && vflag) {
          virial[0] += delx*force.x;
          virial[1] += dely*force.y;
          virial[2] += delz*force.z;
          virial[3] += delx*force.y;
          virial[4] += delx*force.z;
          virial[5] += dely*force.z;
        }
      }
    } // for nbor
  } // if ii
  store_answers_tq(f,tor,energy,e_coul,virial,ii,inum,tid,t_per_atom,offset,
                   eflag,vflag,ans,engv);
}

__kernel void k_dipole_lj_sf_fast(const __global numtyp4 *restrict x_,
                                  const __global numtyp4 *restrict lj1_in,
                                  const __global numtyp4 *restrict lj3_in,
                                  const __global numtyp *restrict sp_lj_in,
                                  const __global int *dev_nbor,
                                  const __global int *dev_packed,
                                  __global acctyp4 *restrict ans,
                                  __global acctyp *restrict engv,
                                  const int eflag, const int vflag,
                                  const int inum, const int nbor_pitch,
                                  const __global numtyp *restrict q_,
                                  const __global numtyp4 *restrict mu_,
                                  const __global numtyp *restrict _cutsq,
                                  const numtyp qqrd2e,
                                  const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  __local numtyp4 lj1[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp4 lj3[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp cutsq[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp sp_lj[8];
  int n_stride;
  local_allocate_store_charge();

  if (tid<8)
    sp_lj[tid]=sp_lj_in[tid];
  if (tid<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
    lj1[tid]=lj1_in[tid];
    cutsq[tid]=_cutsq[tid];
    if (EVFLAG && eflag)
      lj3[tid]=lj3_in[tid];
  }

  acctyp4 f, tor;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  tor.x=(acctyp)0; tor.y=(acctyp)0; tor.z=(acctyp)0;
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
    numtyp4 mui; fetch4(mui,i,mu_tex); //mu_[i];
    int iw=ix.w;
    int itype=fast_mul((int)MAX_SHARED_TYPES,iw);

    for ( ; nbor<nbor_end; nbor+=n_stride) {
      int j=dev_packed[nbor];

      numtyp factor_lj, factor_coul;
      factor_lj = sp_lj[sbmask(j)];
      factor_coul = sp_lj[sbmask(j)+4];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      numtyp qj; fetch(qj,j,q_tex);
      numtyp4 muj; fetch4(muj,j,mu_tex); //mu_[j];
      int mtype=itype+jx.w;

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      if (rsq<cutsq[mtype]) {
        numtyp r2inv=ucl_recip(rsq);
        numtyp force_lj, r6inv;
        numtyp rinv, r3inv, r5inv;
        numtyp pre1, pre2, pre4;
        numtyp pdotp, pidotr, pjdotr;
        numtyp presf,afac,bfac,pqfac,qpfac,rcutlj2inv,rcutlj6inv,rcutcoul2inv;
        numtyp4 aforcecoul, bforcecoul;

        acctyp4 forcecoul, ticoul;
        acctyp4 force;

        forcecoul.x = forcecoul.y = forcecoul.z = (acctyp)0;
        ticoul.x = ticoul.y = ticoul.z = (acctyp)0;

        if (rsq < lj1[mtype].z) {
          r6inv = r2inv*r2inv*r2inv;
          numtyp forceljcut = factor_lj*r6inv*(lj1[mtype].x*r6inv-lj1[mtype].y)*r2inv;

          rcutlj2inv = ucl_recip(lj1[mtype].z);
          rcutlj6inv = rcutlj2inv * rcutlj2inv * rcutlj2inv;
          numtyp forceljsf = rcutlj6inv*(lj1[mtype].x*rcutlj6inv-lj1[mtype].y)*rcutlj2inv;

          force_lj = factor_lj * (forceljcut - forceljsf);
        } else force_lj = (numtyp)0.0;

        if (rsq < lj1[mtype].w) {
          rinv = ucl_rsqrt(rsq);
          rcutcoul2inv = ucl_recip(lj1[mtype].w);

          // charge-charge
          if (qtmp != (numtyp)0.0 && qj != (numtyp)0.0) {
            r3inv = r2inv*rinv;
            pre1 = qtmp*qj*rinv*(r2inv-rcutcoul2inv);

            forcecoul.x += pre1*delx;
            forcecoul.y += pre1*dely;
            forcecoul.z += pre1*delz;
          }

          // dipole-dipole
          if (mui.w > (numtyp)0.0 && muj.w > (numtyp)0.0) {
            r3inv = r2inv*rinv;
            r5inv = r3inv*r2inv;

            pdotp  = mui.x*muj.x + mui.y*muj.y + mui.z*muj.z;
            pidotr = mui.x*delx + mui.y*dely + mui.z*delz;
            pjdotr = muj.x*delx + muj.y*dely + muj.z*delz;

            afac = (numtyp)1.0 - rsq*rsq * rcutcoul2inv*rcutcoul2inv;
            pre1 = afac * (pdotp - (numtyp)3.0*r2inv*pidotr*pjdotr);
            aforcecoul.x = pre1*delx;
            aforcecoul.y = pre1*dely;
            aforcecoul.z = pre1*delz;

            bfac = (numtyp)1.0-(numtyp)4.0*rsq*ucl_sqrt(rsq)*rcutcoul2inv*ucl_sqrt(rcutcoul2inv)+
              (numtyp)3.0*rsq*rsq*rcutcoul2inv*rcutcoul2inv;
            presf = (numtyp)2.0*r2inv*pidotr*pjdotr;
            bforcecoul.x = bfac * (pjdotr*mui.x+pidotr*muj.x-presf*delx);
            bforcecoul.y = bfac * (pjdotr*mui.y+pidotr*muj.y-presf*dely);
            bforcecoul.z = bfac * (pjdotr*mui.z+pidotr*muj.z-presf*delz);

            forcecoul.x += (numtyp)3.0*r5inv*(aforcecoul.x + bforcecoul.x);
            forcecoul.y += (numtyp)3.0*r5inv*(aforcecoul.y + bforcecoul.y);
            forcecoul.z += (numtyp)3.0*r5inv*(aforcecoul.z + bforcecoul.z);

            pre2 = (numtyp)3.0*bfac*r5inv*pjdotr;
            pre4 = -bfac*r3inv;

            numtyp crossx = pre4 * (mui.y*muj.z - mui.z*muj.y);
            numtyp crossy = pre4 * (mui.z*muj.x - mui.x*muj.z);
            numtyp crossz = pre4 * (mui.x*muj.y - mui.y*muj.x);

            ticoul.x += crossx + pre2 * (mui.y*delz - mui.z*dely);
            ticoul.y += crossy + pre2 * (mui.z*delx - mui.x*delz);
            ticoul.z += crossz + pre2 * (mui.x*dely - mui.y*delx);
          }

          // dipole-charge
          if (mui.w > (numtyp)0.0 && qj != (numtyp)0.0) {
            r3inv = r2inv*rinv;
            r5inv = r3inv*r2inv;
            pidotr = mui.x*delx + mui.y*dely + mui.z*delz;
            pre1 = (numtyp)3.0*qj*r5inv * pidotr*((numtyp)1.0-rsq*rcutcoul2inv);
            pqfac = (numtyp)1.0 - (numtyp)3.0*rsq*rcutcoul2inv +
              (numtyp)2.0*rsq*ucl_sqrt(rsq)*rcutcoul2inv*ucl_sqrt(rcutcoul2inv);
            pre2 = qj*r3inv * pqfac;

            forcecoul.x += pre2*mui.x - pre1*delx;
            forcecoul.y += pre2*mui.y - pre1*dely;
            forcecoul.z += pre2*mui.z - pre1*delz;
            ticoul.x += pre2 * (mui.y*delz - mui.z*dely);
            ticoul.y += pre2 * (mui.z*delx - mui.x*delz);
            ticoul.z += pre2 * (mui.x*dely - mui.y*delx);
          }

          // charge-dipole
          if (muj.w > (numtyp)0.0 && qtmp != (numtyp)0.0) {
            r3inv = r2inv*rinv;
            r5inv = r3inv*r2inv;
            pjdotr = muj.x*delx + muj.y*dely + muj.z*delz;

            pre1 = (numtyp)3.0*qtmp*r5inv * pjdotr*((numtyp)1.0-rsq*rcutcoul2inv);
            qpfac = (numtyp)1.0 - (numtyp)3.0*rsq*rcutcoul2inv +
              (numtyp)2.0*rsq*ucl_sqrt(rsq)*rcutcoul2inv*ucl_sqrt(rcutcoul2inv);
            pre2 = qtmp*r3inv * qpfac;

            forcecoul.x += pre1*delx - pre2*muj.x;
            forcecoul.y += pre1*dely - pre2*muj.y;
            forcecoul.z += pre1*delz - pre2*muj.z;
          }
        } else {
          forcecoul.x = forcecoul.y = forcecoul.z = (acctyp)0;
          ticoul.x = ticoul.y = ticoul.z = (acctyp)0;
        }

        numtyp fq = factor_coul*qqrd2e;
        force.x = fq*forcecoul.x + delx*force_lj;
        force.y = fq*forcecoul.y + dely*force_lj;
        force.z = fq*forcecoul.z + delz*force_lj;
        f.x+=force.x;
        f.y+=force.y;
        f.z+=force.z;
        tor.x+=fq*ticoul.x;
        tor.y+=fq*ticoul.y;
        tor.z+=fq*ticoul.z;

        if (EVFLAG && eflag) {
          acctyp e = (acctyp)0.0;
          if (rsq < lj1[mtype].w) {
            numtyp fac = (numtyp)1.0-ucl_sqrt(rsq*rcutcoul2inv);
            e = qtmp*qj*rinv*fac*fac;
            if (mui.w > (numtyp)0.0 && muj.w > (numtyp)0.0)
              e += bfac* (r3inv*pdotp - (numtyp)3.0*r5inv*pidotr*pjdotr);
            if (mui.w > (numtyp)0.0 && qj != (numtyp)0.0)
              e += -qj*r3inv*pidotr * pqfac;
            if (muj.w > (numtyp)0.0 && qtmp != (numtyp)0.0)
              e += qtmp*r3inv*pjdotr * qpfac;
            e *= fq;
          } else e = (acctyp)0.0;
          e_coul += e;

          if (rsq < lj1[mtype].z) {
            e=r6inv*(lj3[mtype].x*r6inv-lj3[mtype].y) +
              rcutlj6inv*((numtyp)6.0*lj3[mtype].x*rcutlj6inv -
              (numtyp)3.0*lj3[mtype].y)*rsq*rcutlj2inv +
              rcutlj6inv*((numtyp)(-7.0)*lj3[mtype].x*rcutlj6inv +
              (numtyp)4.0*lj3[mtype].y);
            energy+=factor_lj*e;
          }
        }
        if (EVFLAG && vflag) {
          virial[0] += delx*force.x;
          virial[1] += dely*force.y;
          virial[2] += delz*force.z;
          virial[3] += delx*force.y;
          virial[4] += delx*force.z;
          virial[5] += dely*force.z;
        }
      }

    } // for nbor
  } // if ii
  store_answers_tq(f,tor,energy,e_coul,virial,ii,inum,tid,t_per_atom,offset,
                   eflag,vflag,ans,engv);
}

