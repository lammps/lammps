// **************************************************************************
//                                dipole_lj.cu
//                             -------------------
//                           Trung Dac Nguyen (ORNL)
//
//  Device code for acceleration of the dipole/cut pair style
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                :
//    email                : nguyentd@ornl.gov
// ***************************************************************************/

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

#if (ARCH < 300)

#define store_answers_tq(f, tor, energy, ecoul, virial, ii, inum, tid,      \
                        t_per_atom, offset, eflag, vflag, ans, engv)        \
  if (t_per_atom>1) {                                                       \
    __local acctyp red_acc[8][BLOCK_PAIR];                                  \
    red_acc[0][tid]=f.x;                                                    \
    red_acc[1][tid]=f.y;                                                    \
    red_acc[2][tid]=f.z;                                                    \
    red_acc[3][tid]=tor.x;                                                  \
    red_acc[4][tid]=tor.y;                                                  \
    red_acc[5][tid]=tor.z;                                                  \
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                         \
      if (offset < s) {                                                     \
        for (int r=0; r<6; r++)                                             \
          red_acc[r][tid] += red_acc[r][tid+s];                             \
      }                                                                     \
    }                                                                       \
    f.x=red_acc[0][tid];                                                    \
    f.y=red_acc[1][tid];                                                    \
    f.z=red_acc[2][tid];                                                    \
    tor.x=red_acc[3][tid];                                                  \
    tor.y=red_acc[4][tid];                                                  \
    tor.z=red_acc[5][tid];                                                  \
    if (eflag>0 || vflag>0) {                                               \
      for (int r=0; r<6; r++)                                               \
        red_acc[r][tid]=virial[r];                                          \
      red_acc[6][tid]=energy;                                               \
      red_acc[7][tid]=ecoul;                                                \
      for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                       \
        if (offset < s) {                                                   \
          for (int r=0; r<8; r++)                                           \
            red_acc[r][tid] += red_acc[r][tid+s];                           \
        }                                                                   \
      }                                                                     \
      for (int r=0; r<6; r++)                                               \
        virial[r]=red_acc[r][tid];                                          \
      energy=red_acc[6][tid];                                               \
      ecoul=red_acc[7][tid];                                                \
    }                                                                       \
  }                                                                         \
  if (offset==0) {                                                          \
    int ei=ii;                                                              \
    if (eflag>0) {                                                          \
      engv[ei]=energy*(acctyp)0.5;                                             \
      ei+=inum;                                                           \
      engv[ei]=e_coul*(acctyp)0.5;                                             \
      ei+=inum;                                                           \
    }                                                                       \
    if (vflag>0) {                                                          \
      for (int i=0; i<6; i++) {                                             \
        engv[ei]=virial[i]*(acctyp)0.5;                                        \
        ei+=inum;                                                         \
      }                                                                     \
    }                                                                       \
    ans[ii]=f;                                                              \
    ans[ii+inum]=tor;                                                       \
  }

#else

#define store_answers_tq(f, tor, energy, e_coul, virial, ii, inum, tid,     \
                        t_per_atom, offset, eflag, vflag, ans, engv)        \
  if (t_per_atom>1) {                                                       \
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                         \
        f.x += shfl_xor(f.x, s, t_per_atom);                                \
        f.y += shfl_xor(f.y, s, t_per_atom);                                \
        f.z += shfl_xor(f.z, s, t_per_atom);                                \
        tor.x += shfl_xor(tor.x, s, t_per_atom);                            \
        tor.y += shfl_xor(tor.y, s, t_per_atom);                            \
        tor.z += shfl_xor(tor.z, s, t_per_atom);                            \
        energy += shfl_xor(energy, s, t_per_atom);                          \
        e_coul += shfl_xor(e_coul, s, t_per_atom);                          \
    }                                                                       \
    if (vflag>0) {                                                          \
      for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                       \
          for (int r=0; r<6; r++)                                           \
            virial[r] += shfl_xor(virial[r], s, t_per_atom);                \
      }                                                                     \
    }                                                                       \
  }                                                                         \
  if (offset==0) {                                                          \
    int ei=ii;                                                              \
    if (eflag>0) {                                                          \
      engv[ei]=energy*(acctyp)0.5;                                             \
      ei+=inum;                                                           \
      engv[ei]=e_coul*(acctyp)0.5;                                             \
      ei+=inum;                                                           \
    }                                                                       \
    if (vflag>0) {                                                          \
      for (int i=0; i<6; i++) {                                             \
        engv[ei]=virial[i]*(acctyp)0.5;                                        \
        ei+=inum;                                                         \
      }                                                                     \
    }                                                                       \
    ans[ii]=f;                                                              \
    ans[ii+inum]=tor;                                                       \
  }

#endif

#define MY_PIS (acctyp)1.77245385090551602729

__kernel void k_dipole_long_lj(const __global numtyp4 *restrict x_,
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
                          const __global numtyp *restrict q_,
                          const __global numtyp4 *restrict mu_,
                          const __global numtyp *restrict cutsq,
                          const numtyp cut_coulsq, const numtyp qqrd2e,
                          const numtyp g_ewald, const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  __local numtyp sp_lj[8];
  sp_lj[0]=sp_lj_in[0];
  sp_lj[1]=sp_lj_in[1];
  sp_lj[2]=sp_lj_in[2];
  sp_lj[3]=sp_lj_in[3];
  sp_lj[4]=sp_lj_in[4];
  sp_lj[5]=sp_lj_in[5];
  sp_lj[6]=sp_lj_in[6];
  sp_lj[7]=sp_lj_in[7];

  acctyp energy=(acctyp)0;
  acctyp e_coul=(acctyp)0;
  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp4 tor;
  tor.x=(acctyp)0;
  tor.y=(acctyp)0;
  tor.z=(acctyp)0;
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(acctyp)0;

  numtyp pre1 = numtyp(2.0) * g_ewald / MY_PIS;
  numtyp pre2 = numtyp(4.0) * (g_ewald*g_ewald*g_ewald) / MY_PIS;
  numtyp pre3 = numtyp(8.0) * (g_ewald*g_ewald*g_ewald*g_ewald*g_ewald) / MY_PIS;

  if (ii<inum) {
    int nbor, nbor_end;
    int i, numj;
    __local int n_stride;
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
        numtyp force_lj,rinv,r6inv;
        numtyp pdotp, pidotr, pjdotr, _erfc;
        numtyp g0,g1,g2,b0,b1,b2,b3,d0,d1,d2,d3;
        numtyp zdix,zdiy,zdiz,zdjx,zdjy,zdjz,zaix,zaiy,zaiz,zajx,zajy,zajz;
        numtyp g0b1_g1b2_g2b3,g0d1_g1d2_g2d3,facm1;
        numtyp fdx,fdy,fdz,fax,fay,faz;
        acctyp4 forcecoul, ticoul;
        acctyp4 force;

        forcecoul.x = forcecoul.y = forcecoul.z = (acctyp)0;
        ticoul.x = ticoul.y = ticoul.z = (acctyp)0;

        if (rsq < lj1[mtype].z) {
          r6inv = r2inv*r2inv*r2inv;
          force_lj = factor_lj*r6inv*(lj1[mtype].x*r6inv-lj1[mtype].y)*r2inv;
        } else force_lj = (numtyp)0.0;

        if (rsq < cut_coulsq) {
          rinv = ucl_rsqrt(rsq);
          numtyp r = ucl_rsqrt(r2inv);
          numtyp grij = g_ewald * r;
          numtyp expm2 = ucl_exp(-grij*grij);
          numtyp t = ucl_recip((numtyp)1.0 + EWALD_P*grij);
          _erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;

          pdotp  = mui.x*muj.x + mui.y*muj.y + mui.z*muj.z;
          pidotr = mui.x*delx + mui.y*dely + mui.z*delz;
          pjdotr = muj.x*delx + muj.y*dely + muj.z*delz;

          g0 = qtmp*qj;
          g1 = qtmp*pjdotr - qj*pidotr + pdotp;
          g2 = -pidotr*pjdotr;

          if (factor_coul > (numtyp)0.0) {
            b0 = _erfc * rinv;
            b1 = (b0 + pre1*expm2) * r2inv;
            b2 = ((numtyp)3.0*b1 + pre2*expm2) * r2inv;
            b3 = ((numtyp)5.0*b2 + pre3*expm2) * r2inv;

            g0b1_g1b2_g2b3 = g0*b1 + g1*b2 + g2*b3;
            fdx = delx * g0b1_g1b2_g2b3 -
              b1 * (qtmp*muj.x - qj*mui.x) +
              b2 * (pjdotr*mui.x + pidotr*muj.x);
            fdy = dely * g0b1_g1b2_g2b3 -
              b1 * (qtmp*muj.y - qj*mui.y) +
              b2 * (pjdotr*mui.y + pidotr*muj.y);
            fdz = delz * g0b1_g1b2_g2b3 -
              b1 * (qtmp*muj.z - qj*mui.z) +
              b2 * (pjdotr*mui.z + pidotr*muj.z);

            zdix = delx * (qj*b1 + b2*pjdotr) - b1*muj.x;
            zdiy = dely * (qj*b1 + b2*pjdotr) - b1*muj.y;
            zdiz = delz * (qj*b1 + b2*pjdotr) - b1*muj.z;
            zdjx = delx * (-qtmp*b1 + b2*pidotr) - b1*mui.x;
            zdjy = dely * (-qtmp*b1 + b2*pidotr) - b1*mui.y;
            zdjz = delz * (-qtmp*b1 + b2*pidotr) - b1*mui.z;

            if (factor_coul < (numtyp)1.0) {
              fdx *= factor_coul;
              fdy *= factor_coul;
              fdz *= factor_coul;
              zdix *= factor_coul;
              zdiy *= factor_coul;
              zdiz *= factor_coul;
              zdjx *= factor_coul;
              zdjy *= factor_coul;
              zdjz *= factor_coul;
            }
          } else {
            fdx = fdy = fdz = (numtyp)0.0;
            zdix = zdiy = zdiz = (numtyp)0.0;
            zdjx = zdjy = zdjz = (numtyp)0.0;
          }

          if (factor_coul < (numtyp)1.0) {
            d0 = (_erfc - (numtyp)1.0) * rinv;
            d1 = (d0 + pre1*expm2) * r2inv;
            d2 = ((numtyp)3.0*d1 + pre2*expm2) * r2inv;
            d3 = ((numtyp)5.0*d2 + pre3*expm2) * r2inv;

            g0d1_g1d2_g2d3 = g0*d1 + g1*d2 + g2*d3;
            fax = delx * g0d1_g1d2_g2d3 -
              d1 * (qtmp*muj.x - qj*mui.x) +
              d2 * (pjdotr*mui.x + pidotr*muj.x);
            fay = dely * g0d1_g1d2_g2d3 -
              d1 * (qtmp*muj.y - qj*mui.y) +
              d2 * (pjdotr*mui.y + pidotr*muj.y);
            faz = delz * g0d1_g1d2_g2d3 -
              d1 * (qtmp*muj.z - qj*mui.z) +
              d2 * (pjdotr*mui.z + pidotr*muj.z);

            zaix = delx * (qj*d1 + d2*pjdotr) - d1*muj.x;
            zaiy = dely * (qj*d1 + d2*pjdotr) - d1*muj.y;
            zaiz = delz * (qj*d1 + d2*pjdotr) - d1*muj.z;
            zajx = delx * (-qtmp*d1 + d2*pidotr) - d1*mui.x;
            zajy = dely * (-qtmp*d1 + d2*pidotr) - d1*mui.y;
            zajz = delz * (-qtmp*d1 + d2*pidotr) - d1*mui.z;

            if (factor_coul > (numtyp)0.0) {
              facm1 = (numtyp)1.0 - factor_coul;
              fax *= facm1;
              fay *= facm1;
              faz *= facm1;
              zaix *= facm1;
              zaiy *= facm1;
              zaiz *= facm1;
              zajx *= facm1;
              zajy *= facm1;
              zajz *= facm1;
            }
          } else {
            fax = fay = faz = (numtyp)0.0;
            zaix = zaiy = zaiz = (numtyp)0.0;
            zajx = zajy = zajz = (numtyp)0.0;
          }

          forcecoul.x = fdx + fax;
          forcecoul.y = fdy + fay;
          forcecoul.z = fdz + faz;

          ticoul.x = mui.y*(zdiz + zaiz) - mui.z*(zdiy + zaiy);
          ticoul.y = mui.z*(zdix + zaix) - mui.x*(zdiz + zaiz);
          ticoul.z = mui.x*(zdiy + zaiy) - mui.y*(zdix + zaix);

        } else {
          forcecoul.x = forcecoul.y = forcecoul.z = (numtyp)0.0;
          ticoul.x = ticoul.y = ticoul.z = (numtyp)0.0;
        }

        force.x = qqrd2e*forcecoul.x + delx*force_lj;
        force.y = qqrd2e*forcecoul.y + dely*force_lj;
        force.z = qqrd2e*forcecoul.z + delz*force_lj;
        f.x+=force.x;
        f.y+=force.y;
        f.z+=force.z;
        tor.x+=qqrd2e*ticoul.x;
        tor.y+=qqrd2e*ticoul.y;
        tor.z+=qqrd2e*ticoul.z;

        if (eflag>0) {
          acctyp e = (acctyp)0.0;
          if (rsq < cut_coulsq && factor_coul > (numtyp)0.0) {
            e = qqrd2e*(b0*g0 + b1*g1 + b2*g2);
            if (factor_coul < (numtyp)1.0) {
              e_coul *= factor_coul;
              e_coul += ((numtyp)1.0-factor_coul) * qqrd2e * (d0*g0 + d1*g1 + d2*g2);
            }
          } else e = (acctyp)0.0;
          e_coul += e;

          if (rsq < lj1[mtype].z) {
            e=r6inv*(lj3[mtype].x*r6inv-lj3[mtype].y);
            energy+=factor_lj*(e-lj3[mtype].z);
          }
        }
        if (vflag>0) {
          virial[0] += delx*force.x;
          virial[1] += dely*force.y;
          virial[2] += delz*force.z;
          virial[3] += delx*force.y;
          virial[4] += delx*force.z;
          virial[5] += dely*force.z;
        }
      }

    } // for nbor
    store_answers_tq(f,tor,energy,e_coul,virial,ii,inum,tid,t_per_atom,offset,eflag,
                    vflag,ans,engv);
  } // if ii
}

__kernel void k_dipole_long_lj_fast(const __global numtyp4 *restrict x_,
                               const __global numtyp4 *restrict lj1_in,
                               const __global numtyp4 *restrict lj3_in,
                               const __global numtyp *restrict sp_lj_in,
                               const __global int *dev_nbor,
                               const __global int *dev_packed,
                               __global acctyp4 *restrict ans,
                               __global acctyp *restrict engv,
                               const int eflag, const int vflag, const int inum,
                               const int nbor_pitch,
                               const __global numtyp *restrict q_,
                               const __global numtyp4 *restrict mu_,
                               const __global numtyp *restrict _cutsq,
                               const numtyp cut_coulsq, const numtyp qqrd2e,
                               const numtyp g_ewald, const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  __local numtyp4 lj1[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp4 lj3[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp cutsq[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp sp_lj[8];
  if (tid<8)
    sp_lj[tid]=sp_lj_in[tid];
  if (tid<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
    lj1[tid]=lj1_in[tid];
    cutsq[tid]=_cutsq[tid];
    if (eflag>0)
      lj3[tid]=lj3_in[tid];
  }

  acctyp energy=(acctyp)0;
  acctyp e_coul=(acctyp)0;
  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp4 tor;
  tor.x=(acctyp)0;
  tor.y=(acctyp)0;
  tor.z=(acctyp)0;
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(acctyp)0;

  __syncthreads();

  numtyp pre1 = numtyp(2.0) * g_ewald / MY_PIS;
  numtyp pre2 = numtyp(4.0) * (g_ewald*g_ewald*g_ewald) / MY_PIS;
  numtyp pre3 = numtyp(8.0) * (g_ewald*g_ewald*g_ewald*g_ewald*g_ewald) / MY_PIS;

  if (ii<inum) {
    int nbor, nbor_end;
    int i, numj;
    __local int n_stride;
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
        numtyp force_lj,rinv,r6inv;
        numtyp pdotp, pidotr, pjdotr, _erfc;
        numtyp g0,g1,g2,b0,b1,b2,b3,d0,d1,d2,d3;
        numtyp zdix,zdiy,zdiz,zdjx,zdjy,zdjz,zaix,zaiy,zaiz,zajx,zajy,zajz;
        numtyp g0b1_g1b2_g2b3,g0d1_g1d2_g2d3,facm1;
        numtyp fdx,fdy,fdz,fax,fay,faz;
        acctyp4 forcecoul, ticoul;
        acctyp4 force;

        forcecoul.x = forcecoul.y = forcecoul.z = (acctyp)0;
        ticoul.x = ticoul.y = ticoul.z = (acctyp)0;

        if (rsq < lj1[mtype].z) {
          r6inv = r2inv*r2inv*r2inv;
          force_lj = factor_lj*r6inv*(lj1[mtype].x*r6inv-lj1[mtype].y)*r2inv;
        } else force_lj = (numtyp)0.0;

        if (rsq < cut_coulsq) {
          rinv = ucl_rsqrt(rsq);
          numtyp r = ucl_rsqrt(r2inv);
          numtyp grij = g_ewald * r;
          numtyp expm2 = ucl_exp(-grij*grij);
          numtyp t = ucl_recip((numtyp)1.0 + EWALD_P*grij);
          _erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;

          pdotp  = mui.x*muj.x + mui.y*muj.y + mui.z*muj.z;
          pidotr = mui.x*delx + mui.y*dely + mui.z*delz;
          pjdotr = muj.x*delx + muj.y*dely + muj.z*delz;

          g0 = qtmp*qj;
          g1 = qtmp*pjdotr - qj*pidotr + pdotp;
          g2 = -pidotr*pjdotr;

          if (factor_coul > (numtyp)0.0) {
            b0 = _erfc * rinv;
            b1 = (b0 + pre1*expm2) * r2inv;
            b2 = ((numtyp)3.0*b1 + pre2*expm2) * r2inv;
            b3 = ((numtyp)5.0*b2 + pre3*expm2) * r2inv;

            g0b1_g1b2_g2b3 = g0*b1 + g1*b2 + g2*b3;
            fdx = delx * g0b1_g1b2_g2b3 -
              b1 * (qtmp*muj.x - qj*mui.x) +
              b2 * (pjdotr*mui.x + pidotr*muj.x);
            fdy = dely * g0b1_g1b2_g2b3 -
              b1 * (qtmp*muj.y - qj*mui.y) +
              b2 * (pjdotr*mui.y + pidotr*muj.y);
            fdz = delz * g0b1_g1b2_g2b3 -
              b1 * (qtmp*muj.z - qj*mui.z) +
              b2 * (pjdotr*mui.z + pidotr*muj.z);

            zdix = delx * (qj*b1 + b2*pjdotr) - b1*muj.x;
            zdiy = dely * (qj*b1 + b2*pjdotr) - b1*muj.y;
            zdiz = delz * (qj*b1 + b2*pjdotr) - b1*muj.z;
            zdjx = delx * (-qtmp*b1 + b2*pidotr) - b1*mui.x;
            zdjy = dely * (-qtmp*b1 + b2*pidotr) - b1*mui.y;
            zdjz = delz * (-qtmp*b1 + b2*pidotr) - b1*mui.z;

            if (factor_coul < (numtyp)1.0) {
              fdx *= factor_coul;
              fdy *= factor_coul;
              fdz *= factor_coul;
              zdix *= factor_coul;
              zdiy *= factor_coul;
              zdiz *= factor_coul;
              zdjx *= factor_coul;
              zdjy *= factor_coul;
              zdjz *= factor_coul;
            }
          } else {
            fdx = fdy = fdz = (numtyp)0.0;
            zdix = zdiy = zdiz = (numtyp)0.0;
            zdjx = zdjy = zdjz = (numtyp)0.0;
          }

          if (factor_coul < (numtyp)1.0) {
            d0 = (_erfc - (numtyp)1.0) * rinv;
            d1 = (d0 + pre1*expm2) * r2inv;
            d2 = ((numtyp)3.0*d1 + pre2*expm2) * r2inv;
            d3 = ((numtyp)5.0*d2 + pre3*expm2) * r2inv;

            g0d1_g1d2_g2d3 = g0*d1 + g1*d2 + g2*d3;
            fax = delx * g0d1_g1d2_g2d3 -
              d1 * (qtmp*muj.x - qj*mui.x) +
              d2 * (pjdotr*mui.x + pidotr*muj.x);
            fay = dely * g0d1_g1d2_g2d3 -
              d1 * (qtmp*muj.y - qj*mui.y) +
              d2 * (pjdotr*mui.y + pidotr*muj.y);
            faz = delz * g0d1_g1d2_g2d3 -
              d1 * (qtmp*muj.z - qj*mui.z) +
              d2 * (pjdotr*mui.z + pidotr*muj.z);

            zaix = delx * (qj*d1 + d2*pjdotr) - d1*muj.x;
            zaiy = dely * (qj*d1 + d2*pjdotr) - d1*muj.y;
            zaiz = delz * (qj*d1 + d2*pjdotr) - d1*muj.z;
            zajx = delx * (-qtmp*d1 + d2*pidotr) - d1*mui.x;
            zajy = dely * (-qtmp*d1 + d2*pidotr) - d1*mui.y;
            zajz = delz * (-qtmp*d1 + d2*pidotr) - d1*mui.z;

            if (factor_coul > (numtyp)0.0) {
              facm1 = (numtyp)1.0 - factor_coul;
              fax *= facm1;
              fay *= facm1;
              faz *= facm1;
              zaix *= facm1;
              zaiy *= facm1;
              zaiz *= facm1;
              zajx *= facm1;
              zajy *= facm1;
              zajz *= facm1;
            }
          } else {
            fax = fay = faz = (numtyp)0.0;
            zaix = zaiy = zaiz = (numtyp)0.0;
            zajx = zajy = zajz = (numtyp)0.0;
          }

          forcecoul.x = fdx + fax;
          forcecoul.y = fdy + fay;
          forcecoul.z = fdz + faz;

          ticoul.x = mui.y*(zdiz + zaiz) - mui.z*(zdiy + zaiy);
          ticoul.y = mui.z*(zdix + zaix) - mui.x*(zdiz + zaiz);
          ticoul.z = mui.x*(zdiy + zaiy) - mui.y*(zdix + zaix);

        } else {
          forcecoul.x = forcecoul.y = forcecoul.z = (numtyp)0.0;
          ticoul.x = ticoul.y = ticoul.z = (numtyp)0.0;
        }

        force.x = qqrd2e*forcecoul.x + delx*force_lj;
        force.y = qqrd2e*forcecoul.y + dely*force_lj;
        force.z = qqrd2e*forcecoul.z + delz*force_lj;
        f.x+=force.x;
        f.y+=force.y;
        f.z+=force.z;
        tor.x+=qqrd2e*ticoul.x;
        tor.y+=qqrd2e*ticoul.y;
        tor.z+=qqrd2e*ticoul.z;

        if (eflag>0) {
          acctyp e = (acctyp)0.0;
          if (rsq < cut_coulsq && factor_coul > (numtyp)0.0) {
            e = qqrd2e*(b0*g0 + b1*g1 + b2*g2);
            if (factor_coul < (numtyp)1.0) {
              e_coul *= factor_coul;
              e_coul += ((numtyp)1.0-factor_coul) * qqrd2e * (d0*g0 + d1*g1 + d2*g2);
            }
          } else e = (acctyp)0.0;
          e_coul += e;

          if (rsq < lj1[mtype].z) {
            e=r6inv*(lj3[mtype].x*r6inv-lj3[mtype].y);
            energy+=factor_lj*(e-lj3[mtype].z);
          }
        }
        if (vflag>0) {
          virial[0] += delx*force.x;
          virial[1] += dely*force.y;
          virial[2] += delz*force.z;
          virial[3] += delx*force.y;
          virial[4] += delx*force.z;
          virial[5] += dely*force.z;
        }
      }

    } // for nbor
    store_answers_tq(f,tor,energy,e_coul,virial,ii,inum,tid,t_per_atom,offset,eflag,
                    vflag,ans,engv);
  } // if ii
}

