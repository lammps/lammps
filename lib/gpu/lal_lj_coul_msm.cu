// **************************************************************************
//                               lj_coul_msm.cu
//                             -------------------
//                           Trung Dac Nguyen (ORNL)
//
//  Device code for acceleration of the lj/cut/coul/msm pair style
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
_texture( gcons_tex,float);
_texture( dgcons_tex,float);
#else
_texture_2d( pos_tex,int4);
_texture( q_tex,int2);
_texture( gcons_tex,int2);
_texture( dgcons_tex,int2);
#endif

#if (__CUDACC_VER_MAJOR__ >= 11)
#define gcons_tex gcons
#define dgcons_tex dgcons
#endif

#else
#define pos_tex x_
#define q_tex q_
#define gcons_tex gcons
#define dgcons_tex dgcons
#endif

/* ----------------------------------------------------------------------
   compute gamma for MSM and pair styles
   see Eq 4 from Parallel Computing 35 (2009) 164177
------------------------------------------------------------------------- */

ucl_inline numtyp gamma(const numtyp rho, const int order,
                        const __global numtyp *gcons) {
  if (rho <= (numtyp)1.0) {
    const int split_order = order/2;
    const numtyp rho2 = rho*rho;
    numtyp g; fetch(g,7*split_order+0,gcons_tex);
    numtyp rho_n = rho2;
    for (int n=1; n<=split_order; n++) {
      numtyp tmp; fetch(tmp,7*split_order+n,gcons_tex);
      g += tmp*rho_n;
      rho_n *= rho2;
    }
    return g;
  } else
    return ((numtyp)1.0/rho);
}

/* ----------------------------------------------------------------------
   compute the derivative of gamma for MSM and pair styles
   see Eq 4 from Parallel Computing 35 (2009) 164-177
------------------------------------------------------------------------- */

ucl_inline numtyp dgamma(const numtyp rho, const int order,
                         const __global numtyp *dgcons) {
  if (rho <= (numtyp)1.0) {
    const int split_order = order/2;
    const numtyp rho2 = rho*rho;
    numtyp dg; fetch(dg,6*split_order+0,dgcons_tex);
    dg *= rho;
    numtyp rho_n = rho*rho2;
    for (int n=1; n<split_order; n++) {
      numtyp tmp; fetch(tmp,6*split_order+n,dgcons_tex);
      dg += tmp*rho_n;
      rho_n *= rho2;
    }
    return dg;
  } else
    return ((numtyp)-1.0/rho/rho);
}

__kernel void k_lj_coul_msm(const __global numtyp4 *restrict x_,
                             const __global numtyp4 *restrict lj1,
                             const __global numtyp4 *restrict lj3,
                             const __global numtyp *restrict gcons,
                             const __global numtyp *restrict dgcons,
                             const int lj_types,
                             const __global numtyp *restrict sp_lj_in,
                             const __global int *dev_nbor,
                             const __global int *dev_packed,
                             __global acctyp3 *restrict ans,
                             __global acctyp *restrict engv,
                             const int eflag, const int vflag, const int inum,
                             const int nbor_pitch,
                             const __global numtyp *restrict q_,
                             const numtyp cut_coulsq, const numtyp qqrd2e,
                             const int order, const int t_per_atom) {
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

    numtyp cut_coul = ucl_sqrt(cut_coulsq);

    for ( ; nbor<nbor_end; nbor+=n_stride) {
      ucl_prefetch(dev_packed+nbor+n_stride);
      int j=dev_packed[nbor];

      numtyp factor_lj, factor_coul;
      factor_lj = sp_lj[sbmask(j)];
      factor_coul = (numtyp)1.0-sp_lj[sbmask(j)+4];
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
        numtyp forcecoul, force_lj, force, r6inv, prefactor, egamma, fgamma;

        if (rsq < lj1[mtype].w) {
          r6inv = r2inv*r2inv*r2inv;
          force_lj = factor_lj*r6inv*(lj1[mtype].x*r6inv-lj1[mtype].y);
        } else
          force_lj = (numtyp)0.0;

        if (rsq < cut_coulsq) {
          numtyp r = ucl_rsqrt(r2inv);
          fetch(prefactor,j,q_tex);
          prefactor *= qqrd2e * qtmp/r;
          numtyp rho = r/cut_coul;
          egamma = (numtyp)1.0 - rho*gamma(rho, order, gcons);
          fgamma = (numtyp)1.0 + (rsq/cut_coulsq)*dgamma(rho, order, dgcons);
          forcecoul = prefactor * fgamma;
        } else
          forcecoul = (numtyp)0.0;

        force = (force_lj + forcecoul) * r2inv;

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (EVFLAG && eflag) {
          if (rsq < cut_coulsq)
            e_coul += prefactor*(egamma-factor_coul);
          if (rsq < lj1[mtype].w) {
            numtyp e=r6inv*(lj3[mtype].x*r6inv-lj3[mtype].y);
            energy+=factor_lj*(e-lj3[mtype].z);
          }
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

__kernel void k_lj_coul_msm_fast(const __global numtyp4 *restrict x_,
                                 const __global numtyp4 *restrict lj1_in,
                                 const __global numtyp4 *restrict lj3_in,
                                 const __global numtyp *restrict gcons,
                                 const __global numtyp *restrict dgcons,
                                 const __global numtyp *restrict sp_lj_in,
                                 const __global int *dev_nbor,
                                 const __global int *dev_packed,
                                 __global acctyp3 *restrict ans,
                                 __global acctyp *restrict engv,
                                 const int eflag, const int vflag,
                                 const int inum,  const int nbor_pitch,
                                 const __global numtyp *restrict q_,
                                 const numtyp cut_coulsq, const numtyp qqrd2e,
                                 const int order, const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  __local numtyp4 lj1[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp4 lj3[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp sp_lj[8];
  int n_stride;
  local_allocate_store_charge();

  if (tid<8)
    sp_lj[tid]=sp_lj_in[tid];
  if (tid<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
    lj1[tid]=lj1_in[tid];
    if (EVFLAG && eflag)
      lj3[tid]=lj3_in[tid];
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

    numtyp cut_coul = ucl_sqrt(cut_coulsq);

    for ( ; nbor<nbor_end; nbor+=n_stride) {
      ucl_prefetch(dev_packed+nbor+n_stride);
      int j=dev_packed[nbor];

      numtyp factor_lj, factor_coul;
      factor_lj = sp_lj[sbmask(j)];
      factor_coul = (numtyp)1.0-sp_lj[sbmask(j)+4];
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
        numtyp forcecoul, force_lj, force, r6inv, prefactor, egamma, fgamma;

        if (rsq < lj1[mtype].w) {
          r6inv = r2inv*r2inv*r2inv;
          force_lj = factor_lj*r6inv*(lj1[mtype].x*r6inv-lj1[mtype].y);
        } else
          force_lj = (numtyp)0.0;

        if (rsq < cut_coulsq) {
          numtyp r = ucl_rsqrt(r2inv);
          fetch(prefactor,j,q_tex);
          prefactor *= qqrd2e * qtmp/r;
          numtyp rho = r/cut_coul;
          egamma = (numtyp)1.0 - rho*gamma(rho, order, gcons);
          fgamma = (numtyp)1.0 + (rsq/cut_coulsq)*dgamma(rho, order, dgcons);
          forcecoul = prefactor * fgamma;
        } else
          forcecoul = (numtyp)0.0;

        force = (force_lj + forcecoul) * r2inv;

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (EVFLAG && eflag) {
          if (rsq < cut_coulsq)
            e_coul += prefactor*(egamma-factor_coul);
          if (rsq < lj1[mtype].w) {
            numtyp e=r6inv*(lj3[mtype].x*r6inv-lj3[mtype].y);
            energy+=factor_lj*(e-lj3[mtype].z);
          }
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

