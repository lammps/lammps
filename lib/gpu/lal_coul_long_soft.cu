// **************************************************************************
//                                 coul_long_soft.cu
//                             -------------------
//      Claude Opus 5.5 (Anthropic), under the direction of Axel Kohlmeyer
//
//  Device code for acceleration of the coul/long/soft pair style
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

__kernel void k_coul_long_soft(const __global numtyp4 *restrict x_,
                          const __global numtyp4 *restrict scl,
                          const int lj_types,
                          const __global numtyp *restrict sp_cl_in,
                          const __global int *dev_nbor,
                          const __global int *dev_packed,
                          __global acctyp3 *restrict ans,
                          __global acctyp *restrict engv,
                          const int eflag, const int vflag, const int inum,
                          const int nbor_pitch,
                          const __global numtyp *restrict q_,
                          const numtyp cut_coulsq, const numtyp qqrd2e,
                          const numtyp g_ewald, const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  __local numtyp sp_cl[4];
  int n_stride;
  local_allocate_store_charge();

  sp_cl[0]=sp_cl_in[0];
  sp_cl[1]=sp_cl_in[1];
  sp_cl[2]=sp_cl_in[2];
  sp_cl[3]=sp_cl_in[3];

  acctyp3 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp e_coul, virial[6];
  if (EVFLAG) {
    e_coul=(acctyp)0;
    for (int i=0; i<6; i++) virial[i]=(acctyp)0;
  }

  if (ii<inum) {
    int nbor, nbor_end;
    int i, numj;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int itype=ix.w;
    numtyp qtmp; fetch(qtmp,i,q_tex);

    for ( ; nbor<nbor_end; nbor+=n_stride) {
      ucl_prefetch(dev_packed+nbor+n_stride);
      int j=dev_packed[nbor];

      numtyp factor_coul;
      factor_coul = (numtyp)1.0-sp_cl[sbmask(j)];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int jtype=jx.w;

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      int mtype=itype*lj_types+jtype;
      if (rsq < cut_coulsq) {
        numtyp r2inv=ucl_recip(rsq);
        numtyp force, prefactor, _erfc;

        numtyp r = ucl_rsqrt(r2inv);
        numtyp grij = g_ewald * r;
        numtyp expm2 = ucl_exp(-grij*grij);
        numtyp t = ucl_recip((numtyp)1.0 + EWALD_P*grij);
        _erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;

        // soft-core Coulomb: scl.x = scale, scl.y = lambda^n, scl.z = alpha_c (1-lambda)^2.
        // the soft core replaces only the 1/r factor, while the Ewald damping keeps the
        // true distance, so the exponential term is scaled by denc^2/r^2
        const numtyp4 sclij = scl[mtype];
        const numtyp denc = ucl_sqrt(sclij.z + rsq);
        fetch(prefactor,j,q_tex);
        prefactor *= qqrd2e * sclij.x * sclij.y * qtmp/(denc*denc*denc);
        force = prefactor * (_erfc + EWALD_F*grij*expm2*denc*denc*r2inv-factor_coul);

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (EVFLAG && eflag) {
          e_coul += prefactor*denc*denc*(_erfc-factor_coul);
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
  acctyp energy;
  if (EVFLAG) energy=(acctyp)0.0;
  store_answers_q(f,energy,e_coul,virial,ii,inum,tid,t_per_atom,offset,eflag,
                  vflag,ans,engv);
}

__kernel void k_coul_long_soft_fast(const __global numtyp4 *restrict x_,
                               const __global numtyp4 *restrict scl_in,
                               const __global numtyp *restrict sp_cl_in,
                               const __global int *dev_nbor,
                               const __global int *dev_packed,
                               __global acctyp3 *restrict ans,
                               __global acctyp *restrict engv,
                               const int eflag, const int vflag, const int inum,
                               const int nbor_pitch,
                               const __global numtyp *restrict q_,
                               const numtyp cut_coulsq, const numtyp qqrd2e,
                               const numtyp g_ewald, const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  __local numtyp4 scl[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp sp_cl[4];
  int n_stride;
  local_allocate_store_charge();

  if (tid<4)
    sp_cl[tid]=sp_cl_in[tid];
  if (tid<MAX_SHARED_TYPES*MAX_SHARED_TYPES)
    scl[tid]=scl_in[tid];

  acctyp3 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp e_coul, virial[6];
  if (EVFLAG) {
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

      numtyp factor_coul;
      factor_coul = (numtyp)1.0-sp_cl[sbmask(j)];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int mtype=itype+jx.w;

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      if (rsq < cut_coulsq) {
        numtyp r2inv=ucl_recip(rsq);
        numtyp force, prefactor, _erfc;

        numtyp r = ucl_rsqrt(r2inv);
        numtyp grij = g_ewald * r;
        numtyp expm2 = ucl_exp(-grij*grij);
        numtyp t = ucl_recip((numtyp)1.0 + EWALD_P*grij);
        _erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;

        // soft-core Coulomb: scl.x = scale, scl.y = lambda^n, scl.z = alpha_c (1-lambda)^2.
        // the soft core replaces only the 1/r factor, while the Ewald damping keeps the
        // true distance, so the exponential term is scaled by denc^2/r^2
        const numtyp4 sclij = scl[mtype];
        const numtyp denc = ucl_sqrt(sclij.z + rsq);
        fetch(prefactor,j,q_tex);
        prefactor *= qqrd2e * sclij.x * sclij.y * qtmp/(denc*denc*denc);
        force = prefactor * (_erfc + EWALD_F*grij*expm2*denc*denc*r2inv-factor_coul);

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (EVFLAG && eflag) {
          e_coul += prefactor*denc*denc*(_erfc-factor_coul);
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
  acctyp energy;
  if (EVFLAG) energy=(acctyp)0.0;
  store_answers_q(f,energy,e_coul,virial,ii,inum,tid,t_per_atom,offset,eflag,
                   vflag,ans,engv);
}

