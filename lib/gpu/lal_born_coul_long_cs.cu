// **************************************************************************
//                            born_coul_long_cs.cu
//                             -------------------
//                           Trung Dac Nguyen (Northwestern)
//
//  Device code for acceleration of the born/coul/long/cs pair style
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                : June 2018
//    email                : ndactrung@gmail.com
// ***************************************************************************/

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

// Note: EWALD_P is different from that in lal_preprocessor.h
//       acctyp is needed for these parameters
#define CS_EWALD_P (acctyp)9.95473818e-1
#define B0        (acctyp)-0.1335096380159268
#define B1        (acctyp)-2.57839507e-1
#define B2        (acctyp)-1.37203639e-1
#define B3        (acctyp)-8.88822059e-3
#define B4        (acctyp)-5.80844129e-3
#define B5        (acctyp)1.14652755e-1

#define EPSILON (acctyp)(1.0e-20)
#define EPS_EWALD (acctyp)(1.0e-6)
#define EPS_EWALD_SQR (acctyp)(1.0e-12)

__kernel void k_born_coul_long_cs(const __global numtyp4 *restrict x_,
                          const __global numtyp4 *restrict coeff1,
                          const __global numtyp4 *restrict coeff2,
                          const int lj_types,
                          const __global numtyp *restrict sp_lj_in,
                          const __global int *dev_nbor,
                          const __global int *dev_packed,
                          __global acctyp4 *restrict ans,
                          __global acctyp *restrict engv,
                          const int eflag, const int vflag, const int inum,
                          const int nbor_pitch,
                          const __global numtyp *restrict q_,
                          const __global numtyp4 *restrict cutsq_sigma,
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
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(acctyp)0;

  if (ii<inum) {
    int nbor, nbor_end;
    int i, numj;
    __local int n_stride;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    numtyp qtmp; fetch(qtmp,i,q_tex);
    int itype=ix.w;

    for ( ; nbor<nbor_end; nbor+=n_stride) {
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
      if (rsq<cutsq_sigma[mtype].x) { // cutsq 
        numtyp forcecoul,forceborn,force,r6inv,prefactor,_erfc,rexp;

        rsq += EPSILON; // Add Epsilon for case: r = 0; Interaction must be removed by special bond;
        numtyp r2inv = ucl_recip(rsq);

        if (rsq < cut_coulsq) {
          numtyp r = ucl_sqrt(rsq);
          fetch(prefactor,j,q_tex);
          prefactor *= qqrd2e * qtmp;
          if (factor_coul<(numtyp)1.0) {
            numtyp grij = g_ewald * (r+EPS_EWALD);
            numtyp expm2 = ucl_exp(-grij*grij);
            acctyp t = ucl_recip((numtyp)1.0 + CS_EWALD_P*grij);
            numtyp u = (numtyp)1.0 - t;
            _erfc = t * ((numtyp)1.0 + u*(B0+u*(B1+u*(B2+u*(B3+u*(B4+u*B5)))))) * expm2;
            prefactor /= (r+EPS_EWALD);
            forcecoul = prefactor * (_erfc + EWALD_F*grij*expm2 - ((numtyp)1.0-factor_coul));
            // Additionally r2inv needs to be accordingly modified since the later
            // scaling of the overall force shall be consistent
            r2inv = ucl_recip(rsq + EPS_EWALD_SQR);
          } else {
            numtyp grij = g_ewald * r;
            numtyp expm2 = ucl_exp(-grij*grij);
            acctyp t = ucl_recip((numtyp)1.0 + CS_EWALD_P*grij);
            numtyp u = (numtyp)1.0 - t;
            _erfc = t * ((numtyp)1.0 + u*(B0+u*(B1+u*(B2+u*(B3+u*(B4+u*B5)))))) * expm2;
            prefactor /= r;
            forcecoul = prefactor*(_erfc + EWALD_F*grij*expm2);
          }
        } else forcecoul = (numtyp)0.0;

        if (rsq < cutsq_sigma[mtype].y) { // cut_ljsq
          numtyp r = ucl_sqrt(rsq);
          rexp = ucl_exp((cutsq_sigma[mtype].z-r)*coeff1[mtype].x);
          r6inv = r2inv*r2inv*r2inv;
          forceborn = (coeff1[mtype].y*r*rexp - coeff1[mtype].z*r6inv
            + coeff1[mtype].w*r2inv*r6inv)*factor_lj;
        } else forceborn = (numtyp)0.0;

        force = (forcecoul + forceborn) * r2inv;

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (eflag>0) {
          if (rsq < cut_coulsq) {
            numtyp e = prefactor*_erfc;
            if (factor_coul<(numtyp)1.0) e -= ((numtyp)1.0-factor_coul)*prefactor;
            e_coul += e;
          }
          if (rsq < cutsq_sigma[mtype].y) {
            numtyp e=coeff2[mtype].x*rexp - coeff2[mtype].y*r6inv
              + coeff2[mtype].z*r2inv*r6inv;
            energy+=factor_lj*(e-coeff2[mtype].w);
          }
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
    store_answers_q(f,energy,e_coul,virial,ii,inum,tid,t_per_atom,offset,eflag,
                    vflag,ans,engv);
  } // if ii
}

__kernel void k_born_coul_long_cs_fast(const __global numtyp4 *restrict x_,
                               const __global numtyp4 *restrict coeff1_in,
                               const __global numtyp4 *restrict coeff2_in,
                               const __global numtyp *restrict sp_lj_in,
                               const __global int *dev_nbor,
                               const __global int *dev_packed,
                               __global acctyp4 *restrict ans,
                               __global acctyp *restrict engv,
                               const int eflag, const int vflag, const int inum,
                               const int nbor_pitch,
                               const __global numtyp *restrict q_,
                               const __global numtyp4 *restrict cutsq_sigma,
                               const numtyp cut_coulsq, const numtyp qqrd2e,
                               const numtyp g_ewald, const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  __local numtyp4 coeff1[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp4 coeff2[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp sp_lj[8];
  if (tid<8)
    sp_lj[tid]=sp_lj_in[tid];
  if (tid<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
    coeff1[tid]=coeff1_in[tid];
    if (eflag>0)
      coeff2[tid]=coeff2_in[tid];
  }

  acctyp energy=(acctyp)0;
  acctyp e_coul=(acctyp)0;
  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp virial[6];
  for (int i=0; i<6; i++)
    virial[i]=(acctyp)0;

  __syncthreads();

  if (ii<inum) {
    int nbor, nbor_end;
    int i, numj;
    __local int n_stride;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    numtyp qtmp; fetch(qtmp,i,q_tex);
    int iw=ix.w;
    int itype=fast_mul((int)MAX_SHARED_TYPES,iw);

    for ( ; nbor<nbor_end; nbor+=n_stride) {
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

      if (rsq<cutsq_sigma[mtype].x) { // cutsq 
        numtyp forcecoul,forceborn,force,r6inv,prefactor,_erfc,rexp;

        rsq += EPSILON; // Add Epsilon for case: r = 0; Interaction must be removed by special bond;
        numtyp r2inv = ucl_recip(rsq);

        if (rsq < cut_coulsq) {
          numtyp r = ucl_sqrt(rsq);
          fetch(prefactor,j,q_tex);
          prefactor *= qqrd2e * qtmp;
          if (factor_coul<(numtyp)1.0) {
            numtyp grij = g_ewald * (r+EPS_EWALD);
            numtyp expm2 = ucl_exp(-grij*grij);
            acctyp t = ucl_recip((numtyp)1.0 + CS_EWALD_P*grij);
            numtyp u = (numtyp)1.0 - t;
            _erfc = t * ((numtyp)1.0 + u*(B0+u*(B1+u*(B2+u*(B3+u*(B4+u*B5)))))) * expm2;
            prefactor /= (r+EPS_EWALD);
            forcecoul = prefactor * (_erfc + EWALD_F*grij*expm2 - ((numtyp)1.0-factor_coul));
            // Additionally r2inv needs to be accordingly modified since the later
            // scaling of the overall force shall be consistent
            r2inv = ucl_recip(rsq + EPS_EWALD_SQR);
          } else {
            numtyp grij = g_ewald * r;
            numtyp expm2 = ucl_exp(-grij*grij);
            acctyp t = ucl_recip((numtyp)1.0 + CS_EWALD_P*grij);
            numtyp u = (numtyp)1.0 - t;
            _erfc = t * ((numtyp)1.0 + u*(B0+u*(B1+u*(B2+u*(B3+u*(B4+u*B5)))))) * expm2;
            prefactor /= r;
            forcecoul = prefactor*(_erfc + EWALD_F*grij*expm2);
          }
        } else forcecoul = (numtyp)0.0;

        if (rsq < cutsq_sigma[mtype].y) { // cut_ljsq
          numtyp r = ucl_sqrt(rsq);
          rexp = ucl_exp((cutsq_sigma[mtype].z-r)*coeff1[mtype].x);
          r6inv = r2inv*r2inv*r2inv;
          forceborn = (coeff1[mtype].y*r*rexp - coeff1[mtype].z*r6inv
            + coeff1[mtype].w*r2inv*r6inv)*factor_lj;
        } else forceborn = (numtyp)0.0;

        force = (forcecoul + forceborn) * r2inv;

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (eflag>0) {
          if (rsq < cut_coulsq) {
            numtyp e = prefactor*_erfc;
            if (factor_coul<(numtyp)1.0) e -= ((numtyp)1.0-factor_coul)*prefactor;
            e_coul += e;
          }
          if (rsq < cutsq_sigma[mtype].y) {
            numtyp e=coeff2[mtype].x*rexp - coeff2[mtype].y*r6inv
              + coeff2[mtype].z*r2inv*r6inv;
            energy+=factor_lj*(e-coeff2[mtype].w);
          }
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
    store_answers_q(f,energy,e_coul,virial,ii,inum,tid,t_per_atom,offset,eflag,
                    vflag,ans,engv);
  } // if ii
}

