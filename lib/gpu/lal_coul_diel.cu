// **************************************************************************
//                                 coul_diel.cu
//                             -------------------
//
//  Device code for acceleration of the coul/diel pair style
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

__kernel void k_coul_diel(const __global numtyp4 *restrict x_,
                          const __global numtyp4 *restrict cutsq,
                          const int lj_types,
                          const __global numtyp *restrict sp_lj_in,
                          const __global int *dev_nbor,
                          const __global int *dev_packed,
                          __global acctyp3 *restrict ans,
                          __global acctyp *restrict engv,
                          const int eflag, const int vflag, const int inum,
                          const int nbor_pitch,
                          const __global numtyp *restrict q_,
                          const numtyp qqrd2e, const numtyp eps_s, const numtyp a_eps,
                          const numtyp b_eps, const int t_per_atom) {
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

  // the dielectric at short distance (5.2)
  const numtyp eps_min = a_eps - b_eps;

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
      numtyp factor_coul;
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
        // Coulomb with a distance dependent dielectric of sigmoidal shape
        numtyp qj; fetch(qj,j,q_tex);
        numtyp r = ucl_sqrt(rsq);
        numtyp rarg = (r-cutsq[mtype].y)/cutsq[mtype].z;
        // epsr = a_eps + b_eps*tanh(rarg) goes from eps_min = a_eps - b_eps to
        // eps_s = a_eps + b_eps.  Written that way, epsr at short distance and eps_s - epsr
        // at long distance are differences of nearly equal numbers, which single precision
        // cannot resolve (and with --use_fast_math, nvcc computes a single precision tanh()
        // with an approximate instruction on sm_75 and newer).  So epsr, deps = eps_s - epsr
        // and 1 - tanh^2(rarg) come from e2 = exp(-2|rarg|), without such differences.
        numtyp e2 = ucl_exp((numtyp)-2.0*ucl_abs(rarg));
        numtyp inv = (numtyp)1.0/((numtyp)1.0 + e2);
        numtyp deps = (numtyp)2.0*b_eps*inv;
        numtyp epsr;
        if (rarg >= (numtyp)0.0) {
          deps *= e2;
          epsr = (eps_s + eps_min*e2)*inv;
        } else {
          epsr = (eps_min + eps_s*e2)*inv;
        }
        numtyp depsdr = (numtyp)4.0*b_eps*e2*inv*inv / cutsq[mtype].z;

        numtyp forcecoul = qqrd2e*qtmp*qj*(epsr*deps + eps_s*r*depsdr)/(epsr*epsr*rsq);
        numtyp force = factor_coul*forcecoul/r;

        if (EVFLAG && eflag) {
          numtyp e = qqrd2e*qtmp*qj*deps/(epsr*r);
          e_coul += factor_coul*e;
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

__kernel void k_coul_diel_fast(const __global numtyp4 *restrict x_,
                               const __global numtyp4 *restrict cutsq_in,
                               const __global numtyp *restrict sp_lj_in,
                               const __global int *dev_nbor,
                               const __global int *dev_packed,
                               __global acctyp3 *restrict ans,
                               __global acctyp *restrict engv,
                               const int eflag, const int vflag, const int inum,
                               const int nbor_pitch,
                               const __global numtyp *restrict q_,
                               const numtyp qqrd2e, const numtyp eps_s, const numtyp a_eps,
                               const numtyp b_eps, const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  __local numtyp4 cutsq[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp sp_lj[8];
  int n_stride;
  local_allocate_store_charge();

  if (tid<8)
    sp_lj[tid]=sp_lj_in[tid];
  if (tid<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
    cutsq[tid]=cutsq_in[tid];
  }

  // the dielectric at short distance (5.2)
  const numtyp eps_min = a_eps - b_eps;

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
      numtyp factor_coul;
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
        // Coulomb with a distance dependent dielectric of sigmoidal shape
        numtyp qj; fetch(qj,j,q_tex);
        numtyp r = ucl_sqrt(rsq);
        numtyp rarg = (r-cutsq[mtype].y)/cutsq[mtype].z;
        // epsr = a_eps + b_eps*tanh(rarg) goes from eps_min = a_eps - b_eps to
        // eps_s = a_eps + b_eps.  Written that way, epsr at short distance and eps_s - epsr
        // at long distance are differences of nearly equal numbers, which single precision
        // cannot resolve (and with --use_fast_math, nvcc computes a single precision tanh()
        // with an approximate instruction on sm_75 and newer).  So epsr, deps = eps_s - epsr
        // and 1 - tanh^2(rarg) come from e2 = exp(-2|rarg|), without such differences.
        numtyp e2 = ucl_exp((numtyp)-2.0*ucl_abs(rarg));
        numtyp inv = (numtyp)1.0/((numtyp)1.0 + e2);
        numtyp deps = (numtyp)2.0*b_eps*inv;
        numtyp epsr;
        if (rarg >= (numtyp)0.0) {
          deps *= e2;
          epsr = (eps_s + eps_min*e2)*inv;
        } else {
          epsr = (eps_min + eps_s*e2)*inv;
        }
        numtyp depsdr = (numtyp)4.0*b_eps*e2*inv*inv / cutsq[mtype].z;

        numtyp forcecoul = qqrd2e*qtmp*qj*(epsr*deps + eps_s*r*depsdr)/(epsr*epsr*rsq);
        numtyp force = factor_coul*forcecoul/r;

        if (EVFLAG && eflag) {
          numtyp e = qqrd2e*qtmp*qj*deps/(epsr*r);
          e_coul += factor_coul*e;
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
