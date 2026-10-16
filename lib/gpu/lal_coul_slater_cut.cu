// **************************************************************************
//                                 coul_slater_cut.cu
//                             -------------------
//
//  Device code for acceleration of the coul/slater/cut pair style
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

__kernel void k_coul_slater_cut(const __global numtyp4 *restrict x_,
                                const __global numtyp *restrict scale,
                                const int lj_types,
                                const __global numtyp *restrict sp_cl_in,
                                const __global int *dev_nbor,
                                const __global int *dev_packed,
                                __global acctyp3 *restrict ans,
                                __global acctyp *restrict engv,
                                const int eflag, const int vflag, const int inum,
                                const int nbor_pitch,
                                const __global numtyp *restrict q_,
                                const __global numtyp *restrict cutsq,
                                const numtyp qqrd2e, const numtyp lamda,
                                const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  __local numtyp sp_cl[8];
  int n_stride;
  local_allocate_store_charge();

  sp_cl[0]=sp_cl_in[0];
  sp_cl[1]=sp_cl_in[1];
  sp_cl[2]=sp_cl_in[2];
  sp_cl[3]=sp_cl_in[3];

  acctyp3 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp energy, e_coul, virial[6];
  if (EVFLAG) {
    energy=(acctyp)0;
    e_coul=(acctyp)0;
    for (int i=0; i<6; i++) virial[i]=(acctyp)0;
  }

  if (ii<inum) {
    int i, numj, nbor, nbor_end;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    numtyp qtmp; fetch(qtmp,i,q_tex);
    int itype=ix.w;

    numtyp factor_coul;
    for ( ; nbor<nbor_end; nbor+=n_stride) {
      ucl_prefetch(dev_packed+nbor+n_stride);

      int j=dev_packed[nbor];
      factor_coul = sp_cl[sbmask(j)];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int jtype=jx.w;

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      int mtype=itype*lj_types+jtype;
      if (rsq<cutsq[mtype]) {
        numtyp r2inv=ucl_recip(rsq);
        numtyp rinv=ucl_rsqrt(rsq);
        numtyp rlamdainv=ucl_recip(rinv*lamda);
        numtyp exprlmdainv=ucl_exp((numtyp)-2.0*rlamdainv);
        numtyp prefactor, force;

        fetch(prefactor,j,q_tex);
        prefactor *= qqrd2e*scale[mtype]*qtmp*rinv*factor_coul;
        numtyp bracket_term = (numtyp)1.0 - exprlmdainv*((numtyp)1.0 +
          (numtyp)2.0*rlamdainv*((numtyp)1.0+rlamdainv));
        force = prefactor * bracket_term * r2inv;

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (EVFLAG && eflag) {
          e_coul += prefactor*((numtyp)1.0 - ((numtyp)1.0+rlamdainv)*exprlmdainv);
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

__kernel void k_coul_slater_cut_fast(const __global numtyp4 *restrict x_,
                                     const __global numtyp *restrict scale,
                                     const __global numtyp *restrict sp_cl_in,
                                     const __global int *dev_nbor,
                                     const __global int *dev_packed,
                                     __global acctyp3 *restrict ans,
                                     __global acctyp *restrict engv,
                                     const int eflag, const int vflag, const int inum,
                                     const int nbor_pitch,
                                     const __global numtyp *restrict q_,
                                     const __global numtyp *restrict _cutsq,
                                     const numtyp qqrd2e, const numtyp lamda,
                                     const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  __local numtyp cutsq[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp sp_cl[4];
  int n_stride;
  local_allocate_store_charge();

  if (tid<4)
    sp_cl[tid]=sp_cl_in[tid];
  if (tid<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
    cutsq[tid]=_cutsq[tid];
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
    int i, numj, nbor, nbor_end;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    numtyp qtmp; fetch(qtmp,i,q_tex);
    int iw=ix.w;
    int itype=fast_mul((int)MAX_SHARED_TYPES,iw);

    numtyp factor_coul;
    for ( ; nbor<nbor_end; nbor+=n_stride) {
      ucl_prefetch(dev_packed+nbor+n_stride);

      int j=dev_packed[nbor];
      factor_coul = sp_cl[sbmask(j)];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int mtype=itype+jx.w;

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      if (rsq<cutsq[mtype]) {
        numtyp r2inv=ucl_recip(rsq);
        numtyp rinv=ucl_rsqrt(rsq);
        numtyp rlamdainv=ucl_recip(rinv*lamda);
        numtyp exprlmdainv=ucl_exp((numtyp)-2.0*rlamdainv);
        numtyp prefactor, force;

        fetch(prefactor,j,q_tex);
        prefactor *= qqrd2e*scale[mtype]*qtmp*rinv*factor_coul;
        numtyp bracket_term = (numtyp)1.0 - exprlmdainv*((numtyp)1.0 +
          (numtyp)2.0*rlamdainv*((numtyp)1.0+rlamdainv));
        force = prefactor * bracket_term * r2inv;

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (EVFLAG && eflag) {
          e_coul += prefactor*((numtyp)1.0 - ((numtyp)1.0+rlamdainv)*exprlmdainv);
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

