/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Mike Brown (ORNL), brownw@ornl.gov
------------------------------------------------------------------------- */

#ifndef CRML_GPU_KERNEL
#define CRML_GPU_KERNEL

#define MAX_BIO_SHARED_TYPES 128

#ifdef _DOUBLE_DOUBLE
#define numtyp double
#define numtyp2 double2
#define numtyp4 double4
#define acctyp double
#define acctyp4 double4
#endif

#ifdef _SINGLE_DOUBLE
#define numtyp float
#define numtyp2 float2
#define numtyp4 float4
#define acctyp double
#define acctyp4 double4
#endif

#ifndef numtyp
#define numtyp float
#define numtyp2 float2
#define numtyp4 float4
#define acctyp float
#define acctyp4 float4
#endif

#define EWALD_F (numtyp)1.12837917
#define EWALD_P (numtyp)0.3275911
#define A1 (numtyp)0.254829592
#define A2 (numtyp)-0.284496736
#define A3 (numtyp)1.421413741
#define A4 (numtyp)-1.453152027
#define A5 (numtyp)1.061405429

#ifdef NV_KERNEL

#include "geryon/ucl_nv_kernel.h"
texture<float4> pos_tex;
texture<float> q_tex;

#ifdef _DOUBLE_DOUBLE
__inline double4 fetch_pos(const int& i, const double4 *pos)
{
  return pos[i];
}
__inline double fetch_q(const int& i, const double *q)
{
  return q[i];
}
#else
__inline float4 fetch_pos(const int& i, const float4 *pos)
{
  return tex1Dfetch(pos_tex, i);
}
__inline float fetch_q(const int& i, const float *q)
{
  return tex1Dfetch(q_tex, i);
}
#endif

#else

#pragma OPENCL EXTENSION cl_khr_fp64: enable
#define GLOBAL_ID_X get_global_id(0)
#define THREAD_ID_X get_local_id(0)
#define BLOCK_ID_X get_group_id(0)
#define BLOCK_SIZE_X get_local_size(0)
#define __syncthreads() barrier(CLK_LOCAL_MEM_FENCE)
#define __inline inline

#define fetch_pos(i,y) x_[i]
#define fetch_q(i,y) q_[i]

#endif

#define SBBITS 30
#define NEIGHMASK 0x3FFFFFFF
__inline int sbmask(int j) { return j >> SBBITS & 3; }

__kernel void kernel_pair(__global numtyp4 *x_, __global numtyp4 *lj1,
                          const int lj_types, 
                          __global numtyp *sp_lj_in, __global int *dev_nbor, 
                          __global acctyp4 *ans, __global acctyp *engv, 
                          const int eflag, const int vflag, const int inum, 
                          const int nall, const int nbor_pitch,
                          __global numtyp *q_, const numtyp cut_coulsq,
                          const numtyp qqrd2e, const numtyp g_ewald,
                          const numtyp denom_lj, const numtyp cut_bothsq, 
                          const numtyp cut_ljsq, const numtyp cut_lj_innersq) {

  // ii indexes the two interacting particles in gi
  int ii=GLOBAL_ID_X;
  __local numtyp sp_lj[8];
  sp_lj[0]=sp_lj_in[0];
  sp_lj[1]=sp_lj_in[1];
  sp_lj[2]=sp_lj_in[2];
  sp_lj[3]=sp_lj_in[3];
  sp_lj[4]=sp_lj_in[4];
  sp_lj[5]=sp_lj_in[5];
  sp_lj[6]=sp_lj_in[6];
  sp_lj[7]=sp_lj_in[7];

  if (ii<inum) {
    acctyp energy=(acctyp)0;
    acctyp e_coul=(acctyp)0;
    acctyp4 f;
    f.x=(acctyp)0;
    f.y=(acctyp)0;
    f.z=(acctyp)0;
    acctyp virial[6];
    for (int i=0; i<6; i++)
      virial[i]=(acctyp)0;

    __global int *nbor=dev_nbor+ii;
    int i=*nbor;
    nbor+=nbor_pitch;
    int numj=*nbor;
    nbor+=nbor_pitch;
    __global int *list_end=nbor+mul24(numj,nbor_pitch);
  
    numtyp4 ix=fetch_pos(i,x_); //x_[i];
    numtyp qtmp=fetch_q(i,q_);
    int itype=ix.w;

    for ( ; nbor<list_end; nbor+=nbor_pitch) {
      int j=*nbor;

      numtyp factor_lj, factor_coul;
      factor_lj = sp_lj[sbmask(j)];
      factor_coul = (numtyp)1.0-sp_lj[sbmask(j)+4];
      j &= NEIGHMASK;

      numtyp4 jx=fetch_pos(j,x_); //x_[j];
      int jtype=jx.w;

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      int mtype=itype*lj_types+jtype;
      if (rsq<cut_bothsq) {
        numtyp r2inv=(numtyp)1.0/rsq;
        numtyp forcecoul, force_lj, force, r6inv, prefactor, _erfc, switch1;

        if (rsq < cut_ljsq) {
          r6inv = r2inv*r2inv*r2inv;
          force_lj = factor_lj*r6inv*(lj1[mtype].x*r6inv-lj1[mtype].y);
          if (rsq > cut_lj_innersq) {
            switch1 = (cut_ljsq-rsq);
            numtyp switch2 = (numtyp)12.0*rsq*switch1*(rsq-cut_lj_innersq)/ 
                             denom_lj;
            switch1 *= switch1;
            switch1 *= (cut_ljsq+(numtyp)2.0*rsq-(numtyp)3.0*cut_lj_innersq)/
                       denom_lj;
            switch2 *= r6inv*(lj1[mtype].z*r6inv-lj1[mtype].w);
            force_lj = force_lj*switch1+switch2;
          }
        } else
          force_lj = (numtyp)0.0;

        if (rsq < cut_coulsq) {
          numtyp r = sqrt(rsq);
          numtyp grij = g_ewald * r;
          numtyp expm2 = exp(-grij*grij);
          numtyp t = (numtyp)1.0 / ((numtyp)1.0 + EWALD_P*grij);
          _erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
          prefactor = qqrd2e * qtmp*fetch_q(j,q_)/r;
          forcecoul = prefactor * (_erfc + EWALD_F*grij*expm2-factor_coul);
        } else {
          forcecoul = (numtyp)0.0;
          prefactor = (numtyp)0.0;
        }

        force = (force_lj + forcecoul) * r2inv;

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (eflag>0) {
          e_coul += prefactor*(_erfc-factor_coul);
          if (rsq < cut_ljsq) {
            numtyp e=r6inv*(lj1[mtype].z*r6inv-lj1[mtype].w);
            if (rsq > cut_lj_innersq)
              e *= switch1;
            energy+=factor_lj*e;
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

    // Store answers
    __global acctyp *ap1=engv+ii;
    if (eflag>0) {
      *ap1=energy;
      ap1+=inum;
      *ap1=e_coul;
      ap1+=inum;
    }
    if (vflag>0) {
      for (int i=0; i<6; i++) {
        *ap1=virial[i];
        ap1+=inum;
      }
    }
    ans[ii]=f;
  } // if ii
}

__kernel void kernel_pair_fast(__global numtyp4 *x_, __global numtyp2 *ljd_in,
                               __global numtyp* sp_lj_in, __global int *dev_nbor, 
                               __global int *dev_list, __global acctyp4 *ans,
                               __global acctyp *engv, const int eflag,
                               const int vflag, const int inum, const int nall,
                               const int nbor_pitch, __global numtyp *q_,
                               const numtyp cut_coulsq, const numtyp qqrd2e,
                               const numtyp g_ewald, const numtyp denom_lj,
                               const numtyp cut_bothsq, const numtyp cut_ljsq,
                               const numtyp cut_lj_innersq,
                               const int t_per_atom) {
  int tid=THREAD_ID_X;
  __local numtyp2 ljd[MAX_BIO_SHARED_TYPES];
  __local numtyp sp_lj[8];

  acctyp energy;
  acctyp e_coul;
  acctyp4 f;
  acctyp virial[6];

  if (tid<8)
    sp_lj[tid]=sp_lj_in[tid];
  ljd[tid]=ljd_in[tid];
  ljd[tid+64]=ljd_in[tid+64];

  int ii=mul24((int)BLOCK_ID_X,(int)(BLOCK_SIZE_X)/t_per_atom);
  ii+=tid/t_per_atom;
  int offset=tid%t_per_atom;
  
  energy=(acctyp)0;
  e_coul=(acctyp)0;
  f.x=(acctyp)0;
  f.y=(acctyp)0;
  f.z=(acctyp)0;
  for (int o=0; o<6; o++)
    virial[o]=(acctyp)0;
  __syncthreads();
  
  if (ii<inum) {
    __global int *nbor=dev_nbor+ii;
    int i=*nbor;
    nbor+=nbor_pitch;
    int numj=*nbor;
    nbor+=nbor_pitch;

    int n_stride;
    __global int *list_end;
    if (dev_nbor==dev_list) {
      list_end=nbor+mul24(numj,nbor_pitch);
      nbor+=mul24(offset,nbor_pitch);
      n_stride=mul24(t_per_atom,nbor_pitch);
    } else {
      nbor=dev_list+*nbor;
      list_end=nbor+numj;
      n_stride=t_per_atom;
      nbor+=offset;
    }
  
    numtyp4 ix=fetch_pos(i,x_); //x_[i];
    numtyp qtmp=fetch_q(i,q_);
    int itype=ix.w;

    for ( ; nbor<list_end; nbor+=n_stride) {
      int j=*nbor;

      numtyp factor_lj, factor_coul;
      factor_lj = sp_lj[sbmask(j)];
      factor_coul = (numtyp)1.0-sp_lj[sbmask(j)+4];
      j &= NEIGHMASK;

      numtyp4 jx=fetch_pos(j,x_); //x_[j];
      int jtype=jx.w;

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      if (rsq<cut_bothsq) {
        numtyp r2inv=(numtyp)1.0/rsq;
        numtyp forcecoul, force_lj, force, prefactor, _erfc, switch1;
        numtyp lj3, lj4;

        if (rsq < cut_ljsq) {
          numtyp eps = sqrt(ljd[itype].x*ljd[jtype].x);
          numtyp sig6 = (numtyp)0.5 * (ljd[itype].y+ljd[jtype].y);

          numtyp sig_r_6 = sig6*sig6*r2inv;
          sig_r_6 = sig_r_6*sig_r_6*sig_r_6;
          lj4 = (numtyp)4.0*eps*sig_r_6;
          lj3 = lj4*sig_r_6;
          force_lj = factor_lj*((numtyp)12.0 * lj3 - (numtyp)6.0 * lj4);
          if (rsq > cut_lj_innersq) {
            switch1 = (cut_ljsq-rsq);
            numtyp switch2 = (numtyp)12.0*rsq*switch1*(rsq-cut_lj_innersq)/ 
                             denom_lj;
            switch1 *= switch1;
            switch1 *= (cut_ljsq+(numtyp)2.0*rsq-(numtyp)3.0*cut_lj_innersq)/
                       denom_lj;
            switch2 *= lj3-lj4;
            force_lj = force_lj*switch1+switch2;
          }
        } else
          force_lj = (numtyp)0.0;

        if (rsq < cut_coulsq) {
          numtyp r = sqrt(rsq);
          numtyp grij = g_ewald * r;
          numtyp expm2 = exp(-grij*grij);
          numtyp t = (numtyp)1.0 / ((numtyp)1.0 + EWALD_P*grij);
          _erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
          prefactor = qqrd2e * qtmp*fetch_q(j,q_)/r;
          forcecoul = prefactor * (_erfc + EWALD_F*grij*expm2-factor_coul);
        } else {
          forcecoul = (numtyp)0.0;
          prefactor = (numtyp)0.0;
        }

        force = (force_lj + forcecoul) * r2inv;

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (eflag>0) {
          e_coul += prefactor*(_erfc-factor_coul);
          if (rsq < cut_ljsq) {
            numtyp e=lj3-lj4;
            if (rsq > cut_lj_innersq)
              e *= switch1;
            energy+=factor_lj*e;
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
  } // if ii

  // Reduce answers
  if (t_per_atom>1) {
    __local acctyp s_energy[64];
    __local acctyp s_e_coul[64];
    __local acctyp4 s_f[64];
    __local acctyp s_virial[6][64];
    
    s_f[tid].x=f.x;
    s_f[tid].y=f.y;
    s_f[tid].z=f.z;
    s_energy[tid]=energy;
    s_e_coul[tid]=e_coul;
    for (int v=0; v<6; v++)
      s_virial[v][tid]=virial[v];

    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {
      if (offset < s) {
        s_f[tid].x += s_f[tid+s].x;
        s_f[tid].y += s_f[tid+s].y;
        s_f[tid].z += s_f[tid+s].z;
        s_energy[tid] += s_energy[tid+s];
        s_e_coul[tid] += s_e_coul[tid+s];
        s_virial[tid][0] += s_virial[tid+s][0];
        s_virial[tid][1] += s_virial[tid+s][1];
        s_virial[tid][2] += s_virial[tid+s][2];
        s_virial[tid][3] += s_virial[tid+s][3];
        s_virial[tid][4] += s_virial[tid+s][4];
        s_virial[tid][5] += s_virial[tid+s][5];
      }
    }
    
    f.x=s_f[tid].x;
    f.y=s_f[tid].y;
    f.z=s_f[tid].z;
    energy=s_energy[tid];
    e_coul=s_e_coul[tid];
    for (int v=0; v<6; v++)
      virial[v]=s_virial[v][tid];
  }

  // Store answers
  __global acctyp *ap1=engv+ii;
  if (ii<inum && offset==0) {
    if (eflag>0) {
      *ap1=energy;
      ap1+=inum;
      *ap1=e_coul;
      ap1+=inum;
    }
    if (vflag>0) {
      for (int v=0; v<6; v++) {
        *ap1=virial[v];
        ap1+=inum;
      }
    }
    ans[ii]=f;
  }
}

#endif
