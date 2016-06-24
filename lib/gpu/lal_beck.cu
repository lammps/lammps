// **************************************************************************
//                                   beck.cu
//                             -------------------
//                           Trung Dac Nguyen (ORNL)
//
//  Device code for acceleration of the beck pair style
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                :
//    email                : nguyentd@ornl.gov
// ***************************************************************************/

#ifdef NV_KERNEL
#include "lal_aux_fun1.h"
#ifndef _DOUBLE_DOUBLE
texture<float4> pos_tex;
#else
texture<int4,1> pos_tex;
#endif
#else
#define pos_tex x_
#endif

__kernel void k_beck(const __global numtyp4 *restrict x_,
                     const __global numtyp4 *restrict beck1,
                     const __global numtyp4 *restrict beck2,
                     const int lj_types,
                     const __global numtyp *restrict sp_lj_in,
                     const __global int *dev_nbor,
                     const __global int *dev_packed,
                     __global acctyp4 *restrict ans,
                     __global acctyp *restrict engv,
                     const int eflag, const int vflag, const int inum,
                     const int nbor_pitch, const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  __local numtyp sp_lj[4];
  sp_lj[0]=sp_lj_in[0];
  sp_lj[1]=sp_lj_in[1];
  sp_lj[2]=sp_lj_in[2];
  sp_lj[3]=sp_lj_in[3];

  acctyp energy=(acctyp)0;
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
    int itype=ix.w;

    numtyp factor_lj;
    for ( ; nbor<nbor_end; nbor+=n_stride) {

      int j=dev_packed[nbor];
      factor_lj = sp_lj[sbmask(j)];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int jtype=jx.w;

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      int mtype=itype*lj_types+jtype;
      if (rsq<beck2[mtype].z) {
        numtyp r = ucl_sqrt(rsq);
        numtyp r5 = rsq*rsq*r;
        numtyp aaij = beck1[mtype].x;
        numtyp alphaij = beck1[mtype].y;
        numtyp betaij = beck1[mtype].z;
        numtyp term1 = aaij*aaij + rsq;
        numtyp term2 = pow(term1,(numtyp)-5.0);
        numtyp term3 = (numtyp)21.672 + (numtyp)30.0*aaij*aaij + (numtyp)6.0*rsq;
        numtyp term4 = alphaij + r5*betaij;
        numtyp term5 = alphaij + (numtyp)6.0*r5*betaij;
        numtyp rinv  = ucl_recip(r);
        numtyp force = beck2[mtype].x*ucl_exp((numtyp)-1.0*r*term4)*term5;
        force -= beck2[mtype].y*r*term2*term3;
        force *= factor_lj*rinv;

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (eflag>0) {
          numtyp term6 = pow(term1,(numtyp)-3);
          numtyp term1inv = ucl_recip(term1);
          numtyp e = beck2[mtype].x*ucl_exp((numtyp)-1.0*r*term4);
          e -= beck2[mtype].y*term6*((numtyp)1.0+((numtyp)2.709+(numtyp)3.0*aaij*aaij)*term1inv);
          energy+=factor_lj*e;
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
    store_answers(f,energy,virial,ii,inum,tid,t_per_atom,offset,eflag,vflag,
                  ans,engv);
  } // if ii
}

__kernel void k_beck_fast(const __global numtyp4 *restrict x_,
                          const __global numtyp4 *restrict beck1_in,
                          const __global numtyp4 *restrict beck2_in,
                          const __global numtyp *restrict sp_lj_in,
                          const __global int *dev_nbor,
                          const __global int *dev_packed,
                          __global acctyp4 *restrict ans,
                          __global acctyp *restrict engv,
                          const int eflag, const int vflag, const int inum,
                          const int nbor_pitch, const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  __local numtyp4 beck1[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp4 beck2[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp sp_lj[4];
  if (tid<4)
    sp_lj[tid]=sp_lj_in[tid];
  if (tid<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
    beck1[tid]=beck1_in[tid];
    beck2[tid]=beck2_in[tid];
  }

  acctyp energy=(acctyp)0;
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
    int iw=ix.w;
    int itype=fast_mul((int)MAX_SHARED_TYPES,iw);

    numtyp factor_lj;
    for ( ; nbor<nbor_end; nbor+=n_stride) {

      int j=dev_packed[nbor];
      factor_lj = sp_lj[sbmask(j)];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int mtype=itype+jx.w;

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      if (rsq<beck2[mtype].z) {
        numtyp r = ucl_sqrt(rsq);
        numtyp r5 = rsq*rsq*r;
        numtyp aaij = beck1[mtype].x;
        numtyp alphaij = beck1[mtype].y;
        numtyp betaij = beck1[mtype].z;
        numtyp term1 = aaij*aaij + rsq;
        numtyp term2 = pow(term1,(numtyp)-5.0);
        numtyp term3 = (numtyp)21.672 + (numtyp)30.0*aaij*aaij + (numtyp)6.0*rsq;
        numtyp term4 = alphaij + r5*betaij;
        numtyp term5 = alphaij + (numtyp)6.0*r5*betaij;
        numtyp rinv  = ucl_recip(r);
        numtyp force = beck2[mtype].x*ucl_exp((numtyp)-1.0*r*term4)*term5;
        force -= beck2[mtype].y*r*term2*term3;
        force *= factor_lj*rinv;

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (eflag>0) {
          numtyp term6 = pow(term1,(numtyp)-3);
          numtyp term1inv = ucl_recip(term1);
          numtyp e = beck2[mtype].x*ucl_exp((numtyp)-1.0*r*term4);
          e -= beck2[mtype].y*term6*((numtyp)1.0+((numtyp)2.709+(numtyp)3.0*aaij*aaij)*term1inv);
          energy+=factor_lj*e;
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
    store_answers(f,energy,virial,ii,inum,tid,t_per_atom,offset,eflag,vflag,
                  ans,engv);
  } // if ii
}

