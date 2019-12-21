// **************************************************************************
//                                   zbl.cu
//                             -------------------
//                              Trung Dac Nguyen
//
//  Device code for acceleration of the zbl pair style
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                :
//    email                : ndactrung@gmail.com
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

// ZBL constants

#define c1 (numtyp)0.02817
#define c2 (numtyp)0.28022
#define c3 (numtyp)0.50986
#define c4 (numtyp)0.18175

/* ----------------------------------------------------------------------
   compute ZBL pair energy
------------------------------------------------------------------------- */

ucl_inline numtyp e_zbl(numtyp r, numtyp d1aij, numtyp d2aij,
                      numtyp d3aij, numtyp d4aij, numtyp zzeij) {

  numtyp rinv = ucl_recip(r);

  numtyp sum = c1*ucl_exp(-d1aij*r);
  sum += c2*ucl_exp(-d2aij*r);
  sum += c3*ucl_exp(-d3aij*r);
  sum += c4*ucl_exp(-d4aij*r);

  numtyp result = zzeij*sum*rinv;

  return result;
};

/* ----------------------------------------------------------------------
   compute ZBL first derivative
------------------------------------------------------------------------- */

ucl_inline numtyp dzbldr(numtyp r, numtyp d1aij, numtyp d2aij,
                         numtyp d3aij, numtyp d4aij, numtyp zzeij) {
  numtyp rinv = ucl_recip(r);

  numtyp e1 = ucl_exp(-d1aij*r);
  numtyp e2 = ucl_exp(-d2aij*r);
  numtyp e3 = ucl_exp(-d3aij*r);
  numtyp e4 = ucl_exp(-d4aij*r);

  numtyp sum = c1*e1;
  sum += c2*e2;
  sum += c3*e3;
  sum += c4*e4;

  numtyp sum_p = -c1*d1aij*e1;
  sum_p -= c2*d2aij*e2;
  sum_p -= c3*d3aij*e3;
  sum_p -= c4*d4aij*e4;

  numtyp result = zzeij*(sum_p - sum*rinv)*rinv;

  return result;
};

__kernel void k_zbl(const __global numtyp4 *restrict x_,
                    const __global numtyp4 *restrict coeff1,
                    const __global numtyp4 *restrict coeff2,
                    const __global numtyp4 *restrict coeff3,
                    const double cut_globalsq,
                    const double cut_innersq,
                    const double cut_inner,
                    const int lj_types,
                    const __global int *dev_nbor,
                    const __global int *dev_packed,
                    __global acctyp4 *restrict ans,
                    __global acctyp *restrict engv,
                    const int eflag, const int vflag, const int inum,
                    const int nbor_pitch, const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

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

    for ( ; nbor<nbor_end; nbor+=n_stride) {

      int j=dev_packed[nbor];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int jtype=jx.w;

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      int mtype=itype*lj_types+jtype;
      if (rsq<cut_globalsq) {
        numtyp r, t, force;
        r = ucl_sqrt(rsq);
        force = dzbldr(r, coeff2[mtype].x, coeff2[mtype].y,
                       coeff2[mtype].z, coeff2[mtype].w, coeff1[mtype].z);
        if (rsq>cut_innersq) {
          t = r - cut_inner;
          force = t*t * (coeff1[mtype].x + coeff1[mtype].y*t);
        }
        force *= (numtyp)-1.0*ucl_recip(r);

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (eflag>0) {
          numtyp e=e_zbl(r, coeff2[mtype].x, coeff2[mtype].y,
                         coeff2[mtype].z, coeff2[mtype].w, coeff1[mtype].z);
          e += coeff3[mtype].z;
          if (rsq > cut_innersq) {
            e += t*t*t * (coeff3[mtype].x + coeff3[mtype].y*t);
          }
          energy+=e;
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

__kernel void k_zbl_fast(const __global numtyp4 *restrict x_,
                         const __global numtyp4 *restrict coeff1_in,
                         const __global numtyp4 *restrict coeff2_in,
                         const __global numtyp4 *restrict coeff3_in,
                         const double cut_globalsq,
                         const double cut_innersq,
                         const double cut_inner,
                         const __global int *dev_nbor,
                         const __global int *dev_packed,
                         __global acctyp4 *restrict ans,
                         __global acctyp *restrict engv,
                         const int eflag, const int vflag, const int inum,
                         const int nbor_pitch, const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  __local numtyp4 coeff1[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp4 coeff2[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp4 coeff3[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  if (tid<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
    coeff1[tid]=coeff1_in[tid];
    coeff2[tid]=coeff2_in[tid];
    coeff3[tid]=coeff3_in[tid];
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

    for ( ; nbor<nbor_end; nbor+=n_stride) {

      int j=dev_packed[nbor];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int mtype=itype+jx.w;

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      if (rsq<cut_globalsq) {
        numtyp r, t, force;
        r = ucl_sqrt(rsq);
        force = dzbldr(r, coeff2[mtype].x, coeff2[mtype].y,
                       coeff2[mtype].z, coeff2[mtype].w, coeff1[mtype].z);
        if (rsq>cut_innersq) {
          t = r - cut_inner;
          force += t*t * (coeff1[mtype].x + coeff1[mtype].y*t);
        }

        force *= (numtyp)-1.0*ucl_recip(r);

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (eflag>0) {
          numtyp e=e_zbl(r, coeff2[mtype].x, coeff2[mtype].y,
                         coeff2[mtype].z, coeff2[mtype].w, coeff1[mtype].z);
          e += coeff3[mtype].z;
          if (rsq > cut_innersq) {
            e += t*t*t * (coeff3[mtype].x + coeff3[mtype].y*t);
          }
          energy+=e;
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

