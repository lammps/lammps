// **************************************************************************
//                                  morse.cu
//                             -------------------
//                           W. Michael Brown (ORNL)
//
//  Device code for acceleration of the morse pair style
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                :
//    email                : brownw@ornl.gov
// ***************************************************************************/

#if defined(NV_KERNEL) || defined(USE_HIP)

#include "lal_aux_fun1.h"
#ifndef _DOUBLE_DOUBLE
_texture( pos_tex,float4);
#else
_texture_2d( pos_tex,int4);
#endif

#else
#define pos_tex x_
#endif

__kernel void k_morse(const __global numtyp4 *restrict x_,
                      const __global numtyp4 *restrict mor1,
                      const __global numtyp2 *restrict mor2,
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
      numtyp r = delx*delx+dely*dely+delz*delz;

      int mtype=itype*lj_types+jtype;
      if (r<mor1[mtype].x) {
        r=ucl_sqrt(r);
        numtyp dexp=r-mor1[mtype].z;
        dexp=ucl_exp(-mor1[mtype].w*dexp);
        numtyp dm=dexp*dexp-dexp;
        numtyp force = mor1[mtype].y*dm/r*factor_lj;

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (eflag>0) {
          numtyp e=mor2[mtype].x*(dexp*dexp - 2.0*dexp) - mor2[mtype].y;
          energy+=e*factor_lj;
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

__kernel void k_morse_fast(const __global numtyp4 *restrict x_,
                           const __global numtyp4 *restrict mor1_in,
                           const __global numtyp2 *restrict mor2_in,
                           const __global numtyp *restrict sp_lj_in,
                           const __global int *dev_nbor,
                           const __global int *dev_packed,
                           __global acctyp4 *restrict ans,
                           __global acctyp *restrict engv,
                           const int eflag, const int vflag, const int inum,
                           const int nbor_pitch, const int t_per_atom) {
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  __local numtyp4 mor1[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp2 mor2[MAX_SHARED_TYPES*MAX_SHARED_TYPES];
  __local numtyp sp_lj[4];
  if (tid<4)
    sp_lj[tid]=sp_lj_in[tid];
  if (tid<MAX_SHARED_TYPES*MAX_SHARED_TYPES) {
    mor1[tid]=mor1_in[tid];
    if (eflag>0)
      mor2[tid]=mor2_in[tid];
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
      numtyp r = delx*delx+dely*dely+delz*delz;

      if (r<mor1[mtype].x) {
        r=ucl_sqrt(r);
        numtyp dexp=r-mor1[mtype].z;
        dexp=ucl_exp(-mor1[mtype].w*dexp);
        numtyp dm=dexp*dexp-dexp;
        numtyp force = mor1[mtype].y*dm/r*factor_lj;

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (eflag>0) {
          numtyp e=mor2[mtype].x*(dm-dexp)-mor2[mtype].y;
          energy+=e*factor_lj;
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

