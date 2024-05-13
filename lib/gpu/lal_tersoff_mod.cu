// **************************************************************************
//                               tersoff_mod.cu
//                             -------------------
//                              Trung Dac Nguyen
//
//  Device code for acceleration of the tersoff pair style
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//       begin                :
//       email                : ndactrung@gmail.com
// ***************************************************************************

#if defined(NV_KERNEL) || defined(USE_HIP)
#include "lal_tersoff_mod_extra.h"

#ifndef _DOUBLE_DOUBLE
_texture( pos_tex,float4);
#else
_texture_2d( pos_tex,int4);
#endif

#else
#define pos_tex x_
#endif

//#define THREE_CONCURRENT

#define TWOTHIRD (numtyp)0.66666666666666666667

#if (SHUFFLE_AVAIL == 0)

#define local_allocate_acc_zeta()                                           \
    __local acctyp red_acc[BLOCK_PAIR];

#define acc_zeta(z, tid, t_per_atom, offset)                                \
  if (t_per_atom>1) {                                                       \
    red_acc[tid]=z;                                                         \
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                         \
      simdsync();                                                           \
      if (offset < s) {                                                     \
        red_acc[tid] += red_acc[tid+s];                                     \
      }                                                                     \
    }                                                                       \
    z=red_acc[tid];                                                         \
  }

#define store_answers_p(f, energy, virial, ii, inum, tid, t_per_atom,       \
                        offset, eflag, vflag, ans, engv, ev_stride)         \
  if (t_per_atom>1) {                                                       \
    simd_reduce_add3(t_per_atom, red_acc, offset, tid, f.x, f.y, f.z);      \
    if (EVFLAG && (vflag==2 || eflag==2)) {                                 \
      if (eflag) {                                                          \
        simdsync();                                                         \
        simd_reduce_add1(t_per_atom, red_acc, offset, tid, energy);         \
      }                                                                     \
      if (vflag) {                                                          \
        simdsync();                                                         \
        simd_reduce_arr(6, t_per_atom, red_acc, offset, tid, virial);       \
      }                                                                     \
    }                                                                       \
  }                                                                         \
  if (offset==0 && ii<inum) {                                               \
    acctyp3 old=ans[ii];                                                    \
    old.x+=f.x;                                                             \
    old.y+=f.y;                                                             \
    old.z+=f.z;                                                             \
    ans[ii]=old;                                                            \
  }                                                                         \
  if (EVFLAG && (eflag || vflag)) {                                         \
    int ei=BLOCK_ID_X;                                                      \
    if (eflag!=2 && vflag!=2) {                                             \
      if (eflag) {                                                          \
        simdsync();                                                         \
        block_reduce_add1(simd_size(), red_acc, tid, energy);               \
        if (vflag) __syncthreads();                                         \
        if (tid==0) {                                                       \
          engv[ei]+=energy*(acctyp)0.5;                                     \
          ei+=ev_stride;                                                    \
        }                                                                   \
      }                                                                     \
      if (vflag) {                                                          \
        simdsync();                                                         \
        block_reduce_arr(6, simd_size(), red_acc, tid, virial);             \
        if (tid==0) {                                                       \
          for (int r=0; r<6; r++) {                                         \
            engv[ei]+=virial[r]*(acctyp)0.5;                                \
            ei+=ev_stride;                                                  \
          }                                                                 \
        }                                                                   \
      }                                                                     \
    } else if (offset==0 && ii<inum) {                                      \
      int ei=ii;                                                            \
      if (EVFLAG && eflag) {                                                \
        engv[ei]+=energy*(acctyp)0.5;                                       \
        ei+=inum;                                                           \
      }                                                                     \
      if (EVFLAG && vflag) {                                                \
        for (int i=0; i<6; i++) {                                           \
          engv[ei]+=virial[i]*(acctyp)0.5;                                  \
          ei+=inum;                                                         \
        }                                                                   \
      }                                                                     \
    }                                                                       \
  }

#else

#define local_allocate_acc_zeta()

#define acc_zeta(z, tid, t_per_atom, offset)                                \
  if (t_per_atom>1) {                                                       \
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                         \
      z += shfl_down(z, s, t_per_atom);                                     \
    }                                                                       \
  }

#if (EVFLAG == 1)

#define store_answers_p(f, energy, virial, ii, inum, tid, t_per_atom,       \
                        offset, eflag, vflag, ans, engv, ev_stride)         \
  if (t_per_atom>1) {                                                       \
    simd_reduce_add3(t_per_atom, f.x, f.y, f.z);                            \
    if (vflag==2 || eflag==2) {                                             \
      if (eflag)                                                            \
        simd_reduce_add1(t_per_atom,energy);                                \
      if (vflag)                                                            \
        simd_reduce_arr(6, t_per_atom,virial);                              \
    }                                                                       \
  }                                                                         \
  if (offset==0 && ii<inum) {                                               \
    acctyp3 old=ans[ii];                                                    \
    old.x+=f.x;                                                             \
    old.y+=f.y;                                                             \
    old.z+=f.z;                                                             \
    ans[ii]=old;                                                            \
  }                                                                         \
  if (eflag || vflag) {                                                     \
    if (eflag!=2 && vflag!=2) {                                             \
      const int vwidth = simd_size();                                       \
      const int voffset = tid & (simd_size() - 1);                          \
      const int bnum = tid/simd_size();                                     \
      int active_subgs = BLOCK_SIZE_X/simd_size();                          \
      for ( ; active_subgs > 1; active_subgs /= vwidth) {                   \
        if (active_subgs < BLOCK_SIZE_X/simd_size()) __syncthreads();       \
        if (bnum < active_subgs) {                                          \
          if (eflag) {                                                      \
            simd_reduce_add1(vwidth, energy);                               \
            if (voffset==0) red_acc[6][bnum] = energy;                      \
          }                                                                 \
          if (vflag) {                                                      \
            simd_reduce_arr(6, vwidth, virial);                             \
            if (voffset==0)                                                 \
              for (int r=0; r<6; r++) red_acc[r][bnum]=virial[r];           \
          }                                                                 \
        }                                                                   \
                                                                            \
        __syncthreads();                                                    \
        if (tid < active_subgs) {                                           \
            if (eflag) energy = red_acc[6][tid];                            \
          if (vflag)                                                        \
            for (int r = 0; r < 6; r++) virial[r] = red_acc[r][tid];        \
        } else {                                                            \
          if (eflag) energy = (acctyp)0;                                    \
          if (vflag) for (int r = 0; r < 6; r++) virial[r] = (acctyp)0;     \
        }                                                                   \
      }                                                                     \
                                                                            \
      if (bnum == 0) {                                                      \
        int ei=BLOCK_ID_X;                                                  \
        if (eflag) {                                                        \
          simd_reduce_add1(vwidth, energy);                                 \
          if (tid==0) {                                                     \
            engv[ei]+=energy*(acctyp)0.5;                                   \
            ei+=ev_stride;                                                  \
          }                                                                 \
        }                                                                   \
        if (vflag) {                                                        \
          simd_reduce_arr(6, vwidth, virial);                               \
          if (tid==0) {                                                     \
            for (int r=0; r<6; r++) {                                       \
              engv[ei]+=virial[r]*(acctyp)0.5;                              \
              ei+=ev_stride;                                                \
            }                                                               \
          }                                                                 \
        }                                                                   \
      }                                                                     \
    } else if (offset==0 && ii<inum) {                                      \
      int ei=ii;                                                            \
      if (eflag) {                                                          \
        engv[ei]+=energy*(acctyp)0.5;                                       \
        ei+=inum;                                                           \
      }                                                                     \
      if (vflag) {                                                          \
        for (int i=0; i<6; i++) {                                           \
          engv[ei]+=virial[i]*(acctyp)0.5;                                  \
          ei+=inum;                                                         \
        }                                                                   \
      }                                                                     \
    }                                                                       \
  }

#else

#define store_answers_p(f, energy, virial, ii, inum, tid, t_per_atom,       \
                        offset, eflag, vflag, ans, engv, ev_stride)         \
  if (t_per_atom>1)                                                         \
    simd_reduce_add3(t_per_atom, f.x, f.y, f.z);                            \
  if (offset==0 && ii<inum) {                                               \
    acctyp3 old=ans[ii];                                                    \
    old.x+=f.x;                                                             \
    old.y+=f.y;                                                             \
    old.z+=f.z;                                                             \
    ans[ii]=old;                                                            \
  }

#endif
#endif

__kernel void k_tersoff_mod_short_nbor(const __global numtyp4 *restrict x_,
                                       const numtyp cutsq, const int ntypes,
                                       __global int * dev_nbor,
                                       const __global int * dev_packed,
                                       const int inum, const int nbor_pitch,
                                       const int t_per_atom) {
  const int ii=GLOBAL_ID_X;

  if (ii<inum) {
    const int i=dev_packed[ii];
    int nbor=ii+nbor_pitch;
    const int numj=dev_packed[nbor];
    nbor+=nbor_pitch;
    const int nbor_end=nbor+fast_mul(numj,nbor_pitch);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int newj=0;

    __global int *out_list=dev_nbor+2*nbor_pitch+ii*t_per_atom;
    const int out_stride=nbor_pitch*t_per_atom-t_per_atom;

    for ( ; nbor<nbor_end; nbor+=nbor_pitch) {
      int sj=dev_packed[nbor];
      int j = sj & NEIGHMASK;
      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      if (rsq<cutsq) {
        *out_list=sj;
        out_list++;
        newj++;
        if ((newj & (t_per_atom-1))==0)
          out_list+=out_stride;
      }
    } // for nbor
    dev_nbor[ii+nbor_pitch]=newj;
  } // if ii
}

// Tersoff is currently used for 3 elements at most: 3*3*3 = 27 entries
// while the block size should never be less than 32.
// SHARED_SIZE = 32 for now to reduce the pressure on the shared memory per block
// must be increased if there will be more than 3 elements in the future.

#define SHARED_SIZE 32

__kernel void k_tersoff_mod_zeta(const __global numtyp4 *restrict x_,
                             const __global numtyp4 *restrict ts1_in,
                             const __global numtyp4 *restrict ts2_in,
                             const __global numtyp4 *restrict ts3_in,
                             const __global numtyp4 *restrict ts4_in,
                             const __global numtyp4 *restrict ts5_in,
                             const __global int *restrict map,
                             const __global int *restrict elem2param,
                             const int nelements, const int nparams,
                             __global acctyp4 * zetaij,
                             const __global int * dev_nbor,
                             const int eflag, const int inum,
                             const int nbor_pitch, const int t_per_atom) {
  const int tpa_sq = fast_mul(t_per_atom,t_per_atom);

  int tid, ii, offset, n_stride;
  atom_info(tpa_sq,ii,tid,offset);

  local_allocate_acc_zeta();

  // must be increased if there will be more than 3 elements in the future.
  __local numtyp4 ts1[SHARED_SIZE];
  __local numtyp4 ts2[SHARED_SIZE];
  __local numtyp4 ts3[SHARED_SIZE];
  __local numtyp4 ts4[SHARED_SIZE];
  __local numtyp4 ts5[SHARED_SIZE];
  if (tid<nparams) {
    ts1[tid]=ts1_in[tid];
    ts2[tid]=ts2_in[tid];
    ts3[tid]=ts3_in[tid];
    ts4[tid]=ts4_in[tid];
    ts5[tid]=ts5_in[tid];
  }

  acctyp z = (acctyp)0;

  __syncthreads();

  if (ii<inum) {
    int nbor_j, nbor_end, i, numj;
    int offset_j=offset/t_per_atom;
    nbor_info_p(dev_nbor,nbor_pitch,t_per_atom,ii,offset_j,i,numj,
                n_stride,nbor_end,nbor_j);
    int offset_k=tid & (t_per_atom-1);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int itype=ix.w;
    itype=map[itype];

    int nborj_start = nbor_j;

    for ( ; nbor_j<nbor_end; nbor_j+=n_stride) {

      int j=dev_nbor[nbor_j];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int jtype=jx.w;
      jtype=map[jtype];
      int ijparam=elem2param[itype*nelements*nelements+jtype*nelements+jtype];

      // Compute rij
      numtyp4 delr1, delr2;
      delr1.x = jx.x-ix.x;
      delr1.y = jx.y-ix.y;
      delr1.z = jx.z-ix.z;
      numtyp rsq1 = delr1.x*delr1.x+delr1.y*delr1.y+delr1.z*delr1.z;

      // compute zeta_ij
      z = (acctyp)0;

      int nbor_k = nborj_start-offset_j+offset_k;
      for ( ; nbor_k < nbor_end; nbor_k+=n_stride) {
        int k=dev_nbor[nbor_k];
        k &= NEIGHMASK;

        if (k == j) continue;

        numtyp4 kx; fetch4(kx,k,pos_tex); //x_[k];
        int ktype=kx.w;
        ktype=map[ktype];
        int ijkparam=elem2param[itype*nelements*nelements+jtype*nelements+ktype];

        // Compute rik
        delr2.x = kx.x-ix.x;
        delr2.y = kx.y-ix.y;
        delr2.z = kx.z-ix.z;
        numtyp rsq2 = delr2.x*delr2.x+delr2.y*delr2.y+delr2.z*delr2.z;

        numtyp4 ts1_ijkparam = ts1[ijkparam]; //fetch4(ts1_ijkparam,ijkparam,ts1_tex);
        numtyp ijkparam_lam3 = ts1_ijkparam.z;
        numtyp ijkparam_powermint = ts1_ijkparam.w;
        numtyp4 ts2_ijkparam = ts2[ijkparam]; //fetch4(ts2_ijkparam,ijkparam,ts2_tex);
        numtyp ijkparam_bigr = ts2_ijkparam.z;
        numtyp ijkparam_bigd = ts2_ijkparam.w;
        numtyp4 ts4_ijkparam = ts4[ijkparam]; //fetch4(ts4_ijkparam,ijkparam,ts4_tex);
        numtyp ijkparam_c1 = ts4_ijkparam.x;
        numtyp ijkparam_c2 = ts4_ijkparam.y;
        numtyp ijkparam_c3 = ts4_ijkparam.z;
        numtyp ijkparam_c4 = ts4_ijkparam.w;
        numtyp4 ts5_ijkparam = ts5[ijkparam]; //fetch4(ts4_ijkparam,ijkparam,ts4_tex);
        numtyp ijkparam_c5 = ts5_ijkparam.x;
        numtyp ijkparam_h = ts5_ijkparam.y;
        z += zeta(ijkparam_powermint, ijkparam_lam3, ijkparam_bigr, ijkparam_bigd,
                  ijkparam_h, ijkparam_c1, ijkparam_c2, ijkparam_c3, ijkparam_c4,
                  ijkparam_c5, rsq1, rsq2, delr1, delr2);
      }

      acc_zeta(z, tid, t_per_atom, offset_k);

      numtyp4 ts1_ijparam = ts1[ijparam]; //fetch4(ts1_ijparam,ijparam,ts1_tex);
      numtyp ijparam_lam2 = ts1_ijparam.y;
      numtyp4 ts2_ijparam = ts2[ijparam]; //fetch4(ts2_ijparam,ijparam,ts2_tex);
      numtyp ijparam_bigb = ts2_ijparam.y;
      numtyp ijparam_bigr = ts2_ijparam.z;
      numtyp ijparam_bigd = ts2_ijparam.w;
      numtyp4 ts3_ijparam = ts3[ijparam]; //fetch4(ts3_ijparam,ijparam,ts3_tex);
      numtyp ijparam_beta = ts3_ijparam.x;
      numtyp ijparam_powern = ts3_ijparam.y;
      numtyp ijparam_powern_del = ts3_ijparam.z;
      numtyp ijparam_ca1 = ts3_ijparam.w;
      numtyp ijparam_ca4 = ucl_recip(ts3_ijparam.w);

      if (offset_k == 0) {
        numtyp fpfeng[4];
        force_zeta(ijparam_bigb, ijparam_bigr, ijparam_bigd, ijparam_lam2,
                   ijparam_beta, ijparam_powern, ijparam_powern_del, ijparam_ca1,
                   ijparam_ca4, rsq1, z, eflag, fpfeng);
        acctyp4 zij;
        zij.x = fpfeng[0];
        zij.y = fpfeng[1];
        zij.z = fpfeng[2];
        zij.w = z;
        zetaij[nbor_j-2*nbor_pitch] = zij;
      }

    } // for nbor
  } // if ii
}

__kernel void k_tersoff_mod_repulsive(const __global numtyp4 *restrict x_,
                                  const __global numtyp4 *restrict ts1_in,
                                  const __global numtyp4 *restrict ts2_in,
                                  const __global int *restrict map,
                                  const __global int *restrict elem2param,
                                  const int nelements, const int nparams,
                                  const __global int * dev_nbor,
                                  __global acctyp3 *restrict ans,
                                  __global acctyp *restrict engv,
                                  const int eflag, const int vflag,
                                  const int inum, const int nbor_pitch,
                                  const int t_per_atom, const int ev_stride) {
  int tid, ii, offset, n_stride;
  atom_info(t_per_atom,ii,tid,offset);

  local_allocate_store_pair();

  __local numtyp4 ts1[SHARED_SIZE];
  __local numtyp4 ts2[SHARED_SIZE];
  if (tid<nparams) {
    ts1[tid]=ts1_in[tid];
    ts2[tid]=ts2_in[tid];
  }

  acctyp3 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp energy, virial[6];
  if (EVFLAG) {
    energy=(acctyp)0;
    for (int i=0; i<6; i++) virial[i]=(acctyp)0;
  }

  __syncthreads();

  if (ii<inum) {
    int nbor, nbor_end, i, numj;
    nbor_info_p(dev_nbor,nbor_pitch,t_per_atom,ii,offset,i,numj,
                n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int itype=ix.w;
    itype=map[itype];

    for ( ; nbor<nbor_end; nbor+=n_stride) {

      int j=dev_nbor[nbor];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int jtype=jx.w;
      jtype=map[jtype];
      int ijparam=elem2param[itype*nelements*nelements+jtype*nelements+jtype];

      // Compute r12

      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      numtyp feng[2];
      numtyp ijparam_lam1 = ts1[ijparam].x;
      numtyp4 ts2_ijparam = ts2[ijparam];
      numtyp ijparam_biga = ts2_ijparam.x;
      numtyp ijparam_bigr = ts2_ijparam.z;
      numtyp ijparam_bigd = ts2_ijparam.w;

      repulsive(ijparam_bigr, ijparam_bigd, ijparam_lam1, ijparam_biga,
                rsq, eflag, feng);

      numtyp force = feng[0];
      f.x+=delx*force;
      f.y+=dely*force;
      f.z+=delz*force;

      if (EVFLAG && eflag)
        energy+=feng[1];
      if (EVFLAG && vflag) {
        virial[0] += delx*delx*force;
        virial[1] += dely*dely*force;
        virial[2] += delz*delz*force;
        virial[3] += delx*dely*force;
        virial[4] += delx*delz*force;
        virial[5] += dely*delz*force;
      }
    } // for nbor
  } // if ii
  store_answers_p(f,energy,virial,ii,inum,tid,t_per_atom,offset,eflag,vflag,
                  ans,engv,ev_stride);
}

__kernel void k_tersoff_mod_three_center(const __global numtyp4 *restrict x_,
                                     const __global numtyp4 *restrict ts1_in,
                                     const __global numtyp4 *restrict ts2_in,
                                     const __global numtyp4 *restrict ts4_in,
                                     const __global numtyp4 *restrict ts5_in,
                                     const __global int *restrict map,
                                     const __global int *restrict elem2param,
                                     const int nelements, const int nparams,
                                     const __global acctyp4 *restrict zetaij,
                                     const __global int * dev_nbor,
                                     __global acctyp3 *restrict ans,
                                     __global acctyp *restrict engv,
                                     const int eflag, const int vflag,
                                     const int inum,  const int nbor_pitch,
                                     const int t_per_atom, const int evatom) {
  const int tpa_sq=fast_mul(t_per_atom,t_per_atom);
  numtyp lam3, powermint, bigr, bigd, c1, c2, c3, c4, c5, h;

  int tid, ii, offset, n_stride;
  atom_info(tpa_sq,ii,tid,offset); // offset ranges from 0 to tpa_sq-1

  local_allocate_store_three();

  __local numtyp4 ts1[SHARED_SIZE];
  __local numtyp4 ts2[SHARED_SIZE];
  __local numtyp4 ts4[SHARED_SIZE];
  __local numtyp4 ts5[SHARED_SIZE];
  if (tid<nparams) {
    ts1[tid]=ts1_in[tid];
    ts2[tid]=ts2_in[tid];
    ts4[tid]=ts4_in[tid];
    ts5[tid]=ts5_in[tid];
  }

  acctyp3 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp energy, virial[6];
  if (EVFLAG) {
    energy=(acctyp)0;
    for (int i=0; i<6; i++) virial[i]=(acctyp)0;
  }
  numtyp tpainv = ucl_recip((numtyp)t_per_atom);

  __syncthreads();

  if (ii<inum) {
    int i, numj, nbor_j, nbor_end;
    int offset_j=offset/t_per_atom;
    nbor_info_p(dev_nbor,nbor_pitch,t_per_atom,ii,offset_j,i,numj,
                n_stride,nbor_end,nbor_j);
    int offset_k=tid & (t_per_atom-1);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int itype=ix.w;
    itype=map[itype];

    int nborj_start = nbor_j;
    for ( ; nbor_j<nbor_end; nbor_j+=n_stride) {

      int j=dev_nbor[nbor_j];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int jtype=jx.w;
      jtype=map[jtype];

      // Compute r12
      numtyp delr1[3];
      delr1[0] = jx.x-ix.x;
      delr1[1] = jx.y-ix.y;
      delr1[2] = jx.z-ix.z;
      numtyp rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];

      numtyp r1 = ucl_sqrt(rsq1);
      numtyp r1inv = ucl_rsqrt(rsq1);

      // look up for zeta_ij
      acctyp4 zeta_ij = zetaij[nbor_j-2*nbor_pitch];
      numtyp force = zeta_ij.x*tpainv;
      numtyp prefactor = zeta_ij.y;
      f.x += delr1[0]*force;
      f.y += delr1[1]*force;
      f.z += delr1[2]*force;

      if (EVFLAG && eflag) {
        energy+=zeta_ij.z*tpainv;
      }
      if (EVFLAG && vflag) {
        numtyp mforce = -force;
        virial[0] += delr1[0]*delr1[0]*mforce;
        virial[1] += delr1[1]*delr1[1]*mforce;
        virial[2] += delr1[2]*delr1[2]*mforce;
        virial[3] += delr1[0]*delr1[1]*mforce;
        virial[4] += delr1[0]*delr1[2]*mforce;
        virial[5] += delr1[1]*delr1[2]*mforce;
      }

      int nbor_k = nborj_start-offset_j+offset_k;
      for ( ; nbor_k<nbor_end; nbor_k+=n_stride) {
        int k=dev_nbor[nbor_k];
        k &= NEIGHMASK;

        if (j == k) continue;

        numtyp4 kx; fetch4(kx,k,pos_tex);
        int ktype=kx.w;
        ktype=map[ktype];
        int ijkparam=elem2param[itype*nelements*nelements+jtype*nelements+ktype];

        numtyp delr2[3];
        delr2[0] = kx.x-ix.x;
        delr2[1] = kx.y-ix.y;
        delr2[2] = kx.z-ix.z;
        numtyp rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];

        numtyp r2 = ucl_sqrt(rsq2);
        numtyp r2inv = ucl_rsqrt(rsq2);

        numtyp fi[3], fj[3], fk[3];
        numtyp4 ts1_ijkparam = ts1[ijkparam]; //fetch4(ts1_ijkparam,ijkparam,ts1_tex);
        lam3 = ts1_ijkparam.z;
        powermint = ts1_ijkparam.w;
        numtyp4 ts2_ijkparam = ts2[ijkparam]; //fetch4(ts2_ijkparam,ijkparam,ts2_tex);
        bigr = ts2_ijkparam.z;
        bigd = ts2_ijkparam.w;
        numtyp4 ts4_ijkparam = ts4[ijkparam]; //fetch4(ts4_ijkparam,ijkparam,ts4_tex);
        c1 = ts4_ijkparam.x;
        c2 = ts4_ijkparam.y;
        c3 = ts4_ijkparam.z;
        c4 = ts4_ijkparam.w;
        numtyp4 ts5_ijkparam = ts5[ijkparam]; //fetch4(ts5_ijkparam,ijkparam,ts5_tex);
        c5 = ts5_ijkparam.x;
        h = ts5_ijkparam.y;
        if (EVFLAG && vflag)
          attractive(bigr, bigd, powermint, lam3, h, c1, c2, c3, c4, c5,
                     prefactor, r1, r1inv, r2, r2inv, delr1, delr2, fi, fj, fk);
        else
          attractive_fi(bigr, bigd, powermint, lam3, h, c1, c2, c3, c4, c5,
                        prefactor, r1, r1inv, r2, r2inv, delr1, delr2, fi);
        f.x += fi[0];
        f.y += fi[1];
        f.z += fi[2];

        if (EVFLAG && vflag) {
          acctyp v[6];
          numtyp pre = (numtyp)2.0;
          if (evatom==1) pre = TWOTHIRD;
          v[0] = pre*(delr1[0]*fj[0] + delr2[0]*fk[0]);
          v[1] = pre*(delr1[1]*fj[1] + delr2[1]*fk[1]);
          v[2] = pre*(delr1[2]*fj[2] + delr2[2]*fk[2]);
          v[3] = pre*(delr1[0]*fj[1] + delr2[0]*fk[1]);
          v[4] = pre*(delr1[0]*fj[2] + delr2[0]*fk[2]);
          v[5] = pre*(delr1[1]*fj[2] + delr2[1]*fk[2]);

          virial[0] += v[0]; virial[1] += v[1]; virial[2] += v[2];
          virial[3] += v[3]; virial[4] += v[4]; virial[5] += v[5];
        }
      } // nbor_k
    } // for nbor_j
  } // if ii
  store_answers(f,energy,virial,ii,inum,tid,tpa_sq,
                offset,eflag,vflag,ans,engv);
}

__kernel void k_tersoff_mod_three_end(const __global numtyp4 *restrict x_,
                                  const __global numtyp4 *restrict ts1_in,
                                  const __global numtyp4 *restrict ts2_in,
                                  const __global numtyp4 *restrict ts4_in,
                                  const __global numtyp4 *restrict ts5_in,
                                  const __global int *restrict map,
                                  const __global int *restrict elem2param,
                                  const int nelements, const int nparams,
                                  const __global acctyp4 *restrict zetaij,
                                  const __global int * dev_nbor,
                                  const __global int * dev_ilist,
                                  __global acctyp3 *restrict ans,
                                  __global acctyp *restrict engv,
                                  const int eflag, const int vflag,
                                  const int inum,  const int nbor_pitch,
                                  const int t_per_atom, const int gpu_nbor) {
  const int tpa_sq=fast_mul(t_per_atom,t_per_atom);
  numtyp lam3, powermint, bigr, bigd, c1, c2, c3, c4, c5, h;

  int tid, ii, offset, n_stride;
  atom_info(tpa_sq,ii,tid,offset);

  local_allocate_store_three();

  __local numtyp4 ts1[SHARED_SIZE];
  __local numtyp4 ts2[SHARED_SIZE];
  __local numtyp4 ts4[SHARED_SIZE];
  __local numtyp4 ts5[SHARED_SIZE];
  if (tid<nparams) {
    ts1[tid]=ts1_in[tid];
    ts2[tid]=ts2_in[tid];
    ts4[tid]=ts4_in[tid];
    ts5[tid]=ts5_in[tid];
  }

  acctyp3 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp energy, virial[6];
  if (EVFLAG) {
    energy=(acctyp)0;
    for (int i=0; i<6; i++) virial[i]=(acctyp)0;
  }

  #ifdef LAL_SIMD_IP_SYNC
  __local int localk[BLOCK_PAIR];
  #endif

  __syncthreads();

  if (ii<inum) {
    int i, numj, nbor_j, nbor_end, k_end;
    int offset_j=offset/t_per_atom;
    nbor_info_p(dev_nbor,nbor_pitch,t_per_atom,ii,offset_j,i,numj,
                n_stride,nbor_end,nbor_j);
    int offset_k=tid & (t_per_atom-1);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int itype=ix.w;
    itype=map[itype];

    numtyp tpainv = ucl_recip((numtyp)t_per_atom);

    for ( ; nbor_j<nbor_end; nbor_j+=n_stride) {

      int j=dev_nbor[nbor_j];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int jtype=jx.w;
      jtype=map[jtype];

      // Compute r12
      numtyp delr1[3];
      delr1[0] = jx.x-ix.x;
      delr1[1] = jx.y-ix.y;
      delr1[2] = jx.z-ix.z;
      numtyp rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];

      numtyp mdelr1[3];
      mdelr1[0] = -delr1[0];
      mdelr1[1] = -delr1[1];
      mdelr1[2] = -delr1[2];

      int nbor_k;
      if (gpu_nbor) nbor_k=j+nbor_pitch;
      else nbor_k=dev_ilist[j]+nbor_pitch;
      const int numk=dev_nbor[nbor_k];
      nbor_k+=nbor_pitch+fast_mul(j,t_per_atom-1);
      k_end=nbor_k+fast_mul(numk/t_per_atom,n_stride)+(numk & (t_per_atom-1));
      nbor_k+=offset_k;

      int nbork_start = nbor_k;
      // look up for zeta_ji: find i in the j's neighbor list
      #ifdef LAL_SIMD_IP_SYNC
      int m = tid / t_per_atom;
      #endif
      int ijnum;
      for ( ; nbor_k<k_end; nbor_k+=n_stride) {
        int k=dev_nbor[nbor_k];
        k &= NEIGHMASK;
        if (k == i) {
          #ifdef LAL_SIMD_IP_SYNC
          localk[m] = nbor_k;
          #else
          ijnum = nbor_k;
          #endif
          break;
        }
      }

      numtyp r1 = ucl_sqrt(rsq1);
      numtyp r1inv = ucl_rsqrt(rsq1);
      #ifdef LAL_SIMD_IP_SYNC
      simdsync();
      ijnum = localk[m];
      #endif

      acctyp4 zeta_ji = zetaij[ijnum-2*nbor_pitch];
      numtyp force = zeta_ji.x*tpainv;
      numtyp prefactor_ji = zeta_ji.y;
      f.x += delr1[0]*force;
      f.y += delr1[1]*force;
      f.z += delr1[2]*force;

      if (EVFLAG && eflag) {
        energy+=zeta_ji.z*tpainv;
      }
      if (EVFLAG && vflag) {
        numtyp mforce = -force;
        virial[0] += mdelr1[0]*mdelr1[0]*mforce;
        virial[1] += mdelr1[1]*mdelr1[1]*mforce;
        virial[2] += mdelr1[2]*mdelr1[2]*mforce;
        virial[3] += mdelr1[0]*mdelr1[1]*mforce;
        virial[4] += mdelr1[0]*mdelr1[2]*mforce;
        virial[5] += mdelr1[1]*mdelr1[2]*mforce;
      }

      // attractive forces
      for (nbor_k = nbork_start ; nbor_k<k_end; nbor_k+=n_stride) {
        int k=dev_nbor[nbor_k];
        k &= NEIGHMASK;

        if (k == i) continue;

        numtyp4 kx; fetch4(kx,k,pos_tex);
        int ktype=kx.w;
        ktype=map[ktype];
        int jikparam=elem2param[jtype*nelements*nelements+itype*nelements+ktype];

        numtyp delr2[3];
        delr2[0] = kx.x-jx.x;
        delr2[1] = kx.y-jx.y;
        delr2[2] = kx.z-jx.z;
        numtyp rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];

        numtyp r2 = ucl_sqrt(rsq2);
        numtyp r2inv = ucl_rsqrt(rsq2);
        numtyp4 ts1_param, ts2_param, ts4_param, ts5_param;
        numtyp fi[3];

        ts1_param = ts1[jikparam]; //fetch4(ts1_jikparam,jikparam,ts1_tex);
        lam3 = ts1_param.z;
        powermint = ts1_param.w;
        ts2_param = ts2[jikparam]; //fetch4(ts2_jikparam,jikparam,ts2_tex);
        bigr = ts2_param.z;
        bigd = ts2_param.w;
        ts4_param = ts4[jikparam]; //fetch4(ts4_jikparam,jikparam,ts4_tex);
        c1 = ts4_param.x;
        c2 = ts4_param.y;
        c3 = ts4_param.z;
        c4 = ts4_param.w;
        ts5_param = ts5[jikparam]; //fetch4(ts5_jikparam,jikparam,ts5_tex);
        c5 = ts5_param.x;
        h = ts5_param.y;
        attractive_fj(bigr, bigd, powermint, lam3, h, c1, c2, c3, c4, c5,
                      prefactor_ji, r1, r1inv, r2, r2inv, mdelr1, delr2, fi);
        f.x += fi[0];
        f.y += fi[1];
        f.z += fi[2];

        acctyp4 zeta_jk = zetaij[nbor_k-2*nbor_pitch];
        numtyp prefactor_jk = zeta_jk.y;
        int jkiparam=elem2param[jtype*nelements*nelements+ktype*nelements+itype];
        ts1_param = ts1[jkiparam]; //fetch4(ts1_jkiparam,jkiparam,ts1_tex);
        lam3 = ts1_param.z;
        powermint = ts1_param.w;
        ts2_param = ts2[jkiparam]; //fetch4(ts2_jkiparam,jkiparam,ts2_tex);
        bigr = ts2_param.z;
        bigd = ts2_param.w;
        ts4_param = ts4[jkiparam]; //fetch4(ts4_jkiparam,jkiparam,ts4_tex);
        c1 = ts4_param.x;
        c2 = ts4_param.y;
        c3 = ts4_param.z;
        c4 = ts4_param.w;
        ts5_param = ts5[jkiparam]; //fetch4(ts5_ikiparam,jkiparam,ts5_tex);
        c5 = ts5_param.x;
        h = ts5_param.y;
        attractive_fk(bigr, bigd, powermint, lam3, h, c1, c2, c3, c4, c5,
                      prefactor_jk, r2, r2inv, r1, r1inv, delr2, mdelr1, fi);
        f.x += fi[0];
        f.y += fi[1];
        f.z += fi[2];
      } // for nbor_k
    } // for nbor_j
  } // if ii
  #ifdef THREE_CONCURRENT
  store_answers(f,energy,virial,ii,inum,tid,tpa_sq,offset,
                eflag,vflag,ans,engv);
  #else
  store_answers_p(f,energy,virial,ii,inum,tid,tpa_sq,offset,
                  eflag,vflag,ans,engv,NUM_BLOCKS_X);
  #endif
}

__kernel void k_tersoff_mod_three_end_vatom(const __global numtyp4 *restrict x_,
                                        const __global numtyp4 *restrict ts1_in,
                                        const __global numtyp4 *restrict ts2_in,
                                        const __global numtyp4 *restrict ts4_in,
                                        const __global numtyp4 *restrict ts5_in,
                                        const __global int *restrict map,
                                        const __global int *restrict elem2param,
                                        const int nelements, const int nparams,
                                        const __global acctyp4 *restrict zetaij,
                                        const __global int * dev_nbor,
                                        const __global int * dev_ilist,
                                        __global acctyp3 *restrict ans,
                                        __global acctyp *restrict engv,
                                        const int eflag, const int vflag,
                                        const int inum,  const int nbor_pitch,
                                        const int t_per_atom, const int gpu_nbor) {
  const int tpa_sq=fast_mul(t_per_atom,t_per_atom);
  numtyp lam3, powermint, bigr, bigd, c1, c2, c3, c4, c5, h;

  int tid, ii, offset, n_stride;
  atom_info(tpa_sq,ii,tid,offset);

  local_allocate_store_three();

  __local numtyp4 ts1[SHARED_SIZE];
  __local numtyp4 ts2[SHARED_SIZE];
  __local numtyp4 ts4[SHARED_SIZE];
  __local numtyp4 ts5[SHARED_SIZE];
  if (tid<nparams) {
    ts1[tid]=ts1_in[tid];
    ts2[tid]=ts2_in[tid];
    ts4[tid]=ts4_in[tid];
    ts5[tid]=ts5_in[tid];
  }

  acctyp3 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp energy, virial[6];
  if (EVFLAG) {
    energy=(acctyp)0;
    for (int i=0; i<6; i++) virial[i]=(acctyp)0;
  }

  #ifdef LAL_SIMD_IP_SYNC
  __local int localk[BLOCK_PAIR];
  #endif

  __syncthreads();

  if (ii<inum) {
    int i, numj, nbor_j, nbor_end, k_end;
    int offset_j=offset/t_per_atom;
    nbor_info_p(dev_nbor,nbor_pitch,t_per_atom,ii,offset_j,i,numj,
                n_stride,nbor_end,nbor_j);
    int offset_k=tid & (t_per_atom-1);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int itype=ix.w;
    itype=map[itype];

    numtyp tpainv = ucl_recip((numtyp)t_per_atom);

    for ( ; nbor_j<nbor_end; nbor_j+=n_stride) {

      int j=dev_nbor[nbor_j];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int jtype=jx.w;
      jtype=map[jtype];

      // Compute r12
      numtyp delr1[3];
      delr1[0] = jx.x-ix.x;
      delr1[1] = jx.y-ix.y;
      delr1[2] = jx.z-ix.z;
      numtyp rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];

      numtyp mdelr1[3];
      mdelr1[0] = -delr1[0];
      mdelr1[1] = -delr1[1];
      mdelr1[2] = -delr1[2];

      int nbor_k;
      if (gpu_nbor) nbor_k=j+nbor_pitch;
      else nbor_k=dev_ilist[j]+nbor_pitch;
      const int numk=dev_nbor[nbor_k];
      nbor_k+=nbor_pitch+fast_mul(j,t_per_atom-1);
      k_end=nbor_k+fast_mul(numk/t_per_atom,n_stride)+(numk & (t_per_atom-1));
      nbor_k+=offset_k;

      int nbork_start = nbor_k;

      // look up for zeta_ji
      #ifdef LAL_SIMD_IP_SYNC
      int m = tid / t_per_atom;
      #endif
      int ijnum;
      for ( ; nbor_k<k_end; nbor_k+=n_stride) {
        int k=dev_nbor[nbor_k];
        k &= NEIGHMASK;
        if (k == i) {
          #ifdef LAL_SIMD_IP_SYNC
          localk[m] = nbor_k;
          #else
          ijnum = nbor_k;
          #endif
          break;
        }
      }

      numtyp r1 = ucl_sqrt(rsq1);
      numtyp r1inv = ucl_rsqrt(rsq1);
      #ifdef LAL_SIMD_IP_SYNC
      simdsync();
      ijnum = localk[m];
      #endif

      acctyp4 zeta_ji = zetaij[ijnum-2*nbor_pitch];
      numtyp force = zeta_ji.x*tpainv;
      numtyp prefactor_ji = zeta_ji.y;
      f.x += delr1[0]*force;
      f.y += delr1[1]*force;
      f.z += delr1[2]*force;

      if (EVFLAG && eflag) {
        energy+=zeta_ji.z*tpainv;
      }
      if (EVFLAG && vflag) {
        numtyp mforce = -force;
        virial[0] += mdelr1[0]*mdelr1[0]*mforce;
        virial[1] += mdelr1[1]*mdelr1[1]*mforce;
        virial[2] += mdelr1[2]*mdelr1[2]*mforce;
        virial[3] += mdelr1[0]*mdelr1[1]*mforce;
        virial[4] += mdelr1[0]*mdelr1[2]*mforce;
        virial[5] += mdelr1[1]*mdelr1[2]*mforce;
      }

      // attractive forces
      for (nbor_k = nbork_start; nbor_k<k_end; nbor_k+=n_stride) {
        int k=dev_nbor[nbor_k];
        k &= NEIGHMASK;

        if (k == i) continue;

        numtyp4 kx; fetch4(kx,k,pos_tex);
        int ktype=kx.w;
        ktype=map[ktype];
        int jikparam=elem2param[jtype*nelements*nelements+itype*nelements+ktype];

        numtyp delr2[3];
        delr2[0] = kx.x-jx.x;
        delr2[1] = kx.y-jx.y;
        delr2[2] = kx.z-jx.z;
        numtyp rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];
        numtyp r2 = ucl_sqrt(rsq2);
        numtyp r2inv = ucl_rsqrt(rsq2);

        numtyp fi[3], fj[3], fk[3];
        numtyp4 ts1_param, ts2_param, ts4_param, ts5_param;
        ts1_param = ts1[jikparam]; //fetch4(ts1_jikparam,jikparam,ts1_tex);
        lam3 = ts1_param.z;
        powermint = ts1_param.w;
        ts2_param = ts2[jikparam]; //fetch4(ts2_jikparam,jikparam,ts2_tex);
        bigr = ts2_param.z;
        bigd = ts2_param.w;
        ts4_param = ts4[jikparam]; //fetch4(ts4_jikparam,jikparam,ts4_tex);
        c1 = ts4_param.x;
        c2 = ts4_param.y;
        c3 = ts4_param.z;
        c4 = ts4_param.w;
        ts5_param = ts5[jikparam]; //fetch4(ts5_jijparam,jikparam,ts5_tex);
        c5 = ts5_param.x;
        h = ts5_param.y;
        attractive(bigr, bigd, powermint, lam3, h, c1, c2, c3, c4, c5,
                   prefactor_ji, r1, r1inv, r2, r2inv, mdelr1, delr2, fi, fj, fk);
        f.x += fj[0];
        f.y += fj[1];
        f.z += fj[2];

        virial[0] += TWOTHIRD*(mdelr1[0]*fj[0] + delr2[0]*fk[0]);
        virial[1] += TWOTHIRD*(mdelr1[1]*fj[1] + delr2[1]*fk[1]);
        virial[2] += TWOTHIRD*(mdelr1[2]*fj[2] + delr2[2]*fk[2]);
        virial[3] += TWOTHIRD*(mdelr1[0]*fj[1] + delr2[0]*fk[1]);
        virial[4] += TWOTHIRD*(mdelr1[0]*fj[2] + delr2[0]*fk[2]);
        virial[5] += TWOTHIRD*(mdelr1[1]*fj[2] + delr2[1]*fk[2]);

        acctyp4 zeta_jk = zetaij[nbor_k-2*nbor_pitch];
        numtyp prefactor_jk = zeta_jk.y;

        int jkiparam=elem2param[jtype*nelements*nelements+ktype*nelements+itype];
        ts1_param = ts1[jkiparam]; //fetch4(ts1_jkiparam,jkiparam,ts1_tex);
        lam3 = ts1_param.z;
        powermint = ts1_param.w;
        ts2_param = ts2[jkiparam]; //fetch4(ts2_jkiparam,jkiparam,ts2_tex);
        bigr = ts2_param.z;
        bigd = ts2_param.w;
        ts4_param = ts4[jkiparam]; //fetch4(ts4_jkiparam,jkiparam,ts4_tex);
        c1 = ts4_param.x;
        c2 = ts4_param.y;
        c3 = ts4_param.z;
        c4 = ts4_param.w;
        ts5_param = ts5[jkiparam]; //fetch4(ts5_ikiparam,jkiparam,ts5_tex);
        c5 = ts5_param.x;
        h = ts5_param.y;
        attractive(bigr, bigd, powermint, lam3, h, c1, c2, c3, c4, c5,
                   prefactor_jk, r2, r2inv, r1, r1inv, delr2, mdelr1, fi, fj, fk);
        f.x += fk[0];
        f.y += fk[1];
        f.z += fk[2];

        virial[0] += TWOTHIRD*(delr2[0]*fj[0] + mdelr1[0]*fk[0]);
        virial[1] += TWOTHIRD*(delr2[1]*fj[1] + mdelr1[1]*fk[1]);
        virial[2] += TWOTHIRD*(delr2[2]*fj[2] + mdelr1[2]*fk[2]);
        virial[3] += TWOTHIRD*(delr2[0]*fj[1] + mdelr1[0]*fk[1]);
        virial[4] += TWOTHIRD*(delr2[0]*fj[2] + mdelr1[0]*fk[2]);
        virial[5] += TWOTHIRD*(delr2[1]*fj[2] + mdelr1[1]*fk[2]);
      }
    } // for nbor
  } // if ii
  #ifdef THREE_CONCURRENT
  store_answers(f,energy,virial,ii,inum,tid,tpa_sq,offset,
                eflag,vflag,ans,engv);
  #else
  store_answers_p(f,energy,virial,ii,inum,tid,tpa_sq,offset,
                  eflag,vflag,ans,engv,NUM_BLOCKS_X);
  #endif
}

