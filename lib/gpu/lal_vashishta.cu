// **************************************************************************
//                                 vashishta.cu
//                             -------------------
//                           Anders Hafreager (UiO)
//
//  Device code for acceleration of the vashishta pair style
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                : Mon June 12, 2017
//    email                : andershaf@gmail.com
// ***************************************************************************

#if defined(NV_KERNEL) || defined(USE_HIP)
#include "lal_aux_fun1.h"

#ifndef _DOUBLE_DOUBLE
_texture( pos_tex,float4);
_texture( param1_tex,float4);
_texture( param2_tex,float4);
_texture( param3_tex,float4);
_texture( param4_tex,float4);
_texture( param5_tex,float4);
#else
_texture_2d( pos_tex,int4);
_texture( param1_tex,int4);
_texture( param2_tex,int4);
_texture( param3_tex,int4);
_texture( param4_tex,int4);
_texture( param5_tex,int4);
#endif

#if (__CUDACC_VER_MAJOR__ >= 11)
#define param1_tex param1
#define param2_tex param2
#define param3_tex param3
#define param4_tex param4
#define param5_tex param5
#endif

#else
#define pos_tex x_
#define param1_tex param1
#define param2_tex param2
#define param3_tex param3
#define param4_tex param4
#define param5_tex param5
#endif



#define THIRD (numtyp)0.66666666666666666667

//#define THREE_CONCURRENT

#if (SHUFFLE_AVAIL == 0)

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

__kernel void k_vashishta_short_nbor(const __global numtyp4 *restrict x_,
                                     const __global numtyp4 *restrict param4,
                                     const __global int *restrict map,
                                     const __global int *restrict elem2param,
                                     const int nelements, const int nparams,
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
    int itype=ix.w;
    itype=map[itype];
    int newj=0;

    __global int *out_list=dev_nbor+2*nbor_pitch+ii*t_per_atom;
    const int out_stride=nbor_pitch*t_per_atom-t_per_atom;

    for ( ; nbor<nbor_end; nbor+=nbor_pitch) {
      ucl_prefetch(dev_packed+nbor+nbor_pitch);
      int sj=dev_packed[nbor];
      int j = sj & NEIGHMASK;
      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];

      int jtype=jx.w;
      jtype=map[jtype];
      int ijparam=elem2param[itype*nelements*nelements+jtype*nelements+jtype];

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      if (rsq<param4[ijparam].x) {
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

__kernel void k_vashishta(const __global numtyp4 *restrict x_,
                   const __global numtyp4 *restrict param1,
                   const __global numtyp4 *restrict param2,
                   const __global numtyp4 *restrict param3,
                   const __global numtyp4 *restrict param4,
                   const __global numtyp4 *restrict param5,
                   const __global int *restrict map,
                   const __global int *restrict elem2param,
                   const int nelements,
                   const __global int * dev_packed,
                   __global acctyp3 *restrict ans,
                   __global acctyp *restrict engv,
                   const int eflag, const int vflag, const int inum,
                   const int nbor_pitch, const int ev_stride) {
  const int ii=GLOBAL_ID_X;

  local_allocate_store_pair();

  acctyp3 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp energy, virial[6];
  if (EVFLAG) {
    energy=(acctyp)0;
    for (int i=0; i<6; i++) virial[i]=(acctyp)0;
  }

  __syncthreads();

  if (ii<inum) {
    const int i=dev_packed[ii];
    int nbor=ii+nbor_pitch;
    const int numj=dev_packed[nbor];
    nbor+=nbor_pitch;
    const int nbor_end=nbor+fast_mul(numj,nbor_pitch);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    int itype=ix.w;
    itype=map[itype];

    for ( ; nbor<nbor_end; nbor+=nbor_pitch) {
      ucl_prefetch(dev_packed+nbor+nbor_pitch);

      int j=dev_packed[nbor];
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

      if (rsq<param4[ijparam].z) { // cutsq = param4[ijparam].z
        numtyp4 param1_ijparam; fetch4(param1_ijparam,ijparam,param1_tex);
        numtyp param1_eta=param1_ijparam.x;
        numtyp param1_lam1inv=param1_ijparam.y;
        numtyp param1_lam4inv=param1_ijparam.z;
        numtyp param1_zizj=param1_ijparam.w;

        numtyp4 param2_ijparam; fetch4(param2_ijparam,ijparam,param2_tex);
        numtyp param2_mbigd=param2_ijparam.x;
        numtyp param2_dvrc =param2_ijparam.y;
        numtyp param2_big6w=param2_ijparam.z;
        numtyp param2_heta =param2_ijparam.w;

        numtyp4 param3_ijparam; fetch4(param3_ijparam,ijparam,param3_tex);
        numtyp param3_bigh=param3_ijparam.x;
        numtyp param3_bigw=param3_ijparam.y;
        numtyp param3_dvrc=param3_ijparam.z;
        numtyp param3_c0  =param3_ijparam.w;

        numtyp r=ucl_sqrt(rsq);
        numtyp rinvsq=1.0/rsq;
        numtyp r4inv = rinvsq*rinvsq;
        numtyp r6inv = rinvsq*r4inv;

        numtyp reta = ucl_powr(r,-param1_eta);
        numtyp lam1r = r*param1_lam1inv;
        numtyp lam4r = r*param1_lam4inv;
        numtyp vc2 = param1_zizj * ucl_exp(-lam1r)/r;
        numtyp vc3 = param2_mbigd * r4inv*ucl_exp(-lam4r);

        numtyp force = (param2_dvrc*r
            - (4.0*vc3 + lam4r*vc3+param2_big6w*r6inv
               - param2_heta*reta - vc2 - lam1r*vc2)
            ) * rinvsq;

        f.x+=delx*force;
        f.y+=dely*force;
        f.z+=delz*force;

        if (EVFLAG && eflag)
          energy += (param3_bigh*reta+vc2-vc3-param3_bigw*r6inv-r*param3_dvrc+param3_c0);

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
  const int tid=THREAD_ID_X;
  store_answers_p(f,energy,virial,ii,inum,tid,1,0,eflag,vflag,ans,engv,
                  ev_stride);
}

#define threebody(delr1x, delr1y, delr1z, eflag, energy)                     \
{                                                                            \
  numtyp r1 = ucl_sqrt(rsq1);                                                \
  numtyp rinvsq1 = ucl_recip(rsq1);                                          \
  numtyp rainv1 = ucl_recip(r1 - param_r0_ij);                               \
  numtyp gsrainv1 = param_gamma_ij * rainv1;                                 \
  numtyp gsrainvsq1 = gsrainv1*rainv1/r1;                                    \
  numtyp expgsrainv1 = ucl_exp(gsrainv1);                                    \
                                                                             \
  numtyp r2 = ucl_sqrt(rsq2);                                                \
  numtyp rinvsq2 = ucl_recip(rsq2);                                          \
  numtyp rainv2 = ucl_recip(r2 - param_r0_ik);                               \
  numtyp gsrainv2 = param_gamma_ik * rainv2;                                 \
  numtyp gsrainvsq2 = gsrainv2*rainv2/r2;                                    \
  numtyp expgsrainv2 = ucl_exp(gsrainv2);                                    \
                                                                             \
  numtyp rinv12 = ucl_recip(r1*r2);                                          \
  numtyp cs = (delr1x*delr2x + delr1y*delr2y + delr1z*delr2z) * rinv12;      \
  numtyp delcs = cs - param_costheta_ijk;                                    \
  numtyp delcssq = delcs*delcs;                                              \
  numtyp pcsinv = param_bigc_ijk*delcssq+1.0;                                \
  numtyp pcsinvsq = pcsinv*pcsinv;                                           \
  numtyp pcs = delcssq/pcsinv;                                               \
                                                                             \
  numtyp facexp = expgsrainv1*expgsrainv2;                                   \
                                                                             \
  numtyp facrad = param_bigb_ijk * facexp*pcs;                               \
  numtyp frad1 = facrad*gsrainvsq1;                                          \
  numtyp frad2 = facrad*gsrainvsq2;                                          \
  numtyp facang = param_big2b_ijk * facexp*delcs/pcsinvsq;                   \
  numtyp facang12 = rinv12*facang;                                           \
  numtyp csfacang = cs*facang;                                               \
  numtyp csfac1 = rinvsq1*csfacang;                                          \
                                                                             \
  fjx = delr1x*(frad1+csfac1)-delr2x*facang12;                               \
  fjy = delr1y*(frad1+csfac1)-delr2y*facang12;                               \
  fjz = delr1z*(frad1+csfac1)-delr2z*facang12;                               \
                                                                             \
  numtyp csfac2 = rinvsq2*csfacang;                                          \
                                                                             \
  fkx = delr2x*(frad2+csfac2)-delr1x*facang12;                               \
  fky = delr2y*(frad2+csfac2)-delr1y*facang12;                               \
  fkz = delr2z*(frad2+csfac2)-delr1z*facang12;                               \
                                                                             \
  if (EVFLAG && eflag)                                                       \
    energy+=facrad;                                                          \
  if (EVFLAG && vflag) {                                                     \
    virial[0] += delr1x*fjx + delr2x*fkx;                                    \
    virial[1] += delr1y*fjy + delr2y*fky;                                    \
    virial[2] += delr1z*fjz + delr2z*fkz;                                    \
    virial[3] += delr1x*fjy + delr2x*fky;                                    \
    virial[4] += delr1x*fjz + delr2x*fkz;                                    \
    virial[5] += delr1y*fjz + delr2y*fkz;                                    \
  }                                                                          \
}

#define threebody_half(delr1x, delr1y, delr1z)                               \
{                                                                            \
  numtyp r1 = ucl_sqrt(rsq1);                                                \
  numtyp rinvsq1 = ucl_recip(rsq1);                                          \
  numtyp rainv1 = ucl_recip(r1 - param_r0_ij);                               \
  numtyp gsrainv1 = param_gamma_ij * rainv1;                                 \
  numtyp gsrainvsq1 = gsrainv1*rainv1/r1;                                    \
  numtyp expgsrainv1 = ucl_exp(gsrainv1);                                    \
                                                                             \
  numtyp r2 = ucl_sqrt(rsq2);                                                \
  numtyp rainv2 = ucl_recip(r2 - param_r0_ik);                               \
  numtyp gsrainv2 = param_gamma_ik * rainv2;                                 \
  numtyp expgsrainv2 = ucl_exp(gsrainv2);                                    \
                                                                             \
  numtyp rinv12 = ucl_recip(r1*r2);                                          \
  numtyp cs = (delr1x*delr2x + delr1y*delr2y + delr1z*delr2z) * rinv12;      \
  numtyp delcs = cs - param_costheta_ijk;                                    \
  numtyp delcssq = delcs*delcs;                                              \
  numtyp pcsinv = param_bigc_ijk*delcssq+1.0;                                \
  numtyp pcsinvsq = pcsinv*pcsinv;                                           \
  numtyp pcs = delcssq/pcsinv;                                               \
                                                                             \
  numtyp facexp = expgsrainv1*expgsrainv2;                                   \
                                                                             \
  numtyp facrad = param_bigb_ijk * facexp*pcs;                               \
  numtyp frad1 = facrad*gsrainvsq1;                                          \
  numtyp facang = param_big2b_ijk * facexp*delcs/pcsinvsq;                   \
  numtyp facang12 = rinv12*facang;                                           \
  numtyp csfacang = cs*facang;                                               \
  numtyp csfac1 = rinvsq1*csfacang;                                          \
                                                                             \
  fjx = delr1x*(frad1+csfac1)-delr2x*facang12;                               \
  fjy = delr1y*(frad1+csfac1)-delr2y*facang12;                               \
  fjz = delr1z*(frad1+csfac1)-delr2z*facang12;                               \
}

__kernel void k_vashishta_three_center(const __global numtyp4 *restrict x_,
                                const __global numtyp4 *restrict param1,
                                const __global numtyp4 *restrict param2,
                                const __global numtyp4 *restrict param3,
                                const __global numtyp4 *restrict param4,
                                const __global numtyp4 *restrict param5,
                                const __global int *restrict map,
                                const __global int *restrict elem2param,
                                const int nelements,
                                const __global int * dev_nbor,
                                __global acctyp3 *restrict ans,
                                __global acctyp *restrict engv,
                                const int eflag, const int vflag,
                                const int inum,  const int nbor_pitch,
                                const int t_per_atom, const int evatom) {
  int n_stride;
  const int tpa_sq=fast_mul(t_per_atom,t_per_atom);
  numtyp param_gamma_ij, param_r0sq_ij, param_r0_ij, param_gamma_ik, param_r0sq_ik, param_r0_ik;
  numtyp param_costheta_ijk, param_bigc_ijk, param_bigb_ijk, param_big2b_ijk;

  int tid, ii, offset;
  atom_info(tpa_sq,ii,tid,offset);

  local_allocate_store_three();

  acctyp3 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp energy, virial[6];
  if (EVFLAG) {
    energy=(acctyp)0;
    for (int i=0; i<6; i++) virial[i]=(acctyp)0;
  }

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

    for ( ; nbor_j<nbor_end; nbor_j+=n_stride) {

      int j=dev_nbor[nbor_j];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int jtype=jx.w;
      jtype=map[jtype];

      // Compute r12
      numtyp delr1x = jx.x-ix.x;
      numtyp delr1y = jx.y-ix.y;
      numtyp delr1z = jx.z-ix.z;
      numtyp rsq1 = delr1x*delr1x+delr1y*delr1y+delr1z*delr1z;

      int ijparam=elem2param[itype*nelements*nelements+jtype*nelements+jtype];

      numtyp4 param4_ijparam; fetch4(param4_ijparam,ijparam,param4_tex);
      param_r0sq_ij=param4_ijparam.x;
      if (rsq1 > param_r0sq_ij) continue; // still keep this for neigh no and tpa > 1
      param_gamma_ij=param4_ijparam.y;
      param_r0_ij=param4_ijparam.w;

      int nbor_k = nbor_j-offset_j+offset_k;
      if (nbor_k<=nbor_j) nbor_k += n_stride;

      for ( ; nbor_k<nbor_end; nbor_k+=n_stride) {
        int k=dev_nbor[nbor_k];
        k &= NEIGHMASK;

        numtyp4 kx; fetch4(kx,k,pos_tex);
        int ktype=kx.w;
        ktype=map[ktype];
        int ikparam=elem2param[itype*nelements*nelements+ktype*nelements+ktype];
        numtyp4 param4_ikparam; fetch4(param4_ikparam,ikparam,param4_tex);

        numtyp delr2x = kx.x-ix.x;
        numtyp delr2y = kx.y-ix.y;
        numtyp delr2z = kx.z-ix.z;
        numtyp rsq2 = delr2x*delr2x + delr2y*delr2y + delr2z*delr2z;

        param_r0sq_ik=param4_ikparam.x;
        if (rsq2 < param_r0sq_ik) {
          param_gamma_ik=param4_ikparam.y;
          param_r0_ik=param4_ikparam.w;

          int ijkparam=elem2param[itype*nelements*nelements+jtype*nelements+ktype];
          numtyp4 param5_ijkparam; fetch4(param5_ijkparam,ijkparam,param5_tex);
          param_bigc_ijk=param5_ijkparam.x;
          param_bigb_ijk=param5_ijkparam.z;
          param_big2b_ijk=param5_ijkparam.w;
          param_costheta_ijk=param5_ijkparam.y;

          numtyp fjx, fjy, fjz, fkx, fky, fkz;
          threebody(delr1x,delr1y,delr1z,eflag,energy);

          f.x -= fjx + fkx;
          f.y -= fjy + fky;
          f.z -= fjz + fkz;
        }
      }
    } // for nbor

    numtyp pre;
    if (evatom==1)
      pre=THIRD;
    else
      pre=(numtyp)2.0;
    energy*=pre;
    for (int i=0; i<6; i++)
      virial[i]*=pre;
  } // if ii
  store_answers(f,energy,virial,ii,inum,tid,tpa_sq,offset,
                eflag,vflag,ans,engv);
}

__kernel void k_vashishta_three_end(const __global numtyp4 *restrict x_,
                             const __global numtyp4 *restrict param1,
                             const __global numtyp4 *restrict param2,
                             const __global numtyp4 *restrict param3,
                             const __global numtyp4 *restrict param4,
                             const __global numtyp4 *restrict param5,
                             const __global int *restrict map,
                             const __global int *restrict elem2param,
                             const int nelements,
                             const __global int * dev_nbor,
                             const __global int * dev_ilist,
                             __global acctyp3 *restrict ans,
                             __global acctyp *restrict engv,
                             const int eflag, const int vflag,
                             const int inum,  const int nbor_pitch,
                             const int t_per_atom, const int gpu_nbor) {
  int n_stride;
  const int tpa_sq=fast_mul(t_per_atom,t_per_atom);
  numtyp param_gamma_ij, param_r0sq_ij, param_r0_ij, param_gamma_ik, param_r0sq_ik, param_r0_ik;
  numtyp param_costheta_ijk, param_bigc_ijk, param_bigb_ijk, param_big2b_ijk;

  int tid, ii, offset;
  atom_info(tpa_sq,ii,tid,offset);

  local_allocate_store_three();

  acctyp3 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp energy, virial[6];
  if (EVFLAG) {
    energy=(acctyp)0;
    for (int i=0; i<6; i++) virial[i]=(acctyp)0;
  }

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

    for ( ; nbor_j<nbor_end; nbor_j+=n_stride) {
      int j=dev_nbor[nbor_j];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int jtype=jx.w;
      jtype=map[jtype];

      // Compute r12
      numtyp delr1x = ix.x-jx.x;
      numtyp delr1y = ix.y-jx.y;
      numtyp delr1z = ix.z-jx.z;
      numtyp rsq1 = delr1x*delr1x+delr1y*delr1y+delr1z*delr1z;

      int ijparam=elem2param[itype*nelements*nelements+jtype*nelements+jtype];
      numtyp4 param4_ijparam; fetch4(param4_ijparam,ijparam,param4_tex);
      param_r0sq_ij = param4_ijparam.x;
      if (rsq1 > param_r0sq_ij) continue; // still keep this for neigh no and tpa > 1

      param_gamma_ij=param4_ijparam.y;
      param_r0_ij = param4_ijparam.w;

      int nbor_k;
      if (gpu_nbor) nbor_k=j+nbor_pitch;
      else nbor_k=dev_ilist[j]+nbor_pitch;
      const int numk=dev_nbor[nbor_k];
      nbor_k+=nbor_pitch+fast_mul(j,t_per_atom-1);
      k_end=nbor_k+fast_mul(numk/t_per_atom,n_stride)+(numk & (t_per_atom-1));
      nbor_k+=offset_k;

      for ( ; nbor_k<k_end; nbor_k+=n_stride) {
        int k=dev_nbor[nbor_k];
        k &= NEIGHMASK;

        if (k == i) continue;

        numtyp4 kx; fetch4(kx,k,pos_tex);
        int ktype=kx.w;
        ktype=map[ktype];
        int ikparam=elem2param[jtype*nelements*nelements+ktype*nelements+ktype]; //jk

        numtyp delr2x = kx.x - jx.x;
        numtyp delr2y = kx.y - jx.y;
        numtyp delr2z = kx.z - jx.z;
        numtyp rsq2 = delr2x*delr2x + delr2y*delr2y + delr2z*delr2z;
        numtyp4 param4_ikparam; fetch4(param4_ikparam,ikparam,param4_tex);
        param_r0sq_ik=param4_ikparam.x;

        if (rsq2 < param_r0sq_ik) {
          param_gamma_ik=param4_ikparam.y;
          param_r0_ik=param4_ikparam.w;

          int ijkparam=elem2param[jtype*nelements*nelements+itype*nelements+ktype]; //jik
          numtyp4 param5_ijkparam; fetch4(param5_ijkparam,ijkparam,param5_tex);
          param_bigc_ijk=param5_ijkparam.x;
          param_costheta_ijk=param5_ijkparam.y;
          param_bigb_ijk=param5_ijkparam.z;
          param_big2b_ijk=param5_ijkparam.w;

          numtyp fjx, fjy, fjz;
          //if (evatom==0) {
            threebody_half(delr1x,delr1y,delr1z);
          //} else {
          //  numtyp fkx, fky, fkz;
          //  threebody(delr1x,delr1y,delr1z,eflag,energy);
          //}

          f.x += fjx;
          f.y += fjy;
          f.z += fjz;
        }
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

__kernel void k_vashishta_three_end_vatom(const __global numtyp4 *restrict x_,
                             const __global numtyp4 *restrict param1,
                             const __global numtyp4 *restrict param2,
                             const __global numtyp4 *restrict param3,
                             const __global numtyp4 *restrict param4,
                             const __global numtyp4 *restrict param5,
                             const __global int *restrict map,
                             const __global int *restrict elem2param,
                             const int nelements,
                             const __global int * dev_nbor,
                             const __global int * dev_ilist,
                             __global acctyp3 *restrict ans,
                             __global acctyp *restrict engv,
                             const int eflag, const int vflag,
                             const int inum,  const int nbor_pitch,
                             const int t_per_atom, const int gpu_nbor) {
  int n_stride;
  const int tpa_sq=fast_mul(t_per_atom,t_per_atom);
  numtyp param_gamma_ij, param_r0sq_ij, param_r0_ij, param_gamma_ik, param_r0sq_ik, param_r0_ik;
  numtyp param_costheta_ijk, param_bigc_ijk, param_bigb_ijk, param_big2b_ijk;

  int tid, ii, offset;
  atom_info(tpa_sq,ii,tid,offset);

  local_allocate_store_three();

  acctyp3 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp energy, virial[6];
  if (EVFLAG) {
    energy=(acctyp)0;
    for (int i=0; i<6; i++) virial[i]=(acctyp)0;
  }

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

    for ( ; nbor_j<nbor_end; nbor_j+=n_stride) {
      int j=dev_nbor[nbor_j];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      int jtype=jx.w;
      jtype=map[jtype];

      // Compute r12
      numtyp delr1x = ix.x-jx.x;
      numtyp delr1y = ix.y-jx.y;
      numtyp delr1z = ix.z-jx.z;
      numtyp rsq1 = delr1x*delr1x+delr1y*delr1y+delr1z*delr1z;

      int ijparam=elem2param[itype*nelements*nelements+jtype*nelements+jtype];
      numtyp4 param4_ijparam; fetch4(param4_ijparam,ijparam,param4_tex);
      param_r0sq_ij=param4_ijparam.x;
      if (rsq1 > param_r0sq_ij) continue;  // still keep this for neigh no and tpa > 1

      param_gamma_ij=param4_ijparam.y;
      param_r0_ij=param4_ijparam.w;

      int nbor_k;
      if (gpu_nbor) nbor_k=j+nbor_pitch;
      else nbor_k=dev_ilist[j]+nbor_pitch;
      const int numk=dev_nbor[nbor_k];
      nbor_k+=nbor_pitch+fast_mul(j,t_per_atom-1);
      k_end=nbor_k+fast_mul(numk/t_per_atom,n_stride)+(numk & (t_per_atom-1));
      nbor_k+=offset_k;

      for ( ; nbor_k<k_end; nbor_k+=n_stride) {
        int k=dev_nbor[nbor_k];
        k &= NEIGHMASK;

        if (k == i) continue;

        numtyp4 kx; fetch4(kx,k,pos_tex);
        int ktype=kx.w;
        ktype=map[ktype];
        int ikparam=elem2param[jtype*nelements*nelements+ktype*nelements+ktype]; // jk
        numtyp4 param4_ikparam; fetch4(param4_ikparam,ikparam,param4_tex);

        numtyp delr2x = kx.x - jx.x;
        numtyp delr2y = kx.y - jx.y;
        numtyp delr2z = kx.z - jx.z;
        numtyp rsq2 = delr2x*delr2x + delr2y*delr2y + delr2z*delr2z;
        param_r0sq_ik=param4_ikparam.x;

        if (rsq2 < param_r0sq_ik) {
          param_gamma_ik=param4_ikparam.y;
          param_r0_ik=param4_ikparam.w;

          int ijkparam=elem2param[jtype*nelements*nelements+itype*nelements+ktype]; // jik
          numtyp4 param5_ijkparam; fetch4(param5_ijkparam,ijkparam,param5_tex);
          param_bigc_ijk=param5_ijkparam.x;
          param_costheta_ijk=param5_ijkparam.y;
          param_bigb_ijk=param5_ijkparam.z;
          param_big2b_ijk=param5_ijkparam.w;

          numtyp fjx, fjy, fjz, fkx, fky, fkz;
          threebody(delr1x,delr1y,delr1z,eflag,energy);

          f.x += fjx;
          f.y += fjy;
          f.z += fjz;
        }
      }

    } // for nbor
    energy*=THIRD;
    for (int i=0; i<6; i++)
      virial[i]*=THIRD;
  } // if ii
  #ifdef THREE_CONCURRENT
  store_answers(f,energy,virial,ii,inum,tid,tpa_sq,offset,
                eflag,vflag,ans,engv);
  #else
  store_answers_p(f,energy,virial,ii,inum,tid,tpa_sq,offset,
                  eflag,vflag,ans,engv,NUM_BLOCKS_X);
  #endif
}

