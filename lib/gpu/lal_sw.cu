// **************************************************************************
//                                   sw.cu
//                             -------------------
//                           W. Michael Brown (ORNL)
//
//  Device code for acceleration of the sw pair style
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                : Tue March 26, 2013
//    email                : brownw@ornl.gov
// ***************************************************************************

#if defined(NV_KERNEL) || defined(USE_HIP)
#include "lal_aux_fun1.h"

#ifndef _DOUBLE_DOUBLE
_texture( pos_tex,float4);
_texture( sw1_tex,float4);
_texture( sw2_tex,float4);
_texture( sw3_tex,float4);
#else
_texture_2d( pos_tex,int4);
_texture( sw1_tex,int4);
_texture( sw2_tex,int4);
_texture( sw3_tex,int4);
#endif

#else
#define pos_tex x_
#define sw1_tex sw1
#define sw2_tex sw2
#define sw3_tex sw3
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
    acctyp4 old=ans[ii];                                                    \
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
    acctyp4 old=ans[ii];                                                    \
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
    acctyp4 old=ans[ii];                                                    \
    old.x+=f.x;                                                             \
    old.y+=f.y;                                                             \
    old.z+=f.z;                                                             \
    ans[ii]=old;                                                            \
  }

#endif
#endif

__kernel void k_sw_short_nbor(const __global numtyp4 *restrict x_,
                              const __global numtyp * restrict cutsq,
                              const int ntypes, __global int * dev_nbor,
                              const __global int * dev_packed,
                              const int inum, const int nbor_pitch,
                              const int t_per_atom) {
  const int ii=GLOBAL_ID_X;

  #ifdef ONETYPE
  const numtyp sw_cutsq=cutsq[ONETYPE];
  #endif

  if (ii<inum) {
    const int i=dev_packed[ii];
    int nbor=ii+nbor_pitch;
    const int numj=dev_packed[nbor];
    nbor+=nbor_pitch;
    const int nbor_end=nbor+fast_mul(numj,nbor_pitch);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    #ifndef ONETYPE
    const int itype=ix.w*ntypes;
    #endif
    int newj=0;

    __global int *out_list=dev_nbor+2*nbor_pitch+ii*t_per_atom;
    const int out_stride=nbor_pitch*t_per_atom-t_per_atom;

    for ( ; nbor<nbor_end; nbor+=nbor_pitch) {
      int sj=dev_packed[nbor];
      int j = sj & NEIGHMASK;
      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];

      #ifndef ONETYPE
      const int mtype=jx.w+itype;
      const numtyp sw_cutsq=cutsq[mtype];
      #endif

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      if (rsq<sw_cutsq) {
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

__kernel void k_sw(const __global numtyp4 *restrict x_,
                   const __global numtyp4 * restrict sw_pre,
                   const __global numtyp4 * restrict c_14,
                   const __global numtyp2 * restrict c_56,
                   const int ntypes, const __global int * dev_nbor,
                   __global acctyp4 *restrict ans,
                   __global acctyp *restrict engv,
                   const int eflag, const int vflag, const int inum,
                   const int nbor_pitch, const int t_per_atom,
                   const int ev_stride) {
  int n_stride;
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);
  local_allocate_store_pair();

  #ifdef ONETYPE
  const numtyp4 pre_sw=sw_pre[ONETYPE];
  const numtyp4 pre_sw_c14=c_14[ONETYPE];
  numtyp2 pre_sw_c56;
  if (EVFLAG && eflag) pre_sw_c56=c_56[ONETYPE];
  #endif

  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp energy, virial[6];
  if (EVFLAG) {
    energy=(acctyp)0;
    for (int i=0; i<6; i++) virial[i]=(acctyp)0;
  }

  if (ii<inum) {
    int nbor, nbor_end, i, numj;
    nbor_info_p(dev_nbor,nbor_pitch,t_per_atom,ii,offset,i,numj,
                n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    #ifndef ONETYPE
    int itype=ix.w;
    itype*=ntypes;
    #endif

    for ( ; nbor<nbor_end; nbor+=n_stride) {

      int j=dev_nbor[nbor];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      #ifndef ONETYPE
      int mtype=jx.w;
      mtype+=itype;
      #endif

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      #ifndef ONETYPE
      numtyp4 pre_sw=sw_pre[mtype];
      numtyp4 pre_sw_c14=c_14[mtype];
      #endif
      numtyp r=ucl_sqrt(rsq);
      #ifdef SPQ
      numtyp rp=r*r;
      rp=ucl_recip(rp*rp);
      numtyp rq=(numtyp)1.0;
      #else
      numtyp rp=ucl_powr(r,-pre_sw.z);
      numtyp rq=ucl_powr(r,-pre_sw.w);
      #endif
      numtyp rainv=ucl_recip(r-pre_sw.x);
      numtyp expsrainv=ucl_exp(pre_sw.y*rainv);
      rainv*=rainv*r;
      numtyp force = (pre_sw_c14.x*rp-pre_sw_c14.y*rq +
                     (pre_sw_c14.z*rp-pre_sw_c14.w*rq) * rainv)*
                     expsrainv*ucl_recip(rsq);

      f.x+=delx*force;
      f.y+=dely*force;
      f.z+=delz*force;

      if (EVFLAG && eflag) {
        #ifndef ONETYPE
        numtyp2 pre_sw_c56=c_56[mtype];
        #endif
        energy+=(pre_sw_c56.x*rp - pre_sw_c56.y*rq) * expsrainv;
      }

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

#define threebody(delr1x,delr1y,delr1z,delr2x,delr2y,delr2z, eflag, energy)  \
{                                                                            \
  numtyp r1 = ucl_sqrt(rsq1);                                                \
  numtyp rinvsq1 = ucl_recip(rsq1);                                          \
  numtyp rainv1 = ucl_recip(r1 - sw_cut_ij);                                 \
  numtyp gsrainv1 = sw_sigma_gamma_ij * rainv1;                              \
  numtyp gsrainvsq1 = gsrainv1*rainv1/r1;                                    \
  numtyp expgsrainv1 = ucl_exp(gsrainv1);                                    \
                                                                             \
  numtyp r2 = ucl_sqrt(rsq2);                                                \
  numtyp rinvsq2 = ucl_recip(rsq2);                                          \
  numtyp rainv2 = ucl_recip(r2 - sw_cut_ik);                                 \
  numtyp gsrainv2 = sw_sigma_gamma_ik * rainv2;                              \
  numtyp gsrainvsq2 = gsrainv2*rainv2/r2;                                    \
  numtyp expgsrainv2 = ucl_exp(gsrainv2);                                    \
                                                                             \
  numtyp rinv12 = ucl_recip(r1*r2);                                          \
  numtyp cs = (delr1x*delr2x + delr1y*delr2y + delr1z*delr2z) * rinv12;      \
  numtyp delcs = cs - sw_costheta_ijk;                                       \
  numtyp delcssq = delcs*delcs;                                              \
                                                                             \
  numtyp facexp = expgsrainv1*expgsrainv2;                                   \
                                                                             \
  numtyp facrad = sw_lambda_epsilon_ijk * facexp*delcssq;                    \
  numtyp frad1 = facrad*gsrainvsq1;                                          \
  numtyp frad2 = facrad*gsrainvsq2;                                          \
  numtyp facang = (numtyp)2.0 * sw_lambda_epsilon_ijk * facexp*delcs;        \
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

#define threebody_half(delr1x, delr1y, delr1z, delr2x, delr2y, delr2z)       \
{                                                                            \
  numtyp r1 = ucl_sqrt(rsq1);                                                \
  numtyp rinvsq1 = ucl_recip(rsq1);                                          \
  numtyp rainv1 = ucl_recip(r1 - sw_cut_ij);                                 \
  numtyp gsrainv1 = sw_sigma_gamma_ij * rainv1;                              \
  numtyp gsrainvsq1 = gsrainv1*rainv1/r1;                                    \
  numtyp expgsrainv1 = ucl_exp(gsrainv1);                                    \
                                                                             \
  numtyp r2 = ucl_sqrt(rsq2);                                                \
  numtyp rainv2 = ucl_recip(r2 - sw_cut_ik);                                 \
  numtyp gsrainv2 = sw_sigma_gamma_ik * rainv2;                              \
  numtyp expgsrainv2 = ucl_exp(gsrainv2);                                    \
                                                                             \
  numtyp rinv12 = ucl_recip(r1*r2);                                          \
  numtyp cs = (delr1x*delr2x + delr1y*delr2y + delr1z*delr2z) * rinv12;      \
  numtyp delcs = cs - sw_costheta_ijk;                                       \
  numtyp delcssq = delcs*delcs;                                              \
                                                                             \
  numtyp facexp = expgsrainv1*expgsrainv2;                                   \
                                                                             \
  numtyp facrad = sw_lambda_epsilon_ijk * facexp*delcssq;                    \
  numtyp frad1 = facrad*gsrainvsq1;                                          \
  numtyp facang = (numtyp)2.0 * sw_lambda_epsilon_ijk * facexp*delcs;        \
  numtyp facang12 = rinv12*facang;                                           \
  numtyp csfacang = cs*facang;                                               \
  numtyp csfac1 = rinvsq1*csfacang;                                          \
                                                                             \
  fjx = delr1x*(frad1+csfac1)-delr2x*facang12;                               \
  fjy = delr1y*(frad1+csfac1)-delr2y*facang12;                               \
  fjz = delr1z*(frad1+csfac1)-delr2z*facang12;                               \
}

#ifdef ONETYPE
#define sw_cut_ij sw_cut
#define sw_cut_ik sw_cut
#define sw_sigma_gamma_ij sw_sigma_gamma
#define sw_sigma_gamma_ik sw_sigma_gamma
#endif

__kernel void k_sw_three_center(const __global numtyp4 *restrict x_,
                                const __global numtyp2 *restrict cut_sig_gamma,
                                const __global numtyp2 *restrict sw_pre3,
                                const int ntypes,
                                const __global int * dev_nbor,
                                __global acctyp4 *restrict ans,
                                __global acctyp *restrict engv,
                                const int eflag, const int vflag,
                                const int inum,  const int nbor_pitch,
                                const int t_per_atom, const int evatom) {
  int n_stride;
  const int tpa_sq=fast_mul(t_per_atom,t_per_atom);
  local_allocate_store_three();

  int tid, ii, offset;
  atom_info(tpa_sq,ii,tid,offset);

  #ifdef ONETYPE
  const numtyp sw_cut=cut_sig_gamma[ONETYPE].x;
  const numtyp sw_sigma_gamma=cut_sig_gamma[ONETYPE].y;
  const numtyp sw_lambda_epsilon_ijk=sw_pre3[ONETYPE3].x;
  const numtyp sw_costheta_ijk=sw_pre3[ONETYPE3].y;
  #endif

  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp energy, virial[6];
  if (EVFLAG) {
    energy=(acctyp)0;
    for (int i=0; i<6; i++) virial[i]=(acctyp)0;
  }

  if (ii<inum) {
    int i, numj, nbor_j, nbor_end;
    int offset_j=offset/t_per_atom;
    nbor_info_p(dev_nbor,nbor_pitch,t_per_atom,ii,offset_j,i,numj,
                n_stride,nbor_end,nbor_j);
    int offset_k=tid & (t_per_atom-1);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    #ifndef ONETYPE
    int itype=ix.w;
    itype*=ntypes;
    #endif

    for ( ; nbor_j<nbor_end; nbor_j+=n_stride) {

      int j=dev_nbor[nbor_j];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      #ifndef ONETYPE
      int mtypej=jx.w;
      mtypej+=itype;
      #endif

      // Compute r12
      numtyp delr1x = jx.x-ix.x;
      numtyp delr1y = jx.y-ix.y;
      numtyp delr1z = jx.z-ix.z;
      numtyp rsq1 = delr1x*delr1x+delr1y*delr1y+delr1z*delr1z;

      #ifndef ONETYPE
      const numtyp sw_cut_ij=cut_sig_gamma[mtypej].x;
      const numtyp sw_sigma_gamma_ij=cut_sig_gamma[mtypej].y;
      #endif

      int nbor_k;
      nbor_k = nbor_j-offset_j+offset_k;
      if (nbor_k<=nbor_j) nbor_k += n_stride;

      for ( ; nbor_k<nbor_end; nbor_k+=n_stride) {
        int k=dev_nbor[nbor_k];
        k &= NEIGHMASK;

        numtyp4 kx; fetch4(kx,k,pos_tex);
        #ifndef ONETYPE
        const int ktype=kx.w;
        const int mtypek=itype+ktype;
        #endif

        numtyp delr2x = kx.x-ix.x;
        numtyp delr2y = kx.y-ix.y;
        numtyp delr2z = kx.z-ix.z;
        numtyp rsq2 = delr2x*delr2x + delr2y*delr2y + delr2z*delr2z;
        #ifndef ONETYPE
        const numtyp sw_cut_ik=cut_sig_gamma[mtypek].x;
        const numtyp sw_sigma_gamma_ik=cut_sig_gamma[mtypek].y;
        const int mtypejk=ntypes*mtypej+ktype;
        const numtyp sw_lambda_epsilon_ijk=sw_pre3[mtypejk].x;
        const numtyp sw_costheta_ijk=sw_pre3[mtypejk].y;
        #endif

        numtyp fjx, fjy, fjz, fkx, fky, fkz;
        threebody(delr1x,delr1y,delr1z,delr2x,delr2y,delr2z,eflag,energy);

        f.x -= fjx + fkx;
        f.y -= fjy + fky;
        f.z -= fjz + fkz;
      }
    } // for nbor

    if (EVFLAG) {
      numtyp pre;
      if (evatom==1)
        pre=THIRD;
      else
        pre=(numtyp)2.0;
      energy*=pre;
      if (vflag)
      for (int i=0; i<6; i++)
        virial[i]*=pre;
    }
  } // if ii
  store_answers(f,energy,virial,ii,inum,tid,tpa_sq,offset,
                eflag,vflag,ans,engv);
}

__kernel void k_sw_three_end(const __global numtyp4 *restrict x_,
                             const __global numtyp2 *restrict cut_sig_gamma,
                             const __global numtyp2 *restrict sw_pre3,
                             const int ntypes, const __global int * dev_nbor,
                             const __global int * dev_ilist,
                             __global acctyp4 *restrict ans,
                             __global acctyp *restrict engv,
                             const int eflag, const int vflag,
                             const int inum,  const int nbor_pitch,
                             const int t_per_atom, const int gpu_nbor) {
  int n_stride;
  const int tpa_sq=fast_mul(t_per_atom,t_per_atom);
  int tid, ii, offset;
  atom_info(tpa_sq,ii,tid,offset);

  local_allocate_store_three();

  #ifdef ONETYPE
  const numtyp sw_cut=cut_sig_gamma[ONETYPE].x;
  const numtyp sw_sigma_gamma=cut_sig_gamma[ONETYPE].y;
  const numtyp sw_lambda_epsilon_ijk=sw_pre3[ONETYPE3].x;
  const numtyp sw_costheta_ijk=sw_pre3[ONETYPE3].y;
  #endif

  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp energy, virial[6];
  if (EVFLAG) {
    energy=(acctyp)0;
    for (int i=0; i<6; i++) virial[i]=(acctyp)0;
  }

  if (ii<inum) {
    int i, numj, nbor_j, nbor_end, k_end;
    int offset_j=offset/t_per_atom;
    nbor_info_p(dev_nbor,nbor_pitch,t_per_atom,ii,offset_j,i,numj,
                n_stride,nbor_end,nbor_j);
    int offset_k=tid & (t_per_atom-1);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    #ifndef ONETYPE
    int itype=ix.w;
    itype*=ntypes;
    #endif

    for ( ; nbor_j<nbor_end; nbor_j+=n_stride) {
      int j=dev_nbor[nbor_j];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      #ifndef ONETYPE
      const int jtype=jx.w;
      const int mtypej=itype+jtype;
      #endif

      // Compute r12
      numtyp delr1x = ix.x-jx.x;
      numtyp delr1y = ix.y-jx.y;
      numtyp delr1z = ix.z-jx.z;
      numtyp rsq1 = delr1x*delr1x+delr1y*delr1y+delr1z*delr1z;

      #ifndef ONETYPE
      const numtyp sw_cut_ij=cut_sig_gamma[mtypej].x;
      const numtyp sw_sigma_gamma_ij=cut_sig_gamma[mtypej].y;
      #endif

      int nbor_k;
      if (gpu_nbor) nbor_k=j+nbor_pitch;
      else nbor_k=dev_ilist[j]+nbor_pitch;
      const int numk=dev_nbor[nbor_k];
      nbor_k+=nbor_pitch+fast_mul(j,t_per_atom-1);
      k_end=nbor_k+fast_mul(numk/t_per_atom,n_stride)+(numk&(t_per_atom-1));
      nbor_k+=offset_k;

      for ( ; nbor_k<k_end; nbor_k+=n_stride) {
        int k=dev_nbor[nbor_k];
        k &= NEIGHMASK;

        if (k == i) continue;

        numtyp4 kx; fetch4(kx,k,pos_tex);
        #ifndef ONETYPE
        const int ktype=kx.w;
        const int mtypek=jtype*ntypes+ktype;
        #endif

        numtyp delr2x = kx.x - jx.x;
        numtyp delr2y = kx.y - jx.y;
        numtyp delr2z = kx.z - jx.z;
        numtyp rsq2 = delr2x*delr2x + delr2y*delr2y + delr2z*delr2z;

        #ifndef ONETYPE
        const numtyp sw_cut_ik=cut_sig_gamma[mtypek].x;
        const numtyp sw_sigma_gamma_ik=cut_sig_gamma[mtypek].y;
        const int mtypejik=jtype*ntypes*ntypes+itype+ktype;
        const numtyp sw_lambda_epsilon_ijk=sw_pre3[mtypejik].x;
        const numtyp sw_costheta_ijk=sw_pre3[mtypejik].y;
        #endif

        numtyp fjx, fjy, fjz;
        threebody_half(delr1x,delr1y,delr1z,delr2x,delr2y,delr2z);

        f.x += fjx;
        f.y += fjy;
        f.z += fjz;
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

__kernel void k_sw_three_end_vatom(const __global numtyp4 *restrict x_,
                             const __global numtyp2 *restrict cut_sig_gamma,
                             const __global numtyp2 *restrict sw_pre3,
                             const int ntypes, const __global int * dev_nbor,
                             const __global int * dev_ilist,
                             __global acctyp4 *restrict ans,
                             __global acctyp *restrict engv,
                             const int eflag, const int vflag,
                             const int inum,  const int nbor_pitch,
                             const int t_per_atom, const int gpu_nbor) {
  int n_stride;
  const int tpa_sq=fast_mul(t_per_atom,t_per_atom);
  int tid, ii, offset;
  atom_info(tpa_sq,ii,tid,offset);

  local_allocate_store_three();

  #ifdef ONETYPE
  const numtyp sw_cut=cut_sig_gamma[ONETYPE].x;
  const numtyp sw_sigma_gamma=cut_sig_gamma[ONETYPE].y;
  const numtyp sw_lambda_epsilon_ijk=sw_pre3[ONETYPE3].x;
  const numtyp sw_costheta_ijk=sw_pre3[ONETYPE3].y;
  #endif

  acctyp4 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp energy, virial[6];
  if (EVFLAG) {
    energy=(acctyp)0;
    for (int i=0; i<6; i++) virial[i]=(acctyp)0;
  }

  if (ii<inum) {
    int i, numj, nbor_j, nbor_end, k_end;
    int offset_j=offset/t_per_atom;
    nbor_info_p(dev_nbor,nbor_pitch,t_per_atom,ii,offset_j,i,numj,
                n_stride,nbor_end,nbor_j);
    int offset_k=tid & (t_per_atom-1);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];
    #ifndef ONETYPE
    int itype=ix.w;
    itype*=ntypes;
    #endif

    for ( ; nbor_j<nbor_end; nbor_j+=n_stride) {
      int j=dev_nbor[nbor_j];
      j &= NEIGHMASK;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];
      #ifndef ONETYPE
      const int jtype=jx.w;
      const int mtypej=itype+jtype;
      #endif

      // Compute r12
      numtyp delr1x = ix.x-jx.x;
      numtyp delr1y = ix.y-jx.y;
      numtyp delr1z = ix.z-jx.z;
      numtyp rsq1 = delr1x*delr1x+delr1y*delr1y+delr1z*delr1z;

      #ifndef ONETYPE
      const numtyp sw_cut_ij=cut_sig_gamma[mtypej].x;
      const numtyp sw_sigma_gamma_ij=cut_sig_gamma[mtypej].y;
      #endif

      int nbor_k;
      if (gpu_nbor) nbor_k=j+nbor_pitch;
      else nbor_k=dev_ilist[j]+nbor_pitch;
      const int numk=dev_nbor[nbor_k];
      nbor_k+=nbor_pitch+fast_mul(j,t_per_atom-1);
      k_end=nbor_k+fast_mul(numk/t_per_atom,n_stride)+(numk&(t_per_atom-1));
      nbor_k+=offset_k;

      for ( ; nbor_k<k_end; nbor_k+=n_stride) {
        int k=dev_nbor[nbor_k];
        k &= NEIGHMASK;

        if (k == i) continue;

        numtyp4 kx; fetch4(kx,k,pos_tex);
        #ifndef ONETYPE
        const int ktype=kx.w;
        const int mtypek=jtype*ntypes+ktype;
        #endif

        numtyp delr2x = kx.x - jx.x;
        numtyp delr2y = kx.y - jx.y;
        numtyp delr2z = kx.z - jx.z;
        numtyp rsq2 = delr2x*delr2x + delr2y*delr2y + delr2z*delr2z;

        #ifndef ONETYPE
        const numtyp sw_cut_ik=cut_sig_gamma[mtypek].x;
        const numtyp sw_sigma_gamma_ik=cut_sig_gamma[mtypek].y;
        const int mtypejik=jtype*ntypes*ntypes+itype+ktype;
        const numtyp sw_lambda_epsilon_ijk=sw_pre3[mtypejik].x;
        const numtyp sw_costheta_ijk=sw_pre3[mtypejik].y;
        #endif

        numtyp fjx, fjy, fjz, fkx, fky, fkz;
        threebody(delr1x,delr1y,delr1z,delr2x,delr2y,delr2z,eflag,energy);

        f.x += fjx;
        f.y += fjy;
        f.z += fjz;
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
