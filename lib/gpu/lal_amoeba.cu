// **************************************************************************
//                                   amoeba.cu
//                             -------------------
//                          Trung Dac Nguyen (Northwestern)
//
//  Device code for acceleration of the amoeba pair style
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                :
//    email                : trung.nguyen@northwestern.edu
// ***************************************************************************

#if defined(NV_KERNEL) || defined(USE_HIP)

#include "lal_aux_fun1.h"
#ifdef LAMMPS_SMALLBIG
#define tagint int
#endif
#ifdef LAMMPS_BIGBIG
#include "inttypes.h"
#define tagint int64_t
#endif
#ifdef LAMMPS_SMALLSMALL
#define tagint int
#endif
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
#ifdef LAMMPS_SMALLBIG
#define tagint int
#endif
#ifdef LAMMPS_BIGBIG
#define tagint long
#endif
#ifdef LAMMPS_SMALLSMALL
#define tagint int
#endif

#endif // defined(NV_KERNEL) || defined(USE_HIP)


#if (SHUFFLE_AVAIL == 0)

#define local_allocate_store_ufld()                                         \
    __local acctyp red_acc[6][BLOCK_PAIR];

#define store_answers_amoeba_tq(tq, ii, inum,tid, t_per_atom, offset, i,    \
                                tep)                                        \
  if (t_per_atom>1) {                                                       \
    red_acc[0][tid]=tq.x;                                                   \
    red_acc[1][tid]=tq.y;                                                   \
    red_acc[2][tid]=tq.z;                                                   \
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                         \
      simdsync();                                                           \
      if (offset < s) {                                                     \
        for (int r=0; r<3; r++)                                             \
          red_acc[r][tid] += red_acc[r][tid+s];                             \
      }                                                                     \
    }                                                                       \
    tq.x=red_acc[0][tid];                                                   \
    tq.y=red_acc[1][tid];                                                   \
    tq.z=red_acc[2][tid];                                                   \
  }                                                                         \
  if (offset==0 && ii<inum) {                                               \
    tep[i]=tq;                                                              \
  }

#define store_answers_tep(ufld, dufld, ii, inum,tid, t_per_atom, offset,    \
                          i, tep)                                           \
  if (t_per_atom>1) {                                                       \
    red_acc[0][tid]=ufld[0];                                                \
    red_acc[1][tid]=ufld[1];                                                \
    red_acc[2][tid]=ufld[2];                                                \
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                         \
      simdsync();                                                           \
      if (offset < s) {                                                     \
        for (int r=0; r<3; r++)                                             \
          red_acc[r][tid] += red_acc[r][tid+s];                             \
      }                                                                     \
    }                                                                       \
    ufld[0]=red_acc[0][tid];                                                \
    ufld[1]=red_acc[1][tid];                                                \
    ufld[2]=red_acc[2][tid];                                                \
    red_acc[0][tid]=dufld[0];                                               \
    red_acc[1][tid]=dufld[1];                                               \
    red_acc[2][tid]=dufld[2];                                               \
    red_acc[3][tid]=dufld[3];                                               \
    red_acc[4][tid]=dufld[4];                                               \
    red_acc[5][tid]=dufld[5];                                               \
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                         \
      simdsync();                                                           \
      if (offset < s) {                                                     \
        for (int r=0; r<6; r++)                                             \
          red_acc[r][tid] += red_acc[r][tid+s];                             \
      }                                                                     \
    }                                                                       \
    dufld[0]=red_acc[0][tid];                                               \
    dufld[1]=red_acc[1][tid];                                               \
    dufld[2]=red_acc[2][tid];                                               \
    dufld[3]=red_acc[3][tid];                                               \
    dufld[4]=red_acc[4][tid];                                               \
    dufld[5]=red_acc[5][tid];                                               \
  }                                                                         \
  if (offset==0 && ii<inum) {                                               \
    acctyp3 t;                                                              \
    t.x = diz*ufld[1] - diy*ufld[2] + qixz*dufld[1] - qixy*dufld[3] +       \
      (numtyp)2.0*qiyz*(dufld[2]-dufld[5]) + (qizz-qiyy)*dufld[4];          \
    t.y = dix*ufld[2] - diz*ufld[0] - qiyz*dufld[1] + qixy*dufld[4] +       \
      (numtyp)2.0*qixz*(dufld[5]-dufld[0]) + (qixx-qizz)*dufld[3];          \
    t.z = diy*ufld[0] - dix*ufld[1] + qiyz*dufld[3] - qixz*dufld[4] +       \
      (numtyp)2.0*qixy*(dufld[0]-dufld[2]) + (qiyy-qixx)*dufld[1];          \
    tep[i]=t;                                                               \
  }

#define store_answers_fieldp(_fieldp, ii, inum,tid, t_per_atom, offset, i,  \
                              fieldp)                                       \
  if (t_per_atom>1) {                                                       \
    red_acc[0][tid]=_fieldp[0];                                             \
    red_acc[1][tid]=_fieldp[1];                                             \
    red_acc[2][tid]=_fieldp[2];                                             \
    red_acc[3][tid]=_fieldp[3];                                             \
    red_acc[4][tid]=_fieldp[4];                                             \
    red_acc[5][tid]=_fieldp[5];                                             \
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                         \
      simdsync();                                                           \
      if (offset < s) {                                                     \
        for (int r=0; r<6; r++)                                             \
          red_acc[r][tid] += red_acc[r][tid+s];                             \
      }                                                                     \
    }                                                                       \
    _fieldp[0]=red_acc[0][tid];                                             \
    _fieldp[1]=red_acc[1][tid];                                             \
    _fieldp[2]=red_acc[2][tid];                                             \
    _fieldp[3]=red_acc[3][tid];                                             \
    _fieldp[4]=red_acc[4][tid];                                             \
    _fieldp[5]=red_acc[5][tid];                                             \
  }                                                                         \
  if (offset==0 && ii<inum) {                                               \
    acctyp3 f, fp;                                                          \
    f.x = _fieldp[0];                                                       \
    f.y = _fieldp[1];                                                       \
    f.z = _fieldp[2];                                                       \
    fieldp[ii] = f;                                                         \
    fp.x = _fieldp[3];                                                      \
    fp.y = _fieldp[4];                                                      \
    fp.z = _fieldp[5];                                                      \
    fieldp[ii+inum] = fp;                                                   \
  }

#define store_answers_acc(f,energy,e_coul, virial, ii, inum, tid, t_per_atom, \
                        offset, eflag, vflag, ans, engv, ev_stride)         \
  if (t_per_atom>1) {                                                       \
    simd_reduce_add3(t_per_atom, red_acc, offset, tid, f.x, f.y, f.z);      \
    if (EVFLAG && (vflag==2 || eflag==2)) {                                 \
      if (eflag) {                                                          \
        simdsync();                                                         \
        simd_reduce_add2(t_per_atom, red_acc, offset, tid, energy, e_coul); \
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
        block_reduce_add2(simd_size(), red_acc, tid, energy, e_coul);       \
        if (vflag) __syncthreads();                                         \
        if (tid==0) {                                                       \
          engv[ei]+=energy*(acctyp)0.5;                                     \
          ei+=ev_stride;                                                    \
          engv[ei]+=e_coul*(acctyp)0.5;                                     \
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
        engv[ei]+=e_coul*(acctyp)0.5;                                       \
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

#else // SHUFFLE_AVAIL == 1

#define local_allocate_store_ufld()

#define store_answers_amoeba_tq(tq, ii, inum,tid, t_per_atom, offset, i,    \
                          tep)                                              \
  if (t_per_atom>1) {                                                       \
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                         \
      tq.x += shfl_down(tq.x, s, t_per_atom);                               \
      tq.y += shfl_down(tq.y, s, t_per_atom);                               \
      tq.z += shfl_down(tq.z, s, t_per_atom);                               \
    }                                                                       \
  }                                                                         \
  if (offset==0 && ii<inum) {                                               \
    tep[i]=tq;                                                              \
  }

#define store_answers_tep(ufld, dufld, ii, inum,tid, t_per_atom, offset,    \
                          i, tep)                                           \
  if (t_per_atom>1) {                                                       \
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                         \
      ufld[0] += shfl_down(ufld[0], s, t_per_atom);                         \
      ufld[1] += shfl_down(ufld[1], s, t_per_atom);                         \
      ufld[2] += shfl_down(ufld[2], s, t_per_atom);                         \
      dufld[0] += shfl_down(dufld[0], s, t_per_atom);                       \
      dufld[1] += shfl_down(dufld[1], s, t_per_atom);                       \
      dufld[2] += shfl_down(dufld[2], s, t_per_atom);                       \
      dufld[3] += shfl_down(dufld[3], s, t_per_atom);                       \
      dufld[4] += shfl_down(dufld[4], s, t_per_atom);                       \
      dufld[5] += shfl_down(dufld[5], s, t_per_atom);                       \
    }                                                                       \
  }                                                                         \
  if (offset==0 && ii<inum) {                                               \
    acctyp3 t;                                                              \
    t.x = diz*ufld[1] - diy*ufld[2] + qixz*dufld[1] - qixy*dufld[3] +       \
      (numtyp)2.0*qiyz*(dufld[2]-dufld[5]) + (qizz-qiyy)*dufld[4];          \
    t.y = dix*ufld[2] - diz*ufld[0] - qiyz*dufld[1] + qixy*dufld[4] +       \
      (numtyp)2.0*qixz*(dufld[5]-dufld[0]) + (qixx-qizz)*dufld[3];          \
    t.z = diy*ufld[0] - dix*ufld[1] + qiyz*dufld[3] - qixz*dufld[4] +       \
      (numtyp)2.0*qixy*(dufld[0]-dufld[2]) + (qiyy-qixx)*dufld[1];          \
    tep[i]=t;                                                               \
  }

#define store_answers_fieldp(_fieldp, ii, inum, tid, t_per_atom, offset, i, \
                             fieldp)                                        \
  if (t_per_atom>1) {                                                       \
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                         \
      _fieldp[0] += shfl_down(_fieldp[0], s, t_per_atom);                   \
      _fieldp[1] += shfl_down(_fieldp[1], s, t_per_atom);                   \
      _fieldp[2] += shfl_down(_fieldp[2], s, t_per_atom);                   \
      _fieldp[3] += shfl_down(_fieldp[3], s, t_per_atom);                   \
      _fieldp[4] += shfl_down(_fieldp[4], s, t_per_atom);                   \
      _fieldp[5] += shfl_down(_fieldp[5], s, t_per_atom);                   \
    }                                                                       \
  }                                                                         \
  if (offset==0 && ii<inum) {                                               \
    acctyp3 f, fp;                                                          \
    f.x = _fieldp[0];                                                       \
    f.y = _fieldp[1];                                                       \
    f.z = _fieldp[2];                                                       \
    fieldp[ii] = f;                                                         \
    fp.x = _fieldp[3];                                                      \
    fp.y = _fieldp[4];                                                      \
    fp.z = _fieldp[5];                                                      \
    fieldp[ii+inum] = fp;                                                   \
  }

#if (EVFLAG == 1)

#define store_answers_acc(f,energy,e_coul, virial, ii, inum, tid, t_per_atom, \
                        offset, eflag, vflag, ans, engv, ev_stride)         \
  if (t_per_atom>1) {                                                       \
    simd_reduce_add3(t_per_atom, f.x, f.y, f.z);                            \
    if (vflag==2 || eflag==2) {                                             \
      if (eflag)                                                            \
        simd_reduce_add2(t_per_atom,energy,e_coul);                         \
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
            simd_reduce_add2(vwidth, energy, e_coul);                       \
            if (voffset==0) {                                               \
              red_acc[6][bnum] = energy;                                    \
              red_acc[7][bnum] = e_coul;                                    \
            }                                                               \
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
          if (eflag) {                                                      \
            energy = red_acc[6][tid];                                       \
            e_coul = red_acc[7][tid];                                       \
          }                                                                 \
          if (vflag)                                                        \
            for (int r = 0; r < 6; r++) virial[r] = red_acc[r][tid];        \
        } else {                                                            \
          if (eflag) energy = e_coul = (acctyp)0;                           \
          if (vflag) for (int r = 0; r < 6; r++) virial[r] = (acctyp)0;     \
        }                                                                   \
      }                                                                     \
                                                                            \
      if (bnum == 0) {                                                      \
        int ei=BLOCK_ID_X;                                                  \
        if (eflag) {                                                        \
          simd_reduce_add2(vwidth, energy, e_coul);                         \
          if (tid==0) {                                                     \
            engv[ei]+=energy*(acctyp)0.5;                                   \
            ei+=ev_stride;                                                  \
            engv[ei]+=e_coul*(acctyp)0.5;                                   \
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
        engv[ei]+=e_coul*(acctyp)0.5;                                       \
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

// EVFLAG == 0
#else

#define store_answers_acc(f,energy,e_coul, virial, ii, inum, tid, t_per_atom, \
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

#endif // EVFLAG
#endif // SHUFFLE_AVAIL

#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define MY_PIS (acctyp)1.77245385090551602729

/* ----------------------------------------------------------------------
   multipole_real = real-space portion of multipole
   adapted from Tinker emreal1d() routine
------------------------------------------------------------------------- */

__kernel void k_amoeba_multipole(const __global numtyp4 *restrict x_,
                                 const __global numtyp4 *restrict extra,
                                 const __global numtyp4 *restrict coeff,
                                 const __global numtyp4 *restrict sp_amoeba,
                                 const __global int *dev_nbor,
                                 const __global int *dev_packed,
                                 const __global int *dev_short_nbor,
                                 __global acctyp3 *restrict ans,
                                 __global acctyp *restrict engv,
                                 __global acctyp3 *restrict tep,
                                 const int eflag, const int vflag, const int inum,
                                 const int nall, const int nbor_pitch,
                                 const int t_per_atom, const numtyp aewald,
                                 const numtyp felec, const numtyp off2,
                                 const numtyp polar_dscale, const numtyp polar_uscale)
{
  int tid, ii, offset, i;
  atom_info(t_per_atom,ii,tid,offset);

  int n_stride;
  local_allocate_store_charge();

  acctyp3 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp energy, e_coul, virial[6];
  if (EVFLAG) {
    energy=(acctyp)0;
    e_coul=(acctyp)0;
    for (int l=0; l<6; l++) virial[l]=(acctyp)0;
  }

  acctyp3 tq;
  tq.x=(acctyp)0; tq.y=(acctyp)0; tq.z=(acctyp)0;

  const __global numtyp4* polar1 = &extra[0];
  const __global numtyp4* polar2 = &extra[nall];
  const __global numtyp4* polar3 = &extra[2*nall];

  if (ii<inum) {
    int numj, nbor, nbor_end;
    const __global int* nbor_mem=dev_packed;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];

    // recalculate numj and nbor_end for use of the short nbor list
    if (dev_packed==dev_nbor) {
      numj = dev_short_nbor[nbor];
      nbor += n_stride;
      nbor_end = nbor+fast_mul(numj,n_stride);
      nbor_mem = dev_short_nbor;
    }

    const numtyp4 pol1i = polar1[i];
    numtyp ci  = pol1i.x;    // rpole[i][0];
    numtyp dix = pol1i.y;    // rpole[i][1];
    numtyp diy = pol1i.z;    // rpole[i][2];
    numtyp diz = pol1i.w;    // rpole[i][3];
    const numtyp4 pol2i = polar2[i];
    numtyp qixx = pol2i.x;   // rpole[i][4];
    numtyp qixy = pol2i.y;   // rpole[i][5];
    numtyp qixz = pol2i.z;   // rpole[i][6];
    numtyp qiyy = pol2i.w;   // rpole[i][8];
    const numtyp4 pol3i = polar3[i];
    numtyp qiyz = pol3i.x;   // rpole[i][9];
    numtyp qizz = pol3i.y;   // rpole[i][12];

    for ( ; nbor<nbor_end; nbor+=n_stride) {

      int jextra=nbor_mem[nbor];
      int j = jextra & NEIGHMASK15;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];

      // Compute r12
      numtyp xr = jx.x - ix.x;
      numtyp yr = jx.y - ix.y;
      numtyp zr = jx.z - ix.z;
      numtyp r2 = xr*xr + yr*yr + zr*zr;

      numtyp r = ucl_sqrt(r2);
      const numtyp4 pol1j = polar1[j];
      numtyp ck  = pol1j.x;  // rpole[j][0];
      numtyp dkx = pol1j.y;  // rpole[j][1];
      numtyp dky = pol1j.z;  // rpole[j][2];
      numtyp dkz = pol1j.w;  // rpole[j][3];
      const numtyp4 pol2j = polar2[j];
      numtyp qkxx = pol2j.x; // rpole[j][4];
      numtyp qkxy = pol2j.y; // rpole[j][5];
      numtyp qkxz = pol2j.z; // rpole[j][6];
      numtyp qkyy = pol2j.w; // rpole[j][8];
      const numtyp4 pol3j = polar3[j];
      numtyp qkyz = pol3j.x; // rpole[j][9];
      numtyp qkzz = pol3j.y; // rpole[j][12];
      //int jtype = pol3j.z; // amtype[j];
      //int jgroup =  pol3j.w; // amgroup[j];

      const numtyp4 sp_pol = sp_amoeba[sbmask15(jextra)];
      numtyp factor_mpole = sp_pol.w; // sp_mpole[sbmask15(jextra)];

      // intermediates involving moments and separation distance

      numtyp dir = dix*xr + diy*yr + diz*zr;
      numtyp qix = qixx*xr + qixy*yr + qixz*zr;
      numtyp qiy = qixy*xr + qiyy*yr + qiyz*zr;
      numtyp qiz = qixz*xr + qiyz*yr + qizz*zr;
      numtyp qir = qix*xr + qiy*yr + qiz*zr;
      numtyp dkr = dkx*xr + dky*yr + dkz*zr;
      numtyp qkx = qkxx*xr + qkxy*yr + qkxz*zr;
      numtyp qky = qkxy*xr + qkyy*yr + qkyz*zr;
      numtyp qkz = qkxz*xr + qkyz*yr + qkzz*zr;
      numtyp qkr = qkx*xr + qky*yr + qkz*zr;

      numtyp dik = dix*dkx + diy*dky + diz*dkz;
      numtyp qik = qix*qkx + qiy*qky + qiz*qkz;
      numtyp diqk = dix*qkx + diy*qky + diz*qkz;
      numtyp dkqi = dkx*qix + dky*qiy + dkz*qiz;
      numtyp qiqk = (numtyp)2.0*(qixy*qkxy+qixz*qkxz+qiyz*qkyz) +
        qixx*qkxx + qiyy*qkyy + qizz*qkzz;

      // additional intermediates involving moments and distance

      numtyp dirx = diy*zr - diz*yr;
      numtyp diry = diz*xr - dix*zr;
      numtyp dirz = dix*yr - diy*xr;
      numtyp dikx = diy*dkz - diz*dky;
      numtyp diky = diz*dkx - dix*dkz;
      numtyp dikz = dix*dky - diy*dkx;
      numtyp qirx = qiz*yr - qiy*zr;
      numtyp qiry = qix*zr - qiz*xr;
      numtyp qirz = qiy*xr - qix*yr;
      numtyp qikx = qky*qiz - qkz*qiy;
      numtyp qiky = qkz*qix - qkx*qiz;
      numtyp qikz = qkx*qiy - qky*qix;
      numtyp qixk = qixx*qkx + qixy*qky + qixz*qkz;
      numtyp qiyk = qixy*qkx + qiyy*qky + qiyz*qkz;
      numtyp qizk = qixz*qkx + qiyz*qky + qizz*qkz;
      numtyp qkxi = qkxx*qix + qkxy*qiy + qkxz*qiz;
      numtyp qkyi = qkxy*qix + qkyy*qiy + qkyz*qiz;
      numtyp qkzi = qkxz*qix + qkyz*qiy + qkzz*qiz;
      numtyp qikrx = qizk*yr - qiyk*zr;
      numtyp qikry = qixk*zr - qizk*xr;
      numtyp qikrz = qiyk*xr - qixk*yr;
      numtyp diqkx = dix*qkxx + diy*qkxy + diz*qkxz;
      numtyp diqky = dix*qkxy + diy*qkyy + diz*qkyz;
      numtyp diqkz = dix*qkxz + diy*qkyz + diz*qkzz;
      numtyp dkqix = dkx*qixx + dky*qixy + dkz*qixz;
      numtyp dkqiy = dkx*qixy + dky*qiyy + dkz*qiyz;
      numtyp dkqiz = dkx*qixz + dky*qiyz + dkz*qizz;
      numtyp dkqirx = dkqiz*yr - dkqiy*zr;
      numtyp dkqiry = dkqix*zr - dkqiz*xr;
      numtyp dkqirz = dkqiy*xr - dkqix*yr;
      numtyp dqikx = diy*qkz - diz*qky + dky*qiz - dkz*qiy -
        (numtyp)2.0*(qixy*qkxz+qiyy*qkyz+qiyz*qkzz - qixz*qkxy-qiyz*qkyy-qizz*qkyz);
      numtyp dqiky = diz*qkx - dix*qkz + dkz*qix - dkx*qiz -
        (numtyp)2.0*(qixz*qkxx+qiyz*qkxy+qizz*qkxz - qixx*qkxz-qixy*qkyz-qixz*qkzz);
      numtyp dqikz = dix*qky - diy*qkx + dkx*qiy - dky*qix -
        (numtyp)2.0*(qixx*qkxy+qixy*qkyy+qixz*qkyz - qixy*qkxx-qiyy*qkxy-qiyz*qkxz);

      // get reciprocal distance terms for this interaction

      numtyp rinv = ucl_recip(r);
      numtyp r2inv = rinv*rinv;
      numtyp rr1 = felec * rinv;
      numtyp rr3 = rr1 * r2inv;
      numtyp rr5 = (numtyp)3.0 * rr3 * r2inv;
      numtyp rr7 = (numtyp)5.0 * rr5 * r2inv;
      numtyp rr9 = (numtyp)7.0 * rr7 * r2inv;
      numtyp rr11 = (numtyp)9.0 * rr9 * r2inv;

      // calculate the real space Ewald error function terms

      numtyp ralpha = aewald * r;
      numtyp exp2a = ucl_exp(-ralpha*ralpha);
      numtyp bn[6];
      bn[0] = ucl_erfc(ralpha) * rinv;

      numtyp alsq2 = (numtyp)2.0 * aewald*aewald;
      numtyp alsq2n = (numtyp)0.0;
      if (aewald > (numtyp)0.0) alsq2n = (numtyp)1.0 / (MY_PIS*aewald);

      int m;
      for (m = 1; m < 6; m++) {
        numtyp bfac = (numtyp) (m+m-1);
        alsq2n = alsq2 * alsq2n;
        bn[m] = (bfac*bn[m-1]+alsq2n*exp2a) * r2inv;
      }
      for (m = 0; m < 6; m++) bn[m] *= felec;

      numtyp term1,term2,term3;
      numtyp term4,term5,term6;

      term1 = ci*ck;
      term2 = ck*dir - ci*dkr + dik;
      term3 = ci*qkr + ck*qir - dir*dkr + (numtyp)2.0*(dkqi-diqk+qiqk);
      term4 = dir*qkr - dkr*qir - (numtyp)4.0*qik;
      term5 = qir*qkr;
      numtyp scalek = (numtyp)1.0 - factor_mpole;
      rr1 = bn[0] - scalek*rr1;
      rr3 = bn[1] - scalek*rr3;
      rr5 = bn[2] - scalek*rr5;
      rr7 = bn[3] - scalek*rr7;
      rr9 = bn[4] - scalek*rr9;
      rr11 = bn[5] - scalek*rr11;
      numtyp e = term1*rr1 + term2*rr3 + term3*rr5 + term4*rr7 + term5*rr9;

      // find standard multipole intermediates for force and torque

      numtyp de = term1*rr3 + term2*rr5 + term3*rr7 + term4*rr9 + term5*rr11;
      term1 = -ck*rr3 + dkr*rr5 - qkr*rr7;
      term2 = ci*rr3 + dir*rr5 + qir*rr7;
      term3 = (numtyp)2.0 * rr5;
      term4 = (numtyp)2.0 * (-ck*rr5+dkr*rr7-qkr*rr9);
      term5 = (numtyp)2.0 * (-ci*rr5-dir*rr7-qir*rr9);
      term6 = (numtyp)4.0 * rr7;

      energy += e;

      // compute the force components for this interaction

      numtyp frcx = de*xr + term1*dix + term2*dkx + term3*(diqkx-dkqix) +
        term4*qix + term5*qkx + term6*(qixk+qkxi);
      numtyp frcy = de*yr + term1*diy + term2*dky + term3*(diqky-dkqiy) +
        term4*qiy + term5*qky + term6*(qiyk+qkyi);
      numtyp frcz = de*zr + term1*diz + term2*dkz + term3*(diqkz-dkqiz) +
        term4*qiz + term5*qkz + term6*(qizk+qkzi);

      // compute the torque components for this interaction

      numtyp ttmix = -rr3*dikx + term1*dirx + term3*(dqikx+dkqirx) -
        term4*qirx - term6*(qikrx+qikx);
      numtyp ttmiy = -rr3*diky + term1*diry + term3*(dqiky+dkqiry) -
        term4*qiry - term6*(qikry+qiky);
      numtyp ttmiz = -rr3*dikz + term1*dirz + term3*(dqikz+dkqirz) -
        term4*qirz - term6*(qikrz+qikz);

      // increment force-based gradient and torque on first site

      f.x -= frcx;
      f.y -= frcy;
      f.z -= frcz;
      tq.x += ttmix;
      tq.y += ttmiy;
      tq.z += ttmiz;

      if (EVFLAG && vflag) {
        numtyp vxx = -xr * frcx;
        numtyp vxy = (numtyp)-0.5 * (yr*frcx+xr*frcy);
        numtyp vxz = (numtyp)-0.5 * (zr*frcx+xr*frcz);
        numtyp vyy = -yr * frcy;
        numtyp vyz = (numtyp)-0.5 * (zr*frcy+yr*frcz);
        numtyp vzz = -zr * frcz;

        virial[0] -= vxx;
        virial[1] -= vyy;
        virial[2] -= vzz;
        virial[3] -= vxy;
        virial[4] -= vxz;
        virial[5] -= vyz;
      }
    } // nbor

  } // ii<inum

  // accumulate tq
  store_answers_amoeba_tq(tq,ii,inum,tid,t_per_atom,offset,i,tep);

  // accumate force, energy and virial: use _acc if not the first kernel
  store_answers_q(f,energy,e_coul,virial,ii,inum,tid,t_per_atom,
     offset,eflag,vflag,ans,engv);
}

/* ----------------------------------------------------------------------
  udirect2b = Ewald real direct field via list
  udirect2b computes the real space contribution of the permanent
   atomic multipole moments to the field via a neighbor list
------------------------------------------------------------------------- */

__kernel void k_amoeba_udirect2b(const __global numtyp4 *restrict x_,
                                 const __global numtyp4 *restrict extra,
                                 const __global numtyp4 *restrict coeff,
                                 const __global numtyp4 *restrict sp_amoeba,
                                 const __global int *dev_nbor,
                                 const __global int *dev_packed,
                                 const __global int *dev_short_nbor,
                                 __global acctyp3 *restrict fieldp,
                                 const int inum,  const int nall,
                                 const int nbor_pitch, const int t_per_atom,
                                 const numtyp aewald, const numtyp off2,
                                 const numtyp polar_dscale, const numtyp polar_uscale)
{
  int tid, ii, offset, i;
  atom_info(t_per_atom,ii,tid,offset);

  int n_stride;
  local_allocate_store_ufld();

  acctyp _fieldp[6];
  for (int l=0; l<6; l++) _fieldp[l]=(acctyp)0;

  const __global numtyp4* polar1 = &extra[0];
  const __global numtyp4* polar2 = &extra[nall];
  const __global numtyp4* polar3 = &extra[2*nall];

  if (ii<inum) {
    int numj, nbor, nbor_end;
    const __global int* nbor_mem=dev_packed;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];

    // recalculate numj and nbor_end for use of the short nbor list
    if (dev_packed==dev_nbor) {
      numj = dev_short_nbor[nbor];
      nbor += n_stride;
      nbor_end = nbor+fast_mul(numj,n_stride);
      nbor_mem = dev_short_nbor;
    }

    const numtyp4 pol3i = polar3[i];
    int itype  = pol3i.z;    // amtype[i];
    int igroup = pol3i.w;    // amgroup[i];

    numtyp pdi = coeff[itype].x;
    numtyp pti = coeff[itype].y;
    numtyp ddi = coeff[itype].z;

    numtyp aesq2 = (numtyp)2.0 * aewald*aewald;
    numtyp aesq2n = (numtyp)0.0;
    if (aewald > (numtyp)0.0) aesq2n = (numtyp)1.0 / (MY_PIS*aewald);

    for ( ; nbor<nbor_end; nbor+=n_stride) {

      int jextra=nbor_mem[nbor];
      int j = jextra & NEIGHMASK15;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];

      // Compute r12
      numtyp xr = jx.x - ix.x;
      numtyp yr = jx.y - ix.y;
      numtyp zr = jx.z - ix.z;
      numtyp r2 = xr*xr + yr*yr + zr*zr;

      numtyp r = ucl_sqrt(r2);
      numtyp rinv = ucl_rsqrt(r2);
      numtyp r2inv = rinv*rinv;
      numtyp rr1 = rinv;
      numtyp rr3 = rr1 * r2inv;
      numtyp rr5 = (numtyp)3.0 * rr3 * r2inv;
      numtyp rr7 = (numtyp)5.0 * rr5 * r2inv;

      const numtyp4 pol1j = polar1[j];
      numtyp ck  = pol1j.x;  // rpole[j][0];
      numtyp dkx = pol1j.y;  // rpole[j][1];
      numtyp dky = pol1j.z;  // rpole[j][2];
      numtyp dkz = pol1j.w;  // rpole[j][3];
      const numtyp4 pol2j = polar2[j];
      numtyp qkxx = pol2j.x; // rpole[j][4];
      numtyp qkxy = pol2j.y; // rpole[j][5];
      numtyp qkxz = pol2j.z; // rpole[j][6];
      numtyp qkyy = pol2j.w; // rpole[j][8];
      const numtyp4 pol3j = polar3[j];
      numtyp qkyz = pol3j.x; // rpole[j][9];
      numtyp qkzz = pol3j.y; // rpole[j][12];
      int jtype = pol3j.z; // amtype[j];
      int jgroup =  pol3j.w; // amgroup[j];

      numtyp factor_dscale, factor_pscale;
      const numtyp4 sp_pol = sp_amoeba[sbmask15(jextra)];
      if (igroup == jgroup) {
        factor_pscale = sp_pol.y; // sp_amoeba_piscale[sbmask15(jextra)];
        factor_dscale = polar_dscale;
      } else {
        factor_pscale = sp_pol.z; // sp_amoeba_pscale[sbmask15(jextra)];
        factor_dscale = (numtyp)1.0;
      }

      // intermediates involving moments and separation distance

      numtyp dkr = dkx*xr + dky*yr + dkz*zr;
      numtyp qkx = qkxx*xr + qkxy*yr + qkxz*zr;
      numtyp qky = qkxy*xr + qkyy*yr + qkyz*zr;
      numtyp qkz = qkxz*xr + qkyz*yr + qkzz*zr;
      numtyp qkr = qkx*xr + qky*yr + qkz*zr;

      // calculate the real space Ewald error function terms

      numtyp ralpha = aewald * r;
      numtyp exp2a = ucl_exp(-ralpha*ralpha);
      numtyp bn[4], bcn[3];
      bn[0] = ucl_erfc(ralpha) * rinv;

      numtyp aefac = aesq2n;
      for (int m = 1; m <= 3; m++) {
        numtyp bfac = (numtyp) (m+m-1);
        aefac = aesq2 * aefac;
        bn[m] = (bfac*bn[m-1]+aefac*exp2a) * r2inv;
      }

      // find the field components for Thole polarization damping

      numtyp scale3 = (numtyp)1.0;
      numtyp scale5 = (numtyp)1.0;
      numtyp scale7 = (numtyp)1.0;
      numtyp damp = pdi * coeff[jtype].x; // pdamp[jtype]
      if (damp != (numtyp)0.0) {
        numtyp pgamma = MIN(ddi,coeff[jtype].z); // dirdamp[jtype]
        if (pgamma != (numtyp)0.0) {
          numtyp tmp = r*ucl_recip(damp);
          damp = pgamma * ucl_sqrt(tmp*tmp*tmp);
          if (damp < (numtyp)50.0) {
            numtyp expdamp = ucl_exp(-damp) ;
            scale3 = (numtyp)1.0 - expdamp ;
            scale5 = (numtyp)1.0 - expdamp*((numtyp)1.0+(numtyp)0.5*damp);
            scale7 = (numtyp)1.0 - expdamp*((numtyp)1.0+(numtyp)0.65*damp + (numtyp)0.15*damp*damp);
          }
        } else {
          pgamma = MIN(pti,coeff[jtype].y); // thole[jtype]
          numtyp tmp = r*ucl_recip(damp);
          damp = pgamma * (tmp*tmp*tmp);
          if (damp < (numtyp)50.0) {
            numtyp expdamp = ucl_exp(-damp);
            scale3 = (numtyp)1.0 - expdamp;
            scale5 = (numtyp)1.0 - expdamp*((numtyp)1.0+damp);
            scale7 = (numtyp)1.0 - expdamp*((numtyp)1.0+damp + (numtyp)0.6*damp*damp);
          }
        }
      } else { // damp == 0: ???
      }

      numtyp scalek = factor_dscale;
      bcn[0] = bn[1] - ((numtyp)1.0-scalek*scale3)*rr3;
      bcn[1] = bn[2] - ((numtyp)1.0-scalek*scale5)*rr5;
      bcn[2] = bn[3] - ((numtyp)1.0-scalek*scale7)*rr7;

      numtyp fid[3];
      fid[0] = -xr*(bcn[0]*ck-bcn[1]*dkr+bcn[2]*qkr) - bcn[0]*dkx + (numtyp)2.0*bcn[1]*qkx;
      fid[1] = -yr*(bcn[0]*ck-bcn[1]*dkr+bcn[2]*qkr) - bcn[0]*dky + (numtyp)2.0*bcn[1]*qky;
      fid[2] = -zr*(bcn[0]*ck-bcn[1]*dkr+bcn[2]*qkr) - bcn[0]*dkz + (numtyp)2.0*bcn[1]*qkz;

      scalek = factor_pscale;
      bcn[0] = bn[1] - ((numtyp)1.0-scalek*scale3)*rr3;
      bcn[1] = bn[2] - ((numtyp)1.0-scalek*scale5)*rr5;
      bcn[2] = bn[3] - ((numtyp)1.0-scalek*scale7)*rr7;
      numtyp fip[3];
      fip[0] = -xr*(bcn[0]*ck-bcn[1]*dkr+bcn[2]*qkr) - bcn[0]*dkx + (numtyp)2.0*bcn[1]*qkx;
      fip[1] = -yr*(bcn[0]*ck-bcn[1]*dkr+bcn[2]*qkr) - bcn[0]*dky + (numtyp)2.0*bcn[1]*qky;
      fip[2] = -zr*(bcn[0]*ck-bcn[1]*dkr+bcn[2]*qkr) - bcn[0]*dkz + (numtyp)2.0*bcn[1]*qkz;

      _fieldp[0] += fid[0];
      _fieldp[1] += fid[1];
      _fieldp[2] += fid[2];
      _fieldp[3] += fip[0];
      _fieldp[4] += fip[1];
      _fieldp[5] += fip[2];
    }  // nbor

  } // ii<inum

  // accumulate field and fieldp

  store_answers_fieldp(_fieldp,ii,inum,tid,t_per_atom,offset,i,fieldp);
}

/* ----------------------------------------------------------------------
  umutual2b = Ewald real mutual field via list
   umutual2b computes the real space contribution of the induced
   atomic dipole moments to the field via a neighbor list
------------------------------------------------------------------------- */

__kernel void k_amoeba_umutual2b(const __global numtyp4 *restrict x_,
                                 const __global numtyp4 *restrict extra,
                                 const __global numtyp4 *restrict coeff,
                                 const __global numtyp4 *restrict sp_amoeba,
                                 const __global int *dev_nbor,
                                 const __global int *dev_packed,
                                 const __global int *dev_short_nbor,
                                 __global acctyp3 *restrict fieldp,
                                 const int inum,  const int nall,
                                 const int nbor_pitch, const int t_per_atom,
                                 const numtyp aewald, const numtyp off2,
                                 const numtyp polar_dscale, const numtyp polar_uscale)
{
  int tid, ii, offset, i;
  atom_info(t_per_atom,ii,tid,offset);

  int n_stride;
  local_allocate_store_ufld();

  acctyp _fieldp[6];
  for (int l=0; l<6; l++) _fieldp[l]=(acctyp)0;

  const __global numtyp4* polar3 = &extra[2*nall];
  const __global numtyp4* polar4 = &extra[3*nall];
  const __global numtyp4* polar5 = &extra[4*nall];

  if (ii<inum) {
    int numj, nbor, nbor_end;
    const __global int* nbor_mem=dev_packed;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];

    // recalculate numj and nbor_end for use of the short nbor list
    if (dev_packed==dev_nbor) {
      numj = dev_short_nbor[nbor];
      nbor += n_stride;
      nbor_end = nbor+fast_mul(numj,n_stride);
      nbor_mem = dev_short_nbor;
    }

    int itype,igroup;
    itype  = polar3[i].z; // amtype[i];
    igroup = polar3[i].w; // amgroup[i];

    numtyp pdi = coeff[itype].x;
    numtyp pti = coeff[itype].y;

    numtyp aesq2 = (numtyp)2.0 * aewald*aewald;
    numtyp aesq2n = (numtyp)0.0;
    if (aewald > (numtyp)0.0) aesq2n = (numtyp)1.0 / (MY_PIS*aewald);

    for ( ; nbor<nbor_end; nbor+=n_stride) {

      int jextra=nbor_mem[nbor];
      int j = jextra & NEIGHMASK15;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];

      // Compute r12
      numtyp xr = jx.x - ix.x;
      numtyp yr = jx.y - ix.y;
      numtyp zr = jx.z - ix.z;
      numtyp r2 = xr*xr + yr*yr + zr*zr;

      numtyp r = ucl_sqrt(r2);
      numtyp rinv = ucl_rsqrt(r2);
      numtyp r2inv = rinv*rinv;
      numtyp rr1 = rinv;
      numtyp rr3 = rr1 * r2inv;
      numtyp rr5 = (numtyp)3.0 * rr3 * r2inv;

      const numtyp4 pol3j = polar3[j];
      int jtype = pol3j.z; // amtype[j];
      int jgroup =  pol3j.w; // amgroup[j];
      const numtyp4 pol4j = polar4[j];
      numtyp ukx = pol4j.x;  // uind[j][0];
      numtyp uky = pol4j.y;  // uind[j][1];
      numtyp ukz = pol4j.z;  // uind[j][2];
      const numtyp4 pol5j = polar5[j];
      numtyp ukxp = pol5j.x; // uinp[j][0];
      numtyp ukyp = pol5j.y; // uinp[j][1];
      numtyp ukzp = pol5j.z; // uinp[j][2];

      numtyp factor_uscale;
      if (igroup == jgroup) factor_uscale = polar_uscale;
      else factor_uscale = (numtyp)1.0;

      // calculate the real space Ewald error function terms

      numtyp ralpha = aewald * r;
      numtyp exp2a = ucl_exp(-ralpha*ralpha);
      numtyp bn[4];
      bn[0] = ucl_erfc(ralpha) * rinv;

      numtyp aefac = aesq2n;
      for (int m = 1; m <= 3; m++) {
        numtyp bfac = (numtyp) (m+m-1);
        aefac = aesq2 * aefac;
        bn[m] = (bfac*bn[m-1]+aefac*exp2a) * r2inv;
      }

      // find terms needed later to compute mutual polarization
      // if (poltyp != DIRECT)
      numtyp scale3 = (numtyp)1.0;
      numtyp scale5 = (numtyp)1.0;
      numtyp damp = pdi * coeff[jtype].x; // pdamp[jtype]
      if (damp != (numtyp)0.0) {
        numtyp pgamma = MIN(pti,coeff[jtype].y); // thole[jtype]
        damp = pgamma * ucl_powr(r/damp,(numtyp)3.0);
        if (damp < (numtyp)50.0) {
          numtyp expdamp = ucl_exp(-damp);
          scale3 = (numtyp)1.0 - expdamp;
          scale5 = (numtyp)1.0 - expdamp*((numtyp)1.0+damp);
        }

      } else { // damp == 0: ???
      }

      numtyp scalek = factor_uscale;
      numtyp bcn[3];
      bcn[0] = bn[1] - ((numtyp)1.0-scalek*scale3)*rr3;
      bcn[1] = bn[2] - ((numtyp)1.0-scalek*scale5)*rr5;

      numtyp tdipdip[6];
      tdipdip[0] = -bcn[0] + bcn[1]*xr*xr;
      tdipdip[1] = bcn[1]*xr*yr;
      tdipdip[2] = bcn[1]*xr*zr;
      tdipdip[3] = -bcn[0] + bcn[1]*yr*yr;
      tdipdip[4] = bcn[1]*yr*zr;
      tdipdip[5] = -bcn[0] + bcn[1]*zr*zr;

      numtyp fid[3];
      fid[0] = tdipdip[0]*ukx + tdipdip[1]*uky + tdipdip[2]*ukz;
      fid[1] = tdipdip[1]*ukx + tdipdip[3]*uky + tdipdip[4]*ukz;
      fid[2] = tdipdip[2]*ukx + tdipdip[4]*uky + tdipdip[5]*ukz;

      numtyp fip[3];
      fip[0] = tdipdip[0]*ukxp + tdipdip[1]*ukyp + tdipdip[2]*ukzp;
      fip[1] = tdipdip[1]*ukxp + tdipdip[3]*ukyp + tdipdip[4]*ukzp;
      fip[2] = tdipdip[2]*ukxp + tdipdip[4]*ukyp + tdipdip[5]*ukzp;

      _fieldp[0] += fid[0];
      _fieldp[1] += fid[1];
      _fieldp[2] += fid[2];
      _fieldp[3] += fip[0];
      _fieldp[4] += fip[1];
      _fieldp[5] += fip[2];
    }  // nbor

  } // ii<inum

  // accumulate field and fieldp

  store_answers_fieldp(_fieldp,ii,inum,tid,t_per_atom,offset,i,fieldp);
}

/* ----------------------------------------------------------------------
   polar_real = real-space portion of induced dipole polarization
   adapted from Tinker epreal1d() routine
------------------------------------------------------------------------- */

__kernel void k_amoeba_polar(const __global numtyp4 *restrict x_,
                             const __global numtyp4 *restrict extra,
                             const __global numtyp4 *restrict coeff,
                             const __global numtyp4 *restrict sp_amoeba,
                             const __global int *dev_nbor,
                             const __global int *dev_packed,
                             const __global int *dev_short_nbor,
                             __global acctyp3 *restrict ans,
                             __global acctyp *restrict engv,
                             __global acctyp3 *restrict tep,
                             const int eflag, const int vflag, const int inum,
                             const int nall, const int nbor_pitch, const int t_per_atom,
                             const numtyp aewald, const numtyp felec,
                             const numtyp off2, const numtyp polar_dscale,
                             const numtyp polar_uscale)
{
  int tid, ii, offset, i;
  atom_info(t_per_atom,ii,tid,offset);

  int n_stride;
  local_allocate_store_charge();

  acctyp3 f;
  f.x=(acctyp)0; f.y=(acctyp)0; f.z=(acctyp)0;
  acctyp energy, e_coul, virial[6];
  if (EVFLAG) {
    energy=(acctyp)0;
    e_coul=(acctyp)0;
    for (int l=0; l<6; l++) virial[l]=(acctyp)0;
  }

  acctyp ufld[3];
  ufld[0] = (acctyp)0; ufld[1]=(acctyp)0; ufld[2]=(acctyp)0;
  acctyp dufld[6];
  for (int l=0; l<6; l++) dufld[l]=(acctyp)0;

  numtyp dix,diy,diz,qixx,qixy,qixz,qiyy,qiyz,qizz;

  const __global numtyp4* polar1 = &extra[0];
  const __global numtyp4* polar2 = &extra[nall];
  const __global numtyp4* polar3 = &extra[2*nall];
  const __global numtyp4* polar4 = &extra[3*nall];
  const __global numtyp4* polar5 = &extra[4*nall];

  if (ii<inum) {
    int itype,igroup;
    numtyp ci,uix,uiy,uiz,uixp,uiyp,uizp;

    int numj, nbor, nbor_end;
    const __global int* nbor_mem=dev_packed;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];

    // recalculate numj and nbor_end for use of the short nbor list
    if (dev_packed==dev_nbor) {
      numj = dev_short_nbor[nbor];
      nbor += n_stride;
      nbor_end = nbor+fast_mul(numj,n_stride);
      nbor_mem = dev_short_nbor;
    }

    const numtyp4 pol1i = polar1[i];
    ci  = pol1i.x;    // rpole[i][0];
    dix = pol1i.y;    // rpole[i][1];
    diy = pol1i.z;    // rpole[i][2];
    diz = pol1i.w;    // rpole[i][3];
    const numtyp4 pol2i = polar2[i];
    qixx = pol2i.x;   // rpole[i][4];
    qixy = pol2i.y;   // rpole[i][5];
    qixz = pol2i.z;   // rpole[i][6];
    qiyy = pol2i.w;   // rpole[i][8];
    const numtyp4 pol3i = polar3[i];
    qiyz = pol3i.x;   // rpole[i][9];
    qizz = pol3i.y;   // rpole[i][12];
    itype  = pol3i.z;    // amtype[i];
    igroup = pol3i.w;    // amgroup[i];
    const numtyp4 pol4i = polar4[i];
    uix = pol4i.x;    // uind[i][0];
    uiy = pol4i.y;    // uind[i][1];
    uiz = pol4i.z;    // uind[i][2];
    const numtyp4 pol5i = polar5[i];
    uixp = pol5i.x;   // uinp[i][0];
    uiyp = pol5i.y;   // uinp[i][1];
    uizp = pol5i.z;   // uinp[i][2];

    numtyp pdi = coeff[itype].x;
    numtyp pti = coeff[itype].y;

    for ( ; nbor<nbor_end; nbor+=n_stride) {

      int jextra=nbor_mem[nbor];
      int j = jextra & NEIGHMASK15;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];

      // Compute r12
      numtyp xr = jx.x - ix.x;
      numtyp yr = jx.y - ix.y;
      numtyp zr = jx.z - ix.z;
      numtyp r2 = xr*xr + yr*yr + zr*zr;
      numtyp r = ucl_sqrt(r2);

      const numtyp4 pol1j = polar1[j];
      numtyp ck = pol1j.x;   // rpole[j][0];
      numtyp dkx = pol1j.y;  // rpole[j][1];
      numtyp dky = pol1j.z;  // rpole[j][2];
      numtyp dkz = pol1j.w;  // rpole[j][3];
      const numtyp4 pol2j = polar2[j];
      numtyp qkxx = pol2j.x; // rpole[j][4];
      numtyp qkxy = pol2j.y; // rpole[j][5];
      numtyp qkxz = pol2j.z; // rpole[j][6];
      numtyp qkyy = pol2j.w; // rpole[j][8];
      const numtyp4 pol3j = polar3[j];
      numtyp qkyz = pol3j.x; // rpole[j][9];
      numtyp qkzz = pol3j.y; // rpole[j][12];
      int jtype =   pol3j.z; // amtype[j];
      int jgroup =  pol3j.w; // amgroup[j];
      const numtyp4 pol4j = polar4[j];
      numtyp ukx = pol4j.x;  // uind[j][0];
      numtyp uky = pol4j.y;  // uind[j][1];
      numtyp ukz = pol4j.z;  // uind[j][2];
      const numtyp4 pol5j = polar5[j];
      numtyp ukxp = pol5j.x; // uinp[j][0];
      numtyp ukyp = pol5j.y; // uinp[j][1];
      numtyp ukzp = pol5j.z; // uinp[j][2];

      numtyp factor_dscale, factor_pscale, factor_uscale;
      const numtyp4 sp_pol = sp_amoeba[sbmask15(jextra)];
      if (igroup == jgroup) {
        factor_pscale = sp_pol.y; // sp_amoeba_piscale[sbmask15(jextra)];
        factor_dscale = polar_dscale;
        factor_uscale = polar_uscale;
      } else {
        factor_pscale = sp_pol.z; // sp_amoeba_pscale[sbmask15(jextra)];
        factor_dscale = factor_uscale = (numtyp)1.0;
      }

      // intermediates involving moments and separation distance

      numtyp dir = dix*xr + diy*yr + diz*zr;
      numtyp qix = qixx*xr + qixy*yr + qixz*zr;
      numtyp qiy = qixy*xr + qiyy*yr + qiyz*zr;
      numtyp qiz = qixz*xr + qiyz*yr + qizz*zr;
      numtyp qir = qix*xr + qiy*yr + qiz*zr;
      numtyp dkr = dkx*xr + dky*yr + dkz*zr;
      numtyp qkx = qkxx*xr + qkxy*yr + qkxz*zr;
      numtyp qky = qkxy*xr + qkyy*yr + qkyz*zr;
      numtyp qkz = qkxz*xr + qkyz*yr + qkzz*zr;
      numtyp qkr = qkx*xr + qky*yr + qkz*zr;
      numtyp uir = uix*xr + uiy*yr + uiz*zr;
      numtyp ukr = ukx*xr + uky*yr + ukz*zr;
      numtyp ukrp = ukxp*xr + ukyp*yr + ukzp*zr;

      // get reciprocal distance terms for this interaction

      numtyp rinv = ucl_recip(r);
      numtyp r2inv = rinv*rinv;
      numtyp rr1 = felec * rinv;
      numtyp rr3 = rr1 * r2inv;
      numtyp rr5 = (numtyp)3.0 * rr3 * r2inv;
      numtyp rr7 = (numtyp)5.0 * rr5 * r2inv;
      numtyp rr9 = (numtyp)7.0 * rr7 * r2inv;

      // calculate the real space Ewald error function terms

      int k,m;
      numtyp psc3,psc5,psc7;
      numtyp dsc3,dsc5,dsc7;
      numtyp usc3,usc5;
      numtyp psr3,psr5,psr7;
      numtyp dsr3,dsr5,dsr7;
      numtyp usr5;
      numtyp term1,term2,term3;
      numtyp term4,term5;
      numtyp term6,term7;
      numtyp rc3[3],rc5[3],rc7[3];
      numtyp prc3[3],prc5[3],prc7[3];
      numtyp drc3[3],drc5[3],drc7[3];
      numtyp urc3[3],urc5[3];

      numtyp ralpha = aewald * r;
      numtyp exp2a = ucl_exp(-ralpha*ralpha);
      numtyp bn[5];
      bn[0] = ucl_erfc(ralpha) * rinv;

      numtyp alsq2 = (numtyp)2.0 * aewald*aewald;
      numtyp alsq2n = (numtyp)0.0;
      if (aewald > (numtyp)0.0) alsq2n = (numtyp)1.0 / (MY_PIS*aewald);

      for (m = 1; m <= 4; m++) {
        numtyp bfac = (numtyp) (m+m-1);
        alsq2n = alsq2 * alsq2n;
        bn[m] = (bfac*bn[m-1]+alsq2n*exp2a) * r2inv;
      }
      for (m = 0; m < 5; m++) bn[m] *= felec;

      // apply Thole polarization damping to scale factors

      numtyp sc3 = (numtyp)1.0;
      numtyp sc5 = (numtyp)1.0;
      numtyp sc7 = (numtyp)1.0;
      for (k = 0; k < 3; k++) {
        rc3[k] = (numtyp)0.0;
        rc5[k] = (numtyp)0.0;
        rc7[k] = (numtyp)0.0;
      }

      // apply Thole polarization damping to scale factors

      numtyp damp = pdi * coeff[jtype].x; // pdamp[jtype]
      if (damp != (numtyp)0.0) {
        numtyp pgamma = MIN(pti,coeff[jtype].y); // thole[jtype]
        numtyp tmp = r*ucl_recip(damp);
        damp = pgamma * (tmp*tmp*tmp);
        if (damp < (numtyp)50.0) {
          numtyp expdamp = ucl_exp(-damp);
          sc3 = (numtyp)1.0 - expdamp;
          sc5 = (numtyp)1.0 - ((numtyp)1.0+damp)*expdamp;
          sc7 = (numtyp)1.0 - ((numtyp)1.0+damp+(numtyp)0.6*damp*damp) * expdamp;
          numtyp temp3 = (numtyp)3.0 * damp * expdamp * r2inv;
          numtyp temp5 = damp;
          numtyp temp7 = (numtyp)-0.2 + (numtyp)0.6*damp;
          rc3[0] = xr * temp3;
          rc3[1] = yr * temp3;
          rc3[2] = zr * temp3;
          rc5[0] = rc3[0] * temp5;
          rc5[1] = rc3[1] * temp5;
          rc5[2] = rc3[2] * temp5;
          rc7[0] = rc5[0] * temp7;
          rc7[1] = rc5[1] * temp7;
          rc7[2] = rc5[2] * temp7;
        }

        psc3 = (numtyp)1.0 - sc3*factor_pscale;
        psc5 = (numtyp)1.0 - sc5*factor_pscale;
        psc7 = (numtyp)1.0 - sc7*factor_pscale;
        dsc3 = (numtyp)1.0 - sc3*factor_dscale;
        dsc5 = (numtyp)1.0 - sc5*factor_dscale;
        dsc7 = (numtyp)1.0 - sc7*factor_dscale;
        usc3 = (numtyp)1.0 - sc3*factor_uscale;
        usc5 = (numtyp)1.0 - sc5*factor_uscale;
        psr3 = bn[1] - psc3*rr3;
        psr5 = bn[2] - psc5*rr5;
        psr7 = bn[3] - psc7*rr7;
        dsr3 = bn[1] - dsc3*rr3;
        dsr5 = bn[2] - dsc5*rr5;
        dsr7 = bn[3] - dsc7*rr7;
        usr5 = bn[2] - usc5*rr5;
        for (k = 0; k < 3; k++) {
          prc3[k] = rc3[k] * factor_pscale;
          prc5[k] = rc5[k] * factor_pscale;
          prc7[k] = rc7[k] * factor_pscale;
          drc3[k] = rc3[k] * factor_dscale;
          drc5[k] = rc5[k] * factor_dscale;
          drc7[k] = rc7[k] * factor_dscale;
          urc3[k] = rc3[k] * factor_uscale;
          urc5[k] = rc5[k] * factor_uscale;
        }
      } else { // damp == 0: ???
      }

      // get the induced dipole field used for dipole torques

      numtyp tix3 = psr3*ukx + dsr3*ukxp;
      numtyp tiy3 = psr3*uky + dsr3*ukyp;
      numtyp tiz3 = psr3*ukz + dsr3*ukzp;
      numtyp tuir = -psr5*ukr - dsr5*ukrp;

      ufld[0] += tix3 + xr*tuir;
      ufld[1] += tiy3 + yr*tuir;
      ufld[2] += tiz3 + zr*tuir;

      // get induced dipole field gradient used for quadrupole torques

      numtyp tix5 = (numtyp)2.0 * (psr5*ukx+dsr5*ukxp);
      numtyp tiy5 = (numtyp)2.0 * (psr5*uky+dsr5*ukyp);
      numtyp tiz5 = (numtyp)2.0 * (psr5*ukz+dsr5*ukzp);
      tuir = -psr7*ukr - dsr7*ukrp;

      dufld[0] += xr*tix5 + xr*xr*tuir;
      dufld[1] += xr*tiy5 + yr*tix5 + (numtyp)2.0*xr*yr*tuir;
      dufld[2] += yr*tiy5 + yr*yr*tuir;
      dufld[3] += xr*tiz5 + zr*tix5 + (numtyp)2.0*xr*zr*tuir;
      dufld[4] += yr*tiz5 + zr*tiy5 + (numtyp)2.0*yr*zr*tuir;
      dufld[5] += zr*tiz5 + zr*zr*tuir;

      // get the dEd/dR terms used for direct polarization force

      term1 = bn[2] - dsc3*rr5;
      term2 = bn[3] - dsc5*rr7;
      term3 = -dsr3 + term1*xr*xr - rr3*xr*drc3[0];
      term4 = rr3*drc3[0] - term1*xr - dsr5*xr;
      term5 = term2*xr*xr - dsr5 - rr5*xr*drc5[0];
      term6 = (bn[4]-dsc7*rr9)*xr*xr - bn[3] - rr7*xr*drc7[0];
      term7 = rr5*drc5[0] - (numtyp)2.0*bn[3]*xr + (dsc5+(numtyp)1.5*dsc7)*rr7*xr;
      numtyp tixx = ci*term3 + dix*term4 + dir*term5 +
        (numtyp)2.0*dsr5*qixx + (qiy*yr+qiz*zr)*dsc7*rr7 + (numtyp)2.0*qix*term7 + qir*term6;
      numtyp tkxx = ck*term3 - dkx*term4 - dkr*term5 +
        (numtyp)2.0*dsr5*qkxx + (qky*yr+qkz*zr)*dsc7*rr7 + (numtyp)2.0*qkx*term7 + qkr*term6;

      term3 = -dsr3 + term1*yr*yr - rr3*yr*drc3[1];
      term4 = rr3*drc3[1] - term1*yr - dsr5*yr;
      term5 = term2*yr*yr - dsr5 - rr5*yr*drc5[1];
      term6 = (bn[4]-dsc7*rr9)*yr*yr - bn[3] - rr7*yr*drc7[1];
      term7 = rr5*drc5[1] - (numtyp)2.0*bn[3]*yr + (dsc5+(numtyp)1.5*dsc7)*rr7*yr;
      numtyp tiyy = ci*term3 + diy*term4 + dir*term5 +
        (numtyp)2.0*dsr5*qiyy + (qix*xr+qiz*zr)*dsc7*rr7 + (numtyp)2.0*qiy*term7 + qir*term6;
      numtyp tkyy = ck*term3 - dky*term4 - dkr*term5 +
        (numtyp)2.0*dsr5*qkyy + (qkx*xr+qkz*zr)*dsc7*rr7 + (numtyp)2.0*qky*term7 + qkr*term6;

      term3 = -dsr3 + term1*zr*zr - rr3*zr*drc3[2];
      term4 = rr3*drc3[2] - term1*zr - dsr5*zr;
      term5 = term2*zr*zr - dsr5 - rr5*zr*drc5[2];
      term6 = (bn[4]-dsc7*rr9)*zr*zr - bn[3] - rr7*zr*drc7[2];
      term7 = rr5*drc5[2] - (numtyp)2.0*bn[3]*zr + (dsc5+(numtyp)1.5*dsc7)*rr7*zr;
      numtyp tizz = ci*term3 + diz*term4 + dir*term5 +
        (numtyp)2.0*dsr5*qizz + (qix*xr+qiy*yr)*dsc7*rr7 + (numtyp)2.0*qiz*term7 + qir*term6;
      numtyp tkzz = ck*term3 - dkz*term4 - dkr*term5 +
        (numtyp)2.0*dsr5*qkzz + (qkx*xr+qky*yr)*dsc7*rr7 + (numtyp)2.0*qkz*term7 + qkr*term6;

      term3 = term1*xr*yr - rr3*yr*drc3[0];
      term4 = rr3*drc3[0] - term1*xr;
      term5 = term2*xr*yr - rr5*yr*drc5[0];
      term6 = (bn[4]-dsc7*rr9)*xr*yr - rr7*yr*drc7[0];
      term7 = rr5*drc5[0] - term2*xr;
      numtyp tixy = ci*term3 - dsr5*dix*yr + diy*term4 + dir*term5 +
        (numtyp)2.0*dsr5*qixy - (numtyp)2.0*dsr7*yr*qix + (numtyp)2.0*qiy*term7 + qir*term6;
      numtyp tkxy = ck*term3 + dsr5*dkx*yr - dky*term4 - dkr*term5 +
        (numtyp)2.0*dsr5*qkxy - (numtyp)2.0*dsr7*yr*qkx +(numtyp) 2.0*qky*term7 + qkr*term6;

      term3 = term1*xr*zr - rr3*zr*drc3[0];
      term5 = term2*xr*zr - rr5*zr*drc5[0];
      term6 = (bn[4]-dsc7*rr9)*xr*zr - rr7*zr*drc7[0];
      numtyp tixz = ci*term3 - dsr5*dix*zr + diz*term4 + dir*term5 +
        (numtyp)2.0*dsr5*qixz - (numtyp)2.0*dsr7*zr*qix + (numtyp)2.0*qiz*term7 + qir*term6;
      numtyp tkxz = ck*term3 + dsr5*dkx*zr - dkz*term4 - dkr*term5 +
        (numtyp)2.0*dsr5*qkxz - (numtyp)2.0*dsr7*zr*qkx + (numtyp)2.0*qkz*term7 + qkr*term6;

      term3 = term1*yr*zr - rr3*zr*drc3[1];
      term4 = rr3*drc3[1] - term1*yr;
      term5 = term2*yr*zr - rr5*zr*drc5[1];
      term6 = (bn[4]-dsc7*rr9)*yr*zr - rr7*zr*drc7[1];
      term7 = rr5*drc5[1] - term2*yr;
      numtyp tiyz = ci*term3 - dsr5*diy*zr + diz*term4 + dir*term5 +
        (numtyp)2.0*dsr5*qiyz - (numtyp)2.0*dsr7*zr*qiy + (numtyp)2.0*qiz*term7 + qir*term6;
      numtyp tkyz = ck*term3 + dsr5*dky*zr - dkz*term4 - dkr*term5 +
        (numtyp)2.0*dsr5*qkyz - (numtyp)2.0*dsr7*zr*qky + (numtyp)2.0*qkz*term7 + qkr*term6;

      numtyp depx = tixx*ukxp + tixy*ukyp + tixz*ukzp - tkxx*uixp - tkxy*uiyp - tkxz*uizp;
      numtyp depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp - tkxy*uixp - tkyy*uiyp - tkyz*uizp;
      numtyp depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp - tkxz*uixp - tkyz*uiyp - tkzz*uizp;

      numtyp frcx = depx;
      numtyp frcy = depy;
      numtyp frcz = depz;

      // get the dEp/dR terms used for direct polarization force

      // tixx and tkxx
      term1 = bn[2] - psc3*rr5;
      term2 = bn[3] - psc5*rr7;
      term3 = -psr3 + term1*xr*xr - rr3*xr*prc3[0];
      term4 = rr3*prc3[0] - term1*xr - psr5*xr;
      term5 = term2*xr*xr - psr5 - rr5*xr*prc5[0];
      term6 = (bn[4]-psc7*rr9)*xr*xr - bn[3] - rr7*xr*prc7[0];
      term7 = rr5*prc5[0] - (numtyp)2.0*bn[3]*xr + (psc5+(numtyp)1.5*psc7)*rr7*xr;
      tixx = ci*term3 + dix*term4 + dir*term5 +
        (numtyp)2.0*psr5*qixx + (qiy*yr+qiz*zr)*psc7*rr7 + (numtyp)2.0*qix*term7 + qir*term6;
      tkxx = ck*term3 - dkx*term4 - dkr*term5 +
        (numtyp)2.0*psr5*qkxx + (qky*yr+qkz*zr)*psc7*rr7 + (numtyp)2.0*qkx*term7 + qkr*term6;

      // tiyy and tkyy
      term3 = -psr3 + term1*yr*yr - rr3*yr*prc3[1];
      term4 = rr3*prc3[1] - term1*yr - psr5*yr;
      term5 = term2*yr*yr - psr5 - rr5*yr*prc5[1];
      term6 = (bn[4]-psc7*rr9)*yr*yr - bn[3] - rr7*yr*prc7[1];
      term7 = rr5*prc5[1] - (numtyp)2.0*bn[3]*yr + (psc5+(numtyp)1.5*psc7)*rr7*yr;
      tiyy = ci*term3 + diy*term4 + dir*term5 +
        (numtyp)2.0*psr5*qiyy + (qix*xr+qiz*zr)*psc7*rr7 + (numtyp)2.0*qiy*term7 + qir*term6;
      tkyy = ck*term3 - dky*term4 - dkr*term5 +
        (numtyp)2.0*psr5*qkyy + (qkx*xr+qkz*zr)*psc7*rr7 + (numtyp)2.0*qky*term7 + qkr*term6;

      // tizz and tkzz
      term3 = -psr3 + term1*zr*zr - rr3*zr*prc3[2];
      term4 = rr3*prc3[2] - term1*zr - psr5*zr;
      term5 = term2*zr*zr - psr5 - rr5*zr*prc5[2];
      term6 = (bn[4]-psc7*rr9)*zr*zr - bn[3] - rr7*zr*prc7[2];
      term7 = rr5*prc5[2] - (numtyp)2.0*bn[3]*zr + (psc5+(numtyp)1.5*psc7)*rr7*zr;
      tizz = ci*term3 + diz*term4 + dir*term5 +
        (numtyp)2.0*psr5*qizz + (qix*xr+qiy*yr)*psc7*rr7 + (numtyp)2.0*qiz*term7 + qir*term6;
      tkzz = ck*term3 - dkz*term4 - dkr*term5 +
        (numtyp)2.0*psr5*qkzz + (qkx*xr+qky*yr)*psc7*rr7 + (numtyp)2.0*qkz*term7 + qkr*term6;

      // tixy and tkxy
      term3 = term1*xr*yr - rr3*yr*prc3[0];
      term4 = rr3*prc3[0] - term1*xr;
      term5 = term2*xr*yr - rr5*yr*prc5[0];
      term6 = (bn[4]-psc7*rr9)*xr*yr - rr7*yr*prc7[0];
      term7 = rr5*prc5[0] - term2*xr;
      tixy = ci*term3 - psr5*dix*yr + diy*term4 + dir*term5 +
        (numtyp)2.0*psr5*qixy - (numtyp)2.0*psr7*yr*qix + (numtyp)2.0*qiy*term7 + qir*term6;
      tkxy = ck*term3 + psr5*dkx*yr - dky*term4 - dkr*term5 +
        (numtyp)2.0*psr5*qkxy - (numtyp)2.0*psr7*yr*qkx + (numtyp)2.0*qky*term7 + qkr*term6;

      // tixz and tkxz
      term3 = term1*xr*zr - rr3*zr*prc3[0];
      term5 = term2*xr*zr - rr5*zr*prc5[0];
      term6 = (bn[4]-psc7*rr9)*xr*zr - rr7*zr*prc7[0];
      tixz = ci*term3 - psr5*dix*zr + diz*term4 + dir*term5 +
        (numtyp)2.0*psr5*qixz - (numtyp)2.0*psr7*zr*qix + (numtyp)2.0*qiz*term7 + qir*term6;
      tkxz = ck*term3 + psr5*dkx*zr - dkz*term4 - dkr*term5 +
        (numtyp)2.0*psr5*qkxz - (numtyp)2.0*psr7*zr*qkx + (numtyp)2.0*qkz*term7 + qkr*term6;

      // tiyz and tkyz
      term3 = term1*yr*zr - rr3*zr*prc3[1];
      term4 = rr3*prc3[1] - term1*yr;
      term5 = term2*yr*zr - rr5*zr*prc5[1];
      term6 = (bn[4]-psc7*rr9)*yr*zr - rr7*zr*prc7[1];
      term7 = rr5*prc5[1] - term2*yr;
      tiyz = ci*term3 - psr5*diy*zr + diz*term4 + dir*term5 +
        (numtyp)2.0*psr5*qiyz - (numtyp)2.0*psr7*zr*qiy + (numtyp)2.0*qiz*term7 + qir*term6;
      tkyz = ck*term3 + psr5*dky*zr - dkz*term4 - dkr*term5 +
        (numtyp)2.0*psr5*qkyz - (numtyp)2.0*psr7*zr*qky + (numtyp)2.0*qkz*term7 + qkr*term6;

      depx = tixx*ukx + tixy*uky + tixz*ukz - tkxx*uix - tkxy*uiy - tkxz*uiz;
      depy = tixy*ukx + tiyy*uky + tiyz*ukz - tkxy*uix - tkyy*uiy - tkyz*uiz;
      depz = tixz*ukx + tiyz*uky + tizz*ukz - tkxz*uix - tkyz*uiy - tkzz*uiz;

      frcx = frcx + depx;
      frcy = frcy + depy;
      frcz = frcz + depz;

      // get the dtau/dr terms used for mutual polarization force
      // poltyp == MUTUAL  && amoeba

      term1 = bn[2] - usc3*rr5;
      term2 = bn[3] - usc5*rr7;
      term3 = usr5 + term1;
      term4 = rr3 * factor_uscale;
      term5 = -xr*term3 + rc3[0]*term4;
      term6 = -usr5 + xr*xr*term2 - rr5*xr*urc5[0];
      tixx = uix*term5 + uir*term6;
      tkxx = ukx*term5 + ukr*term6;

      term5 = -yr*term3 + rc3[1]*term4;
      term6 = -usr5 + yr*yr*term2 - rr5*yr*urc5[1];
      tiyy = uiy*term5 + uir*term6;
      tkyy = uky*term5 + ukr*term6;

      term5 = -zr*term3 + rc3[2]*term4;
      term6 = -usr5 + zr*zr*term2 - rr5*zr*urc5[2];
      tizz = uiz*term5 + uir*term6;
      tkzz = ukz*term5 + ukr*term6;

      term4 = -usr5 * yr;
      term5 = -xr*term1 + rr3*urc3[0];
      term6 = xr*yr*term2 - rr5*yr*urc5[0];
      tixy = uix*term4 + uiy*term5 + uir*term6;
      tkxy = ukx*term4 + uky*term5 + ukr*term6;

      term4 = -usr5 * zr;
      term6 = xr*zr*term2 - rr5*zr*urc5[0];
      tixz = uix*term4 + uiz*term5 + uir*term6;
      tkxz = ukx*term4 + ukz*term5 + ukr*term6;

      term5 = -yr*term1 + rr3*urc3[1];
      term6 = yr*zr*term2 - rr5*zr*urc5[1];
      tiyz = uiy*term4 + uiz*term5 + uir*term6;
      tkyz = uky*term4 + ukz*term5 + ukr*term6;

      depx = tixx*ukxp + tixy*ukyp + tixz*ukzp
        + tkxx*uixp + tkxy*uiyp + tkxz*uizp;
      depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp
        + tkxy*uixp + tkyy*uiyp + tkyz*uizp;
      depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp
        + tkxz*uixp + tkyz*uiyp + tkzz*uizp;

      frcx = frcx + depx;
      frcy = frcy + depy;
      frcz = frcz + depz;

      f.x += frcx;
      f.y += frcy;
      f.z += frcz;

      if (EVFLAG && vflag) {
        numtyp vxx = xr * frcx;
        numtyp vxy = (numtyp)0.5 * (yr*frcx+xr*frcy);
        numtyp vxz = (numtyp)0.5 * (zr*frcx+xr*frcz);
        numtyp vyy = yr * frcy;
        numtyp vyz = (numtyp)0.5 * (zr*frcy+yr*frcz);
        numtyp vzz = zr * frcz;

        virial[0] -= vxx;
        virial[1] -= vyy;
        virial[2] -= vzz;
        virial[3] -= vxy;
        virial[4] -= vxz;
        virial[5] -= vyz;
      }
    } // nbor

  } // ii<inum

  // accumulate ufld and dufld to compute tep
  store_answers_tep(ufld,dufld,ii,inum,tid,t_per_atom,offset,i,tep);

  // accumate force, energy and virial
  store_answers_acc(f,energy,e_coul,virial,ii,inum,tid,t_per_atom,
     offset,eflag,vflag,ans,engv,NUM_BLOCKS_X);
}

/* ----------------------------------------------------------------------
   fphi_uind = induced potential from grid
   fphi_uind extracts the induced dipole potential from the particle mesh Ewald grid
------------------------------------------------------------------------- */

__kernel void k_amoeba_fphi_uind(const __global numtyp4 *restrict thetai1,
                          const __global numtyp4 *restrict thetai2,
                          const __global numtyp4 *restrict thetai3,
                          const __global int *restrict igrid,
                          const __global numtyp2 *restrict grid,
                          __global acctyp *restrict fdip_phi1,
                          __global acctyp *restrict fdip_phi2,
                          __global acctyp *restrict fdip_sum_phi,
                          const int bsorder, const int inum,
                          const int nzlo_out, const int nylo_out,
                          const int nxlo_out, const int ngridxy,
                          const int ngridx)
{
  int tid=THREAD_ID_X;
  int ii=tid+BLOCK_ID_X*BLOCK_SIZE_X;

  if (ii<inum) {

    const int nlpts = (bsorder-1) / 2;

    int istart = fast_mul(ii,4);
    const int igridx = igrid[istart];
    const int igridy = igrid[istart+1];
    const int igridz = igrid[istart+2];

    // now istart is used to index thetai1, thetai2 and thetai3
    istart = fast_mul(ii,bsorder);

    // extract the permanent multipole field at each site

    numtyp tuv100_1 = (numtyp)0.0;
    numtyp tuv010_1 = (numtyp)0.0;
    numtyp tuv001_1 = (numtyp)0.0;
    numtyp tuv200_1 = (numtyp)0.0;
    numtyp tuv020_1 = (numtyp)0.0;
    numtyp tuv002_1 = (numtyp)0.0;
    numtyp tuv110_1 = (numtyp)0.0;
    numtyp tuv101_1 = (numtyp)0.0;
    numtyp tuv011_1 = (numtyp)0.0;
    numtyp tuv100_2 = (numtyp)0.0;
    numtyp tuv010_2 = (numtyp)0.0;
    numtyp tuv001_2 = (numtyp)0.0;
    numtyp tuv200_2 = (numtyp)0.0;
    numtyp tuv020_2 = (numtyp)0.0;
    numtyp tuv002_2 = (numtyp)0.0;
    numtyp tuv110_2 = (numtyp)0.0;
    numtyp tuv101_2 = (numtyp)0.0;
    numtyp tuv011_2 = (numtyp)0.0;
    numtyp tuv000 = (numtyp)0.0;
    numtyp tuv001 = (numtyp)0.0;
    numtyp tuv010 = (numtyp)0.0;
    numtyp tuv100 = (numtyp)0.0;
    numtyp tuv200 = (numtyp)0.0;
    numtyp tuv020 = (numtyp)0.0;
    numtyp tuv002 = (numtyp)0.0;
    numtyp tuv110 = (numtyp)0.0;
    numtyp tuv101 = (numtyp)0.0;
    numtyp tuv011 = (numtyp)0.0;
    numtyp tuv300 = (numtyp)0.0;
    numtyp tuv030 = (numtyp)0.0;
    numtyp tuv003 = (numtyp)0.0;
    numtyp tuv210 = (numtyp)0.0;
    numtyp tuv201 = (numtyp)0.0;
    numtyp tuv120 = (numtyp)0.0;
    numtyp tuv021 = (numtyp)0.0;
    numtyp tuv102 = (numtyp)0.0;
    numtyp tuv012 = (numtyp)0.0;
    numtyp tuv111 = (numtyp)0.0;

    int k = (igridz - nzlo_out) - nlpts;
    for (int kb = 0; kb < bsorder; kb++) {
      const int mz = fast_mul(k, ngridxy);
      const int i3 = istart + kb;
      const numtyp4 tha3 = thetai3[i3];
      const numtyp v0 = tha3.x; // thetai3[m][kb][0];
      const numtyp v1 = tha3.y; // thetai3[m][kb][1];
      const numtyp v2 = tha3.z; // thetai3[m][kb][2];
      const numtyp v3 = tha3.w; // thetai3[m][kb][3];
      numtyp tu00_1 = (numtyp)0.0;
      numtyp tu01_1 = (numtyp)0.0;
      numtyp tu10_1 = (numtyp)0.0;
      numtyp tu20_1 = (numtyp)0.0;
      numtyp tu11_1 = (numtyp)0.0;
      numtyp tu02_1 = (numtyp)0.0;
      numtyp tu00_2 = (numtyp)0.0;
      numtyp tu01_2 = (numtyp)0.0;
      numtyp tu10_2 = (numtyp)0.0;
      numtyp tu20_2 = (numtyp)0.0;
      numtyp tu11_2 = (numtyp)0.0;
      numtyp tu02_2 = (numtyp)0.0;
      numtyp tu00 = (numtyp)0.0;
      numtyp tu10 = (numtyp)0.0;
      numtyp tu01 = (numtyp)0.0;
      numtyp tu20 = (numtyp)0.0;
      numtyp tu11 = (numtyp)0.0;
      numtyp tu02 = (numtyp)0.0;
      numtyp tu30 = (numtyp)0.0;
      numtyp tu21 = (numtyp)0.0;
      numtyp tu12 = (numtyp)0.0;
      numtyp tu03 = (numtyp)0.0;

      int j = (igridy - nylo_out) - nlpts;
      for (int jb = 0; jb < bsorder; jb++) {
        const int my = mz + fast_mul(j, ngridx);
        const int i2 = istart + jb;
        const numtyp4 tha2 = thetai2[i2];
        const numtyp u0 = tha2.x; // thetai2[m][jb][0];
        const numtyp u1 = tha2.y; // thetai2[m][jb][1];
        const numtyp u2 = tha2.z; // thetai2[m][jb][2];
        const numtyp u3 = tha2.w; // thetai2[m][jb][3];
        numtyp t0_1 = (numtyp)0.0;
        numtyp t1_1 = (numtyp)0.0;
        numtyp t2_1 = (numtyp)0.0;
        numtyp t0_2 = (numtyp)0.0;
        numtyp t1_2 = (numtyp)0.0;
        numtyp t2_2 = (numtyp)0.0;
        numtyp t3 = (numtyp)0.0;

        int i = (igridx - nxlo_out) - nlpts;
        for (int ib = 0; ib < bsorder; ib++) {
          const int i1 = istart + ib;
          const numtyp4 tha1 = thetai1[i1];
          const int gidx = my + i; // k*ngridxy + j*ngridx + i;
          const numtyp2 tq = grid[gidx];
          const numtyp tq_1 = tq.x; //grid[gidx];
          const numtyp tq_2 = tq.y; //grid[gidx+1];
          t0_1 += tq_1*tha1.x;
          t1_1 += tq_1*tha1.y;
          t2_1 += tq_1*tha1.z;
          t0_2 += tq_2*tha1.x;
          t1_2 += tq_2*tha1.y;
          t2_2 += tq_2*tha1.z;
          t3 += (tq_1+tq_2)*tha1.w;
          i++;
        }

        tu00_1 += t0_1*u0;
        tu10_1 += t1_1*u0;
        tu01_1 += t0_1*u1;
        tu20_1 += t2_1*u0;
        tu11_1 += t1_1*u1;
        tu02_1 += t0_1*u2;
        tu00_2 += t0_2*u0;
        tu10_2 += t1_2*u0;
        tu01_2 += t0_2*u1;
        tu20_2 += t2_2*u0;
        tu11_2 += t1_2*u1;
        tu02_2 += t0_2*u2;
        numtyp t0 = t0_1 + t0_2;
        numtyp t1 = t1_1 + t1_2;
        numtyp t2 = t2_1 + t2_2;
        tu00 += t0*u0;
        tu10 += t1*u0;
        tu01 += t0*u1;
        tu20 += t2*u0;
        tu11 += t1*u1;
        tu02 += t0*u2;
        tu30 += t3*u0;
        tu21 += t2*u1;
        tu12 += t1*u2;
        tu03 += t0*u3;
        j++;
      }

      tuv100_1 += tu10_1*v0;
      tuv010_1 += tu01_1*v0;
      tuv001_1 += tu00_1*v1;
      tuv200_1 += tu20_1*v0;
      tuv020_1 += tu02_1*v0;
      tuv002_1 += tu00_1*v2;
      tuv110_1 += tu11_1*v0;
      tuv101_1 += tu10_1*v1;
      tuv011_1 += tu01_1*v1;
      tuv100_2 += tu10_2*v0;
      tuv010_2 += tu01_2*v0;
      tuv001_2 += tu00_2*v1;
      tuv200_2 += tu20_2*v0;
      tuv020_2 += tu02_2*v0;
      tuv002_2 += tu00_2*v2;
      tuv110_2 += tu11_2*v0;
      tuv101_2 += tu10_2*v1;
      tuv011_2 += tu01_2*v1;
      tuv000 += tu00*v0;
      tuv100 += tu10*v0;
      tuv010 += tu01*v0;
      tuv001 += tu00*v1;
      tuv200 += tu20*v0;
      tuv020 += tu02*v0;
      tuv002 += tu00*v2;
      tuv110 += tu11*v0;
      tuv101 += tu10*v1;
      tuv011 += tu01*v1;
      tuv300 += tu30*v0;
      tuv030 += tu03*v0;
      tuv003 += tu00*v3;
      tuv210 += tu21*v0;
      tuv201 += tu20*v1;
      tuv120 += tu12*v0;
      tuv021 += tu02*v1;
      tuv102 += tu10*v2;
      tuv012 += tu01*v2;
      tuv111 += tu11*v1;
      k++;
    }

    int idx;
    acctyp fdip_buf[20];

    fdip_buf[0] = (numtyp)0.0;
    fdip_buf[1] = tuv100_1;
    fdip_buf[2] = tuv010_1;
    fdip_buf[3] = tuv001_1;
    fdip_buf[4] = tuv200_1;
    fdip_buf[5] = tuv020_1;
    fdip_buf[6] = tuv002_1;
    fdip_buf[7] = tuv110_1;
    fdip_buf[8] = tuv101_1;
    fdip_buf[9] = tuv011_1;
    idx = ii;
    for (int m = 0; m < 10; m++) {
      fdip_phi1[idx] = fdip_buf[m];
      idx += inum;
    }

    fdip_buf[0] = (numtyp)0.0;
    fdip_buf[1] = tuv100_2;
    fdip_buf[2] = tuv010_2;
    fdip_buf[3] = tuv001_2;
    fdip_buf[4] = tuv200_2;
    fdip_buf[5] = tuv020_2;
    fdip_buf[6] = tuv002_2;
    fdip_buf[7] = tuv110_2;
    fdip_buf[8] = tuv101_2;
    fdip_buf[9] = tuv011_2;
    idx = ii;
    for (int m = 0; m < 10; m++) {
      fdip_phi2[idx] = fdip_buf[m];
      idx += inum;
    }

    fdip_buf[0] = tuv000;
    fdip_buf[1] = tuv100;
    fdip_buf[2] = tuv010;
    fdip_buf[3] = tuv001;
    fdip_buf[4] = tuv200;
    fdip_buf[5] = tuv020;
    fdip_buf[6] = tuv002;
    fdip_buf[7] = tuv110;
    fdip_buf[8] = tuv101;
    fdip_buf[9] = tuv011;
    fdip_buf[10] = tuv300;
    fdip_buf[11] = tuv030;
    fdip_buf[12] = tuv003;
    fdip_buf[13] = tuv210;
    fdip_buf[14] = tuv201;
    fdip_buf[15] = tuv120;
    fdip_buf[16] = tuv021;
    fdip_buf[17] = tuv102;
    fdip_buf[18] = tuv012;
    fdip_buf[19] = tuv111;
    idx = ii;
    for (int m = 0; m < 20; m++) {
      fdip_sum_phi[idx] = fdip_buf[m];
      idx += inum;
    }
  }
}

/* ----------------------------------------------------------------------
   fphi_mpole = multipole potential from grid
   fphi_mpole extracts the permanent multipole potential from
   the particle mesh Ewald grid
------------------------------------------------------------------------- */

__kernel void k_amoeba_fphi_mpole(const __global numtyp4 *restrict thetai1,
                          const __global numtyp4 *restrict thetai2,
                          const __global numtyp4 *restrict thetai3,
                          const __global int *restrict igrid,
                          const __global numtyp2 *restrict grid,
                          __global acctyp *restrict fphi,
                          const int bsorder, const int inum, const numtyp felec,
                          const int nzlo_out, const int nylo_out,
                          const int nxlo_out, const int ngridxy,
                          const int ngridx)
{
  int tid=THREAD_ID_X;
  int ii=tid+BLOCK_ID_X*BLOCK_SIZE_X;

  if (ii<inum) {

    int nlpts = (bsorder-1) / 2;

    int istart = fast_mul(ii,4);
    int igridx = igrid[istart];
    int igridy = igrid[istart+1];
    int igridz = igrid[istart+2];

    // now istart is used to index thetai1, thetai2 and thetai3
    istart = fast_mul(ii,bsorder);

    // extract the permanent multipole field at each site

    numtyp tuv000 = (numtyp)0.0;
    numtyp tuv001 = (numtyp)0.0;
    numtyp tuv010 = (numtyp)0.0;
    numtyp tuv100 = (numtyp)0.0;
    numtyp tuv200 = (numtyp)0.0;
    numtyp tuv020 = (numtyp)0.0;
    numtyp tuv002 = (numtyp)0.0;
    numtyp tuv110 = (numtyp)0.0;
    numtyp tuv101 = (numtyp)0.0;
    numtyp tuv011 = (numtyp)0.0;
    numtyp tuv300 = (numtyp)0.0;
    numtyp tuv030 = (numtyp)0.0;
    numtyp tuv003 = (numtyp)0.0;
    numtyp tuv210 = (numtyp)0.0;
    numtyp tuv201 = (numtyp)0.0;
    numtyp tuv120 = (numtyp)0.0;
    numtyp tuv021 = (numtyp)0.0;
    numtyp tuv102 = (numtyp)0.0;
    numtyp tuv012 = (numtyp)0.0;
    numtyp tuv111 = (numtyp)0.0;

    int k = (igridz - nzlo_out) - nlpts;
    for (int kb = 0; kb < bsorder; kb++) {
      int i3 = istart + kb;
      numtyp4 tha3 = thetai3[i3];
      numtyp v0 = tha3.x;
      numtyp v1 = tha3.y;
      numtyp v2 = tha3.z;
      numtyp v3 = tha3.w;
      numtyp tu00 = (numtyp)0.0;
      numtyp tu10 = (numtyp)0.0;
      numtyp tu01 = (numtyp)0.0;
      numtyp tu20 = (numtyp)0.0;
      numtyp tu11 = (numtyp)0.0;
      numtyp tu02 = (numtyp)0.0;
      numtyp tu30 = (numtyp)0.0;
      numtyp tu21 = (numtyp)0.0;
      numtyp tu12 = (numtyp)0.0;
      numtyp tu03 = (numtyp)0.0;

      int j = (igridy - nylo_out) - nlpts;
      for (int jb = 0; jb < bsorder; jb++) {
        int i2 = istart + jb;
        numtyp4 tha2 = thetai2[i2];
        numtyp u0 = tha2.x;
        numtyp u1 = tha2.y;
        numtyp u2 = tha2.z;
        numtyp u3 = tha2.w;
        numtyp t0 = (numtyp)0.0;
        numtyp t1 = (numtyp)0.0;
        numtyp t2 = (numtyp)0.0;
        numtyp t3 = (numtyp)0.0;

        int i = (igridx - nxlo_out) - nlpts;
        for (int ib = 0; ib < bsorder; ib++) {
          int i1 = istart + ib;
          numtyp4 tha1 = thetai1[i1];
          int gidx = k*ngridxy + j*ngridx + i;
          numtyp tq = grid[gidx].x;
          t0 += tq*tha1.x;
          t1 += tq*tha1.y;
          t2 += tq*tha1.z;
          t3 += tq*tha1.w;
          i++;
        }

        tu00 += t0*u0;
        tu10 += t1*u0;
        tu01 += t0*u1;
        tu20 += t2*u0;
        tu11 += t1*u1;
        tu02 += t0*u2;
        tu30 += t3*u0;
        tu21 += t2*u1;
        tu12 += t1*u2;
        tu03 += t0*u3;
        j++;
      }

      tuv000 += tu00*v0;
      tuv100 += tu10*v0;
      tuv010 += tu01*v0;
      tuv001 += tu00*v1;
      tuv200 += tu20*v0;
      tuv020 += tu02*v0;
      tuv002 += tu00*v2;
      tuv110 += tu11*v0;
      tuv101 += tu10*v1;
      tuv011 += tu01*v1;
      tuv300 += tu30*v0;
      tuv030 += tu03*v0;
      tuv003 += tu00*v3;
      tuv210 += tu21*v0;
      tuv201 += tu20*v1;
      tuv120 += tu12*v0;
      tuv021 += tu02*v1;
      tuv102 += tu10*v2;
      tuv012 += tu01*v2;
      tuv111 += tu11*v1;
      k++;
    }

    numtyp buf[20];
    buf[0] = tuv000;
    buf[1] = tuv100;
    buf[2] = tuv010;
    buf[3] = tuv001;
    buf[4] = tuv200;
    buf[5] = tuv020;
    buf[6] = tuv002;
    buf[7] = tuv110;
    buf[8] = tuv101;
    buf[9] = tuv011;
    buf[10] = tuv300;
    buf[11] = tuv030;
    buf[12] = tuv003;
    buf[13] = tuv210;
    buf[14] = tuv201;
    buf[15] = tuv120;
    buf[16] = tuv021;
    buf[17] = tuv102;
    buf[18] = tuv012;
    buf[19] = tuv111;

    int idx = ii;
    for (int m = 0; m < 20; m++) {
      fphi[idx] = felec * buf[m];
      idx += inum;
    }
  }
}

/* ----------------------------------------------------------------------
   scan standard neighbor list and make it compatible with 1-5 neighbors
   if IJ entry is a 1-2,1-3,1-4 neighbor then adjust offset to SBBITS15
   else scan special15 to see if a 1-5 neighbor and adjust offset to SBBITS15
   else do nothing to IJ entry
------------------------------------------------------------------------- */

__kernel void k_amoeba_special15(__global int * dev_nbor,
                          const __global int * dev_packed,
                          const __global tagint *restrict tag,
                          const __global int *restrict nspecial15,
                          const __global tagint *restrict special15,
                          const int inum, const int nall, const int nbor_pitch,
                          const int t_per_atom) {
  int tid, ii, offset, n_stride, i;
  atom_info(t_per_atom,ii,tid,offset);

  if (ii<inum) {

    int numj, nbor, nbor_end;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);

    int n15 = nspecial15[ii];

    for ( ; nbor<nbor_end; nbor+=n_stride) {

      int sj=dev_packed[nbor];
      int which = sj >> SBBITS & 3;
      int j = sj & NEIGHMASK;
      tagint jtag = tag[j];

      if (!which) {
        int offset=ii;
        for (int k=0; k<n15; k++) {
          if (special15[offset] == jtag) {
            which = 4;
            break;
          }
          offset += nall;
        }
      }

      if (which) dev_nbor[nbor] = j ^ (which << SBBITS15);
    } // for nbor

  } // if ii
}

__kernel void k_amoeba_short_nbor(const __global numtyp4 *restrict x_,
                                  const __global int * dev_nbor,
                                  const __global int * dev_packed,
                                  __global int * dev_short_nbor,
                                  const numtyp off2,
                                  const int inum, const int nbor_pitch,
                                  const int t_per_atom) {
  __local int n_stride;
  int tid, ii, offset;
  atom_info(t_per_atom,ii,tid,offset);

  if (ii<inum) {
    int nbor, nbor_end;
    int i, numj;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,i,numj,
              n_stride,nbor_end,nbor);

    numtyp4 ix; fetch4(ix,i,pos_tex); //x_[i];

    int ncount = 0;
    int m = nbor;
    dev_short_nbor[m] = 0;
    int nbor_short = nbor+n_stride;

    for ( ; nbor<nbor_end; nbor+=n_stride) {

      int j=dev_packed[nbor];
      int nj = j;
      j &= NEIGHMASK15;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];

      // Compute r12
      numtyp delx = ix.x-jx.x;
      numtyp dely = ix.y-jx.y;
      numtyp delz = ix.z-jx.z;
      numtyp rsq = delx*delx+dely*dely+delz*delz;

      if (rsq<off2) {
        dev_short_nbor[nbor_short] = nj;
        nbor_short += n_stride;
        ncount++;
      }
    } // for nbor

    // store the number of neighbors for each thread
    dev_short_nbor[m] = ncount;

  } // if ii
}
