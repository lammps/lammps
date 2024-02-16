// **************************************************************************
//                                   hippo.cu
//                             -------------------
//                          Trung Dac Nguyen (Northwestern)
//
//  Device code for acceleration of the hippo pair style
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                :
//    email                : trung.nguyen@northwestern.edu
// ***************************************************************************

#if defined(NV_KERNEL) || defined(USE_HIP)

#include "lal_hippo_extra.h"
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

#define store_answers_hippo_tq(tq, ii, inum,tid, t_per_atom, offset, i,    \
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

#define store_answers_hippo_tq(tq, ii, inum,tid, t_per_atom, offset, i,    \
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
   repulsion = Pauli repulsion interactions
   adapted from Tinker erepel1b() routine
------------------------------------------------------------------------- */

__kernel void k_hippo_repulsion(const __global numtyp4 *restrict x_,
                                const __global numtyp4 *restrict extra,
                                const __global numtyp4 *restrict coeff_rep,
                                const __global numtyp4 *restrict sp_nonpolar,
                                const __global int *dev_nbor,
                                const __global int *dev_packed,
                                const __global int *dev_short_nbor,
                                __global acctyp3 *restrict ans,
                                __global acctyp *restrict engv,
                                __global acctyp3 *restrict tep,
                                const int eflag, const int vflag, const int inum,
                                const int nall, const int nbor_pitch,
                                const int t_per_atom, const numtyp aewald,
                                const numtyp off2, const numtyp cut2,
                                const numtyp c0, const numtyp c1, const numtyp c2,
                                const numtyp c3, const numtyp c4, const numtyp c5)
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
    //numtyp ci  = pol1i.x;  // rpole[i][0];
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
    int itype = pol3i.z; // amtype[i];
    numtyp sizi = coeff_rep[itype].x; // sizpr[itype];
    numtyp dmpi = coeff_rep[itype].y; // dmppr[itype];
    numtyp vali = coeff_rep[itype].z; // elepr[itype];

    for ( ; nbor<nbor_end; nbor+=n_stride) {

      int jextra=nbor_mem[nbor];
      int j = jextra & NEIGHMASK15;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];

      // Compute r12
      numtyp xr = jx.x - ix.x;
      numtyp yr = jx.y - ix.y;
      numtyp zr = jx.z - ix.z;
      numtyp r2 = xr*xr + yr*yr + zr*zr;

      const numtyp4 pol1j = polar1[j];
      //numtyp ck  = pol1j.x;  // rpole[j][0];
      numtyp dkx = pol1j.y;    // rpole[j][1];
      numtyp dky = pol1j.z;    // rpole[j][2];
      numtyp dkz = pol1j.w;    // rpole[j][3];
      const numtyp4 pol2j = polar2[j];
      numtyp qkxx = pol2j.x;   // rpole[j][4];
      numtyp qkxy = pol2j.y;   // rpole[j][5];
      numtyp qkxz = pol2j.z;   // rpole[j][6];
      numtyp qkyy = pol2j.w;   // rpole[j][8];
      const numtyp4 pol3j = polar3[j];
      numtyp qkyz = pol3j.x;   // rpole[j][9];
      numtyp qkzz = pol3j.y;   // rpole[j][12];
      int jtype = pol3j.z;     // amtype[j];

      numtyp sizk = coeff_rep[jtype].x; // sizpr[jtype];
      numtyp dmpk = coeff_rep[jtype].y; // dmppr[jtype];
      numtyp valk = coeff_rep[jtype].z; // elepr[jtype];

      const numtyp4 sp_nonpol = sp_nonpolar[sbmask15(jextra)];
      numtyp factor_repel = sp_nonpol.x; // factor_repel = special_repel[sbmask15(j)];
      if (factor_repel == (numtyp)0) continue;

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

      numtyp r = ucl_sqrt(r2);
      numtyp rinv = ucl_recip(r);
      numtyp r2inv = rinv*rinv;
      numtyp rr1 = rinv;
      numtyp rr3 = rr1 * r2inv;
      numtyp rr5 = (numtyp)3.0 * rr3 * r2inv;
      numtyp rr7 = (numtyp)5.0 * rr5 * r2inv;
      numtyp rr9 = (numtyp)7.0 * rr7 * r2inv;
      numtyp rr11 = (numtyp)9.0 * rr9 * r2inv;

      // get damping coefficients for the Pauli repulsion energy
      numtyp dmpik[11];
      damprep(r,r2,rr1,rr3,rr5,rr7,rr9,rr11,11,dmpi,dmpk,dmpik);

      // calculate intermediate terms needed for the energy

      numtyp term1 = vali*valk;
      numtyp term2 = valk*dir - vali*dkr + dik;
      numtyp term3 = vali*qkr + valk*qir - dir*dkr + (numtyp)2.0*(dkqi-diqk+qiqk);
      numtyp term4 = dir*qkr - dkr*qir - (numtyp)4.0*qik;
      numtyp term5 = qir*qkr;
      numtyp eterm = term1*dmpik[0] + term2*dmpik[2] +
        term3*dmpik[4] + term4*dmpik[6] + term5*dmpik[8];

      // compute the Pauli repulsion energy for this interaction

      numtyp sizik = sizi * sizk * factor_repel;
      numtyp e = sizik * eterm * rr1;

      // calculate intermediate terms for force and torque

      numtyp de = term1*dmpik[2] + term2*dmpik[4] + term3*dmpik[6] +
        term4*dmpik[8] + term5*dmpik[10];
      term1 = -valk*dmpik[2] + dkr*dmpik[4] - qkr*dmpik[6];
      term2 = vali*dmpik[2] + dir*dmpik[4] + qir*dmpik[6];
      term3 = (numtyp)2.0 * dmpik[4];
      term4 = (numtyp)2.0 * (-valk*dmpik[4] + dkr*dmpik[6] - qkr*dmpik[8]);
      term5 = (numtyp)2.0 * (-vali*dmpik[4] - dir*dmpik[6] - qir*dmpik[8]);
      numtyp term6 = (numtyp)4.0 * dmpik[6];

      // compute the force components for this interaction

      numtyp frcx = de*xr + term1*dix + term2*dkx + term3*(diqkx-dkqix) +
        term4*qix + term5*qkx + term6*(qixk+qkxi);
      numtyp frcy = de*yr + term1*diy + term2*dky + term3*(diqky-dkqiy) +
        term4*qiy + term5*qky + term6*(qiyk+qkyi);
      numtyp frcz = de*zr + term1*diz + term2*dkz + term3*(diqkz-dkqiz) +
        term4*qiz + term5*qkz + term6*(qizk+qkzi);

      frcx = frcx*rr1 + eterm*rr3*xr;
      frcy = frcy*rr1 + eterm*rr3*yr;
      frcz = frcz*rr1 + eterm*rr3*zr;
      frcx = sizik * frcx;
      frcy = sizik * frcy;
      frcz = sizik * frcz;

      // compute the torque components for this interaction

      numtyp ttmix = -dmpik[2]*dikx + term1*dirx + term3*(dqikx+dkqirx) -
        term4*qirx - term6*(qikrx+qikx);
      numtyp ttmiy = -dmpik[2]*diky + term1*diry + term3*(dqiky+dkqiry) -
        term4*qiry - term6*(qikry+qiky);
      numtyp ttmiz = -dmpik[2]*dikz + term1*dirz + term3*(dqikz+dkqirz) -
        term4*qirz - term6*(qikrz+qikz);
      ttmix = sizik * ttmix * rr1;
      ttmiy = sizik * ttmiy * rr1;
      ttmiz = sizik * ttmiz * rr1;

      // use energy switching if near the cutoff distance

      if (r2 > cut2) {
        numtyp r3 = r2 * r;
        numtyp r4 = r2 * r2;
        numtyp r5 = r2 * r3;
        numtyp taper = c5*r5 + c4*r4 + c3*r3 + c2*r2 + c1*r + c0;
        numtyp dtaper = (numtyp)5.0*c5*r4 + (numtyp)4.0*c4*r3 +
          (numtyp)3.0*c3*r2 + (numtyp)2.0*c2*r + c1;
        dtaper *= e * rr1;
        e *= taper;
        frcx = frcx*taper - dtaper*xr;
        frcy = frcy*taper - dtaper*yr;
        frcz = frcz*taper - dtaper*zr;
        ttmix *= taper;
        ttmiy *= taper;
        ttmiz *= taper;
      }

      energy += e;

      // increment force-based gradient and torque on atom I

      f.x -= frcx;
      f.y -= frcy;
      f.z -= frcz;
      tq.x += ttmix;
      tq.y += ttmiy;
      tq.z += ttmiz;

      // increment the internal virial tensor components
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
  store_answers_hippo_tq(tq,ii,inum,tid,t_per_atom,offset,i,tep);
  // accumate force, energy and virial
  store_answers_q(f,energy,e_coul,virial,ii,inum,tid,t_per_atom,
     offset,eflag,vflag,ans,engv);
}

/* ----------------------------------------------------------------------
   dispersion = real-space portion of Ewald dispersion
   adapted from Tinker edreal1d() routine
------------------------------------------------------------------------- */

__kernel void k_hippo_dispersion(const __global numtyp4 *restrict x_,
                                 const __global numtyp4 *restrict extra,
                                 const __global numtyp4 *restrict coeff_amtype,
                                 const __global numtyp4 *restrict coeff_amclass,
                                 const __global numtyp4 *restrict sp_nonpolar,
                                 const __global int *dev_nbor,
                                 const __global int *dev_packed,
                                 const __global int *dev_short_nbor,
                                 __global acctyp3 *restrict ans,
                                 __global acctyp *restrict engv,
                                 const int eflag, const int vflag, const int inum,
                                 const int nall, const int nbor_pitch,
                                 const int t_per_atom, const numtyp aewald,
                                 const numtyp off2)
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

  const __global numtyp4* polar3 = &extra[2*nall];

  if (ii<inum) {
    int itype,iclass;
    numtyp ci,ai;

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

    itype  = polar3[i].z;            // amtype[i];
    iclass = coeff_amtype[itype].w;  // amtype2class[itype];
    ci = coeff_amclass[iclass].x;    // csix[iclass];
    ai = coeff_amclass[iclass].y;    // adisp[iclass];

    for ( ; nbor<nbor_end; nbor+=n_stride) {

      int jextra=nbor_mem[nbor];
      int j = jextra & NEIGHMASK15;

      numtyp4 jx; fetch4(jx,j,pos_tex); //x_[j];

      // Compute r12
      numtyp xr = ix.x - jx.x;
      numtyp yr = ix.y - jx.y;
      numtyp zr = ix.z - jx.z;
      numtyp r2 = xr*xr + yr*yr + zr*zr;

      int jtype =   polar3[j].z; // amtype[j];
      int jclass = coeff_amtype[jtype].w;  // amtype2class[jtype];
      numtyp ck = coeff_amclass[jclass].x;    // csix[jclass];
      numtyp ak = coeff_amclass[jclass].y;    // adisp[jclass];

      numtyp r6 = r2*r2*r2;
      numtyp ralpha2 = r2 * aewald*aewald;
      numtyp term = (numtyp)1.0 + ralpha2 + (numtyp)0.5*ralpha2*ralpha2;
      numtyp expterm = ucl_exp(-ralpha2);
      numtyp expa = expterm * term;

      // find the damping factor for the dispersion interaction

      numtyp r = ucl_sqrt(r2);
      numtyp r7 = r6 * r;
      numtyp di = ai * r;
      numtyp di2 = di * di;
      numtyp di3 = di * di2;
      numtyp dk = ak * r;
      numtyp expi = ucl_exp(-di);
      numtyp expk = ucl_exp(-dk);

      numtyp ai2,ak2;
      numtyp di4,di5;
      numtyp dk2,dk3;
      numtyp ti,ti2;
      numtyp tk,tk2;
      numtyp damp3,damp5;
      numtyp ddamp;
      const numtyp4 sp_nonpol = sp_nonpolar[sbmask15(jextra)];
      numtyp factor_disp = sp_nonpol.y; // factor_disp = special_disp[sbmask15(j)];

      if (ai != ak) {
        ai2 = ai * ai;
        ak2 = ak * ak;
        dk2 = dk * dk;
        dk3 = dk * dk2;
        ti = ak2 / (ak2-ai2);
        ti2 = ti * ti;
        tk = ai2 / (ai2-ak2);
        tk2 = tk * tk;
        damp3 = (numtyp)1.0 - ti2*((numtyp)1.0 + di + (numtyp)0.5*di2) * expi
          - tk2*((numtyp)1.0 + dk + (numtyp)0.5*dk2) * expk
          - (numtyp)2.0*ti2*tk * ((numtyp)1.0 + di)* expi
          - (numtyp)2.0*tk2*ti * ((numtyp)1.0 + dk) *expk;
        damp5 = (numtyp)1.0 - ti2*((numtyp)1.0 + di + (numtyp)0.5*di2 + di3/(numtyp)6.0) * expi
          - tk2*((numtyp)1.0 + dk + (numtyp)0.5*dk2 + dk3/(numtyp)6.0) * expk
          - (numtyp)2.0*ti2*tk*((numtyp)1.0 + di + di2/(numtyp)3.0) * expi
          - (numtyp)2.0*tk2*ti*((numtyp)1.0 + dk + dk2/(numtyp)3.0) * expk;
        ddamp = (numtyp)0.25 * di2 * ti2 * ai * expi * (r*ai+(numtyp)4.0*tk - (numtyp)1.0) +
          (numtyp)0.25 * dk2 * tk2 * ak * expk * (r*ak+(numtyp)4.0*ti-(numtyp)1.0);

      } else {
        di4 = di2 * di2;
        di5 = di2 * di3;
        damp3 = (numtyp)1.0 - ((numtyp)1.0+di+(numtyp)0.5*di2 + (numtyp)7.0*di3/(numtyp)48.0+di4/(numtyp)48.0)*expi;
        damp5 = (numtyp)1.0 - ((numtyp)1.0+di+(numtyp)0.5*di2 + di3/(numtyp)6.0+di4/(numtyp)24.0+di5/(numtyp)144.0)*expi;
        ddamp = ai * expi * (di5-(numtyp)3.0*di3-(numtyp)3.0*di2) / (numtyp)96.0;
      }

      numtyp damp = (numtyp)1.5*damp5 - (numtyp)0.5*damp3;

      // apply damping and scaling factors for this interaction

      numtyp scale = factor_disp * damp*damp;
      scale = scale - (numtyp)1.0;
      numtyp e = -ci * ck * (expa+scale) / r6;
      numtyp rterm = -ralpha2*ralpha2*ralpha2 * expterm / r;
      numtyp de = (numtyp)-6.0*e/r2 - ci*ck*rterm/r7 - (numtyp)2.0*ci*ck*factor_disp*damp*ddamp/r7;

      energy+= e;

      // increment the damped dispersion derivative components

      numtyp dedx = de * xr;
      numtyp dedy = de * yr;
      numtyp dedz = de * zr;
      f.x -= dedx;
      f.y -= dedy;
      f.z -= dedz;

      // increment the internal virial tensor components

      numtyp vxx = xr * dedx;
      numtyp vyx = yr * dedx;
      numtyp vzx = zr * dedx;
      numtyp vyy = yr * dedy;
      numtyp vzy = zr * dedy;
      numtyp vzz = zr * dedz;

      virial[0] -= vxx;
      virial[1] -= vyy;
      virial[2] -= vzz;
      virial[3] -= vyx;
      virial[4] -= vzx;
      virial[5] -= vzy;
    } // nbor

  } // ii<inum

  store_answers_acc(f,energy,e_coul,virial,ii,inum,tid,t_per_atom,
     offset,eflag,vflag,ans,engv,NUM_BLOCKS_X);
}

/* ----------------------------------------------------------------------
   multipole_real = real-space portion of multipole
   adapted from Tinker emreal1d() routine
------------------------------------------------------------------------- */

__kernel void k_hippo_multipole(const __global numtyp4 *restrict x_,
                                const __global numtyp4 *restrict extra,
                                const __global numtyp4 *restrict coeff_amtype,
                                const __global numtyp4 *restrict coeff_amclass,
                                const __global numtyp4 *restrict sp_polar,
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
  const __global numtyp4* polar6 = &extra[5*nall];

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
    int itype  = pol3i.z;        // amtype[i];
    int iclass = coeff_amtype[itype].w;  // amtype2class[itype];
    numtyp corei = coeff_amclass[iclass].z;  // pcore[iclass];
    numtyp alphai = coeff_amclass[iclass].w; // palpha[iclass];
    numtyp vali = polar6[i].x;

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
      int jclass = coeff_amtype[jtype].w;  // amtype2class[jtype];

      const numtyp4 sp_pol = sp_polar[sbmask15(jextra)];
      numtyp factor_mpole = sp_pol.w; // sp_mpole[sbmask15(jextra)];

      numtyp corek = coeff_amclass[jclass].z;  // pcore[jclass];
      numtyp alphak = coeff_amclass[jclass].w; // palpha[jclass];
      numtyp valk = polar6[j].x;

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
#ifdef INTEL_OCL
      numtyp t = ucl_recip((numtyp)1.0 + EWALD_P*ralpha);
      bn[0] = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * exp2a * rinv;
#else
      bn[0] = ucl_erfc(ralpha) * rinv;
#endif

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

      term1 = corei*corek;
      numtyp term1i = corek*vali;
      numtyp term2i = corek*dir;
      numtyp term3i = corek*qir;
      numtyp term1k = corei*valk;
      numtyp term2k = -corei*dkr;
      numtyp term3k = corei*qkr;
      numtyp term1ik = vali*valk;
      numtyp term2ik = valk*dir - vali*dkr + dik;
      numtyp term3ik = vali*qkr + valk*qir - dir*dkr + 2.0*(dkqi-diqk+qiqk);
      numtyp term4ik = dir*qkr - dkr*qir - 4.0*qik;
      numtyp term5ik = qir*qkr;
      numtyp dmpi[9],dmpj[9];
      numtyp dmpij[11];
      damppole(r,11,alphai,alphak,dmpi,dmpj,dmpij);
      numtyp scalek = factor_mpole;
      numtyp rr1i = bn[0] - ((numtyp)1.0-scalek*dmpi[0])*rr1;
      numtyp rr3i = bn[1] - ((numtyp)1.0-scalek*dmpi[2])*rr3;
      numtyp rr5i = bn[2] - ((numtyp)1.0-scalek*dmpi[4])*rr5;
      numtyp rr7i = bn[3] - ((numtyp)1.0-scalek*dmpi[6])*rr7;
      numtyp rr1k = bn[0] - ((numtyp)1.0-scalek*dmpj[0])*rr1;
      numtyp rr3k = bn[1] - ((numtyp)1.0-scalek*dmpj[2])*rr3;
      numtyp rr5k = bn[2] - ((numtyp)1.0-scalek*dmpj[4])*rr5;
      numtyp rr7k = bn[3] - ((numtyp)1.0-scalek*dmpj[6])*rr7;
      numtyp rr1ik = bn[0] - ((numtyp)1.0-scalek*dmpij[0])*rr1;
      numtyp rr3ik = bn[1] - ((numtyp)1.0-scalek*dmpij[2])*rr3;
      numtyp rr5ik = bn[2] - ((numtyp)1.0-scalek*dmpij[4])*rr5;
      numtyp rr7ik = bn[3] - ((numtyp)1.0-scalek*dmpij[6])*rr7;
      numtyp rr9ik = bn[4] - ((numtyp)1.0-scalek*dmpij[8])*rr9;
      numtyp rr11ik = bn[5] - ((numtyp)1.0-scalek*dmpij[10])*rr11;
      rr1 = bn[0] - ((numtyp)1.0-scalek)*rr1;
      rr3 = bn[1] - ((numtyp)1.0-scalek)*rr3;
      numtyp e = term1*rr1 + term4ik*rr7ik + term5ik*rr9ik +
        term1i*rr1i + term1k*rr1k + term1ik*rr1ik +
        term2i*rr3i + term2k*rr3k + term2ik*rr3ik +
        term3i*rr5i + term3k*rr5k + term3ik*rr5ik;

      // find damped multipole intermediates for force and torque

      numtyp de = term1*rr3 + term4ik*rr9ik + term5ik*rr11ik +
        term1i*rr3i + term1k*rr3k + term1ik*rr3ik +
        term2i*rr5i + term2k*rr5k + term2ik*rr5ik +
        term3i*rr7i + term3k*rr7k + term3ik*rr7ik;
      term1 = -corek*rr3i - valk*rr3ik + dkr*rr5ik - qkr*rr7ik;
      term2 = corei*rr3k + vali*rr3ik + dir*rr5ik + qir*rr7ik;
      term3 = (numtyp)2.0 * rr5ik;
      term4 = (numtyp)-2.0 * (corek*rr5i+valk*rr5ik - dkr*rr7ik+qkr*rr9ik);
      term5 = (numtyp)-2.0 * (corei*rr5k+vali*rr5ik + dir*rr7ik+qir*rr9ik);
      term6 = (numtyp)4.0 * rr7ik;
      rr3 = rr3ik;

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
  store_answers_hippo_tq(tq,ii,inum,tid,t_per_atom,offset,i,tep);

  // accumate force, energy and virial
  store_answers_acc(f,energy,e_coul,virial,ii,inum,tid,t_per_atom,
     offset,eflag,vflag,ans,engv,NUM_BLOCKS_X);
}

/* ----------------------------------------------------------------------
  udirect2b = Ewald real direct field via list
  udirect2b computes the real space contribution of the permanent
   atomic multipole moments to the field via a neighbor list
------------------------------------------------------------------------- */

__kernel void k_hippo_udirect2b(const __global numtyp4 *restrict x_,
                                const __global numtyp4 *restrict extra,
                                const __global numtyp4 *restrict coeff_amtype,
                                const __global numtyp4 *restrict coeff_amclass,
                                const __global numtyp4 *restrict sp_polar,
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
#if (SHUFFLE_AVAIL == 0)
  local_allocate_store_charge();
#endif

  acctyp _fieldp[6];
  for (int l=0; l<6; l++) _fieldp[l]=(acctyp)0;

  const __global numtyp4* polar1 = &extra[0];
  const __global numtyp4* polar2 = &extra[nall];
  const __global numtyp4* polar3 = &extra[2*nall];
  const __global numtyp4* polar6 = &extra[5*nall];

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
    int iclass = coeff_amtype[itype].w;  // amtype2class[itype];

    numtyp alphai = coeff_amclass[iclass].w; // palpha[iclass];

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
      numtyp rinv = ucl_recip(r);
      numtyp r2inv = rinv*rinv;
      numtyp rr1 = rinv;
      numtyp rr3 = rr1 * r2inv;
      numtyp rr5 = (numtyp)3.0 * rr3 * r2inv;
      numtyp rr7 = (numtyp)5.0 * rr5 * r2inv;

      const numtyp4 pol1j = polar1[j];
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
      int jclass = coeff_amtype[jtype].w;  // amtype2class[jtype];

      numtyp corek = coeff_amclass[jclass].z;  // pcore[jclass];
      numtyp alphak = coeff_amclass[jclass].w; // palpha[jclass];
      numtyp valk = polar6[j].x;

      numtyp factor_dscale, factor_pscale;
      const numtyp4 sp_pol = sp_polar[sbmask15(jextra)];
      if (igroup == jgroup) {
        factor_dscale = factor_pscale = sp_pol.y; // special_polar_piscale[sbmask15(jextra)];
      } else {
        factor_dscale = factor_pscale = sp_pol.z; // special_polar_pscale[sbmask15(jextra)];
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
      numtyp bn[4];
#ifdef INTEL_OCL
      numtyp t = ucl_recip((numtyp)1.0 + EWALD_P*ralpha);
      bn[0] = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * exp2a * rinv;
#else
      bn[0] = ucl_erfc(ralpha) * rinv;
#endif

      numtyp aefac = aesq2n;
      for (int m = 1; m <= 3; m++) {
        numtyp bfac = (numtyp) (m+m-1);
        aefac = aesq2 * aefac;
        bn[m] = (bfac*bn[m-1]+aefac*exp2a) * r2inv;
      }

      // find the field components for charge penetration damping
      numtyp dmpi[7],dmpk[7];
      dampdir(r,alphai,alphak,dmpi,dmpk);

      numtyp scalek = factor_dscale;
      numtyp rr3k = bn[1] - ((numtyp)1.0-scalek*dmpk[2])*rr3;
      numtyp rr5k = bn[2] - ((numtyp)1.0-scalek*dmpk[4])*rr5;
      numtyp rr7k = bn[3] - ((numtyp)1.0-scalek*dmpk[6])*rr7;
      rr3 = bn[1] - ((numtyp)1.0-scalek)*rr3;
      numtyp fid[3];
      fid[0] = -xr*(rr3*corek + rr3k*valk - rr5k*dkr + rr7k*qkr) -
        rr3k*dkx + (numtyp)2.0*rr5k*qkx;
      fid[1] = -yr*(rr3*corek + rr3k*valk - rr5k*dkr + rr7k*qkr) -
        rr3k*dky + (numtyp)2.0*rr5k*qky;
      fid[2] = -zr*(rr3*corek + rr3k*valk - rr5k*dkr + rr7k*qkr) -
        rr3k*dkz + (numtyp)2.0*rr5k*qkz;

      scalek = factor_pscale;
      rr3 = r2inv * rr1;
      rr3k = bn[1] - ((numtyp)1.0-scalek*dmpk[2])*rr3;
      rr5k = bn[2] - ((numtyp)1.0-scalek*dmpk[4])*rr5;
      rr7k = bn[3] - ((numtyp)1.0-scalek*dmpk[6])*rr7;
      rr3 = bn[1] - ((numtyp)1.0-scalek)*rr3;
      numtyp fip[3];
      fip[0] = -xr*(rr3*corek + rr3k*valk - rr5k*dkr + rr7k*qkr) -
        rr3k*dkx + (numtyp)2.0*rr5k*qkx;
      fip[1] = -yr*(rr3*corek + rr3k*valk - rr5k*dkr + rr7k*qkr) -
        rr3k*dky + (numtyp)2.0*rr5k*qky;
      fip[2] = -zr*(rr3*corek + rr3k*valk - rr5k*dkr + rr7k*qkr) -
        rr3k*dkz + (numtyp)2.0*rr5k*qkz;

      // find terms needed later to compute mutual polarization

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

__kernel void k_hippo_umutual2b(const __global numtyp4 *restrict x_,
                                const __global numtyp4 *restrict extra,
                                const __global numtyp4 *restrict coeff_amtype,
                                const __global numtyp4 *restrict coeff_amclass,
                                const __global numtyp4 *restrict sp_polar,
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
#if (SHUFFLE_AVAIL == 0)
  local_allocate_store_charge();
#endif

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

    int itype;
    itype  = polar3[i].z; // amtype[i];
    int iclass = coeff_amtype[itype].w;  // amtype2class[itype];
    numtyp alphai = coeff_amclass[iclass].w; // palpha[iclass];

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
      const numtyp4 pol4j = polar4[j];
      numtyp ukx = pol4j.x;  // uind[j][0];
      numtyp uky = pol4j.y;  // uind[j][1];
      numtyp ukz = pol4j.z;  // uind[j][2];
      const numtyp4 pol5j = polar5[j];
      numtyp ukxp = pol5j.x; // uinp[j][0];
      numtyp ukyp = pol5j.y; // uinp[j][1];
      numtyp ukzp = pol5j.z; // uinp[j][2];

      int jclass = coeff_amtype[jtype].w;  // amtype2class[jtype];
      numtyp alphak = coeff_amclass[jclass].w; // palpha[jclass];

      numtyp factor_wscale;
      const numtyp4 sp_pol = sp_polar[sbmask15(jextra)];
      factor_wscale = sp_pol.x; // special_polar_wscale[sbmask15(jextra)];

      // calculate the real space Ewald error function terms

      numtyp ralpha = aewald * r;
      numtyp exp2a = ucl_exp(-ralpha*ralpha);
      numtyp bn[4];
#ifdef INTEL_OCL
      numtyp t = ucl_recip((numtyp)1.0 + EWALD_P*ralpha);
      bn[0] = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * exp2a * rinv;
#else
      bn[0] = ucl_erfc(ralpha) * rinv;
#endif

      numtyp aefac = aesq2n;
      for (int m = 1; m <= 3; m++) {
        numtyp bfac = (numtyp) (m+m-1);
        aefac = aesq2 * aefac;
        bn[m] = (bfac*bn[m-1]+aefac*exp2a) * r2inv;
      }

      // find terms needed later to compute mutual polarization
      // if (poltyp != DIRECT)
      numtyp dmpik[5];
      dampmut(r,alphai,alphak,dmpik);
      numtyp scalek = factor_wscale;
      rr3 = r2inv * rr1;
      numtyp rr3ik = bn[1] - ((numtyp)1.0-scalek*dmpik[2])*rr3;
      numtyp rr5ik = bn[2] - ((numtyp)1.0-scalek*dmpik[4])*rr5;

      numtyp tdipdip[6];
      tdipdip[0] = -rr3ik + rr5ik*xr*xr;
      tdipdip[1] = rr5ik*xr*yr;
      tdipdip[2] = rr5ik*xr*zr;
      tdipdip[3] = -rr3ik + rr5ik*yr*yr;
      tdipdip[4] = rr5ik*yr*zr;
      tdipdip[5] = -rr3ik + rr5ik*zr*zr;

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

__kernel void k_hippo_polar(const __global numtyp4 *restrict x_,
                            const __global numtyp4 *restrict extra,
                            const __global numtyp4 *restrict coeff_amtype,
                            const __global numtyp4 *restrict coeff_amclass,
                            const __global numtyp4 *restrict sp_polar,
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
  const __global numtyp4* polar6 = &extra[5*nall];

  if (ii<inum) {
    int itype,igroup;
    numtyp uix,uiy,uiz,uixp,uiyp,uizp;

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

    numtyp corei = coeff_amclass[itype].z;  // pcore[iclass];
    numtyp alphai = coeff_amclass[itype].w; // palpha[iclass];
    numtyp vali = polar6[i].x;

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

      numtyp factor_wscale, factor_dscale;
      const numtyp4 sp_pol = sp_polar[sbmask15(jextra)];
      factor_wscale = sp_pol.x; // special_polar_wscale[sbmask15(jextra)];
      if (igroup == jgroup) {
        factor_dscale = sp_pol.y; // special_polar_piscale[sbmask15(jextra)];
      } else {
        factor_dscale = sp_pol.z; // special_polar_pscale[sbmask15(jextra)];
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

      // get reciprocal distance terms for this interaction

      numtyp rinv = ucl_recip(r);
      numtyp r2inv = rinv*rinv;
      numtyp rr1 = felec * rinv;
      numtyp rr3 = rr1 * r2inv;
      numtyp rr5 = (numtyp)3.0 * rr3 * r2inv;
      numtyp rr7 = (numtyp)5.0 * rr5 * r2inv;
      numtyp rr9 = (numtyp)7.0 * rr7 * r2inv;

      // calculate the real space Ewald error function terms

      int m;
      numtyp ralpha = aewald * r;
      numtyp exp2a = ucl_exp(-ralpha*ralpha);
      numtyp bn[5];
#ifdef INTEL_OCL
      numtyp t = ucl_recip((numtyp)1.0 + EWALD_P*ralpha);
      bn[0] = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * exp2a * rinv;
#else
      bn[0] = ucl_erfc(ralpha) * rinv;
#endif

      numtyp alsq2 = (numtyp)2.0 * aewald*aewald;
      numtyp alsq2n = (numtyp)0.0;
      if (aewald > (numtyp)0.0) alsq2n = (numtyp)1.0 / (MY_PIS*aewald);

      for (m = 1; m <= 4; m++) {
        numtyp bfac = (numtyp) (m+m-1);
        alsq2n = alsq2 * alsq2n;
        bn[m] = (bfac*bn[m-1]+alsq2n*exp2a) * r2inv;
      }
      for (m = 0; m < 5; m++) bn[m] *= felec;

      // apply charge penetration damping to scale factors

      numtyp corek = coeff_amclass[jtype].z;  // pcore[jclass];
      numtyp alphak = coeff_amclass[jtype].w; // palpha[jclass];
      numtyp valk = polar6[j].x;
      numtyp dmpi[9],dmpk[9];
      numtyp dmpik[9];
      damppole(r,9,alphai,alphak,dmpi,dmpk,dmpik);
      numtyp rr3core = bn[1] - ((numtyp)1.0-factor_dscale)*rr3;
      numtyp rr5core = bn[2] - ((numtyp)1.0-factor_dscale)*rr5;

      numtyp rr3i = bn[1] - ((numtyp)1.0-factor_dscale*dmpi[2])*rr3;
      numtyp rr5i = bn[2] - ((numtyp)1.0-factor_dscale*dmpi[4])*rr5;
      numtyp rr7i = bn[3] - ((numtyp)1.0-factor_dscale*dmpi[6])*rr7;
      numtyp rr9i = bn[4] - ((numtyp)1.0-factor_dscale*dmpi[8])*rr9;
      numtyp rr3k = bn[1] - ((numtyp)1.0-factor_dscale*dmpk[2])*rr3;
      numtyp rr5k = bn[2] - ((numtyp)1.0-factor_dscale*dmpk[4])*rr5;
      numtyp rr7k = bn[3] - ((numtyp)1.0-factor_dscale*dmpk[6])*rr7;
      numtyp rr9k = bn[4] - ((numtyp)1.0-factor_dscale*dmpk[8])*rr9;
      numtyp rr5ik = bn[2] - ((numtyp)1.0-factor_wscale*dmpik[4])*rr5;
      numtyp rr7ik = bn[3] - ((numtyp)1.0-factor_wscale*dmpik[6])*rr7;

      // get the induced dipole field used for dipole torques

      numtyp tix3 = (numtyp)2.0*rr3i*ukx;
      numtyp tiy3 = (numtyp)2.0*rr3i*uky;
      numtyp tiz3 = (numtyp)2.0*rr3i*ukz;
      numtyp tuir = (numtyp)-2.0*rr5i*ukr;

      ufld[0] += tix3 + xr*tuir;
      ufld[1] += tiy3 + yr*tuir;
      ufld[2] += tiz3 + zr*tuir;

      // get induced dipole field gradient used for quadrupole torques

      numtyp tix5 = (numtyp)4.0 * (rr5i*ukx);
      numtyp tiy5 = (numtyp)4.0 * (rr5i*uky);
      numtyp tiz5 = (numtyp)4.0 * (rr5i*ukz);
      tuir = (numtyp)-2.0*rr7i*ukr;

      dufld[0] += xr*tix5 + xr*xr*tuir;
      dufld[1] += xr*tiy5 + yr*tix5 + (numtyp)2.0*xr*yr*tuir;
      dufld[2] += yr*tiy5 + yr*yr*tuir;
      dufld[3] += xr*tiz5 + zr*tix5 + (numtyp)2.0*xr*zr*tuir;
      dufld[4] += yr*tiz5 + zr*tiy5 + (numtyp)2.0*yr*zr*tuir;
      dufld[5] += zr*tiz5 + zr*zr*tuir;

      // get the field gradient for direct polarization force

      numtyp term1i,term2i,term3i,term4i,term5i,term6i,term7i,term8i;
      numtyp term1k,term2k,term3k,term4k,term5k,term6k,term7k,term8k;
      numtyp term1core;
      numtyp tixx,tiyy,tizz,tixy,tixz,tiyz;
      numtyp tkxx,tkyy,tkzz,tkxy,tkxz,tkyz;

      term1i = rr3i - rr5i*xr*xr;
      term1core = rr3core - rr5core*xr*xr;
      term2i = (numtyp)2.0*rr5i*xr ;
      term3i = rr7i*xr*xr - rr5i;
      term4i = (numtyp)2.0*rr5i;
      term5i = (numtyp)5.0*rr7i*xr;
      term6i = rr9i*xr*xr;
      term1k = rr3k - rr5k*xr*xr;
      term2k = (numtyp)2.0*rr5k*xr;
      term3k = rr7k*xr*xr - rr5k;
      term4k = (numtyp)2.0*rr5k;
      term5k = (numtyp)5.0*rr7k*xr;
      term6k = rr9k*xr*xr;
      tixx = vali*term1i + corei*term1core + dix*term2i - dir*term3i -
        qixx*term4i + qix*term5i - qir*term6i + (qiy*yr+qiz*zr)*rr7i;
      tkxx = valk*term1k + corek*term1core - dkx*term2k + dkr*term3k -
        qkxx*term4k + qkx*term5k - qkr*term6k + (qky*yr+qkz*zr)*rr7k;

      term1i = rr3i - rr5i*yr*yr;
      term1core = rr3core - rr5core*yr*yr;
      term2i = (numtyp)2.0*rr5i*yr;
      term3i = rr7i*yr*yr - rr5i;
      term4i = (numtyp)2.0*rr5i;
      term5i = (numtyp)5.0*rr7i*yr;
      term6i = rr9i*yr*yr;
      term1k = rr3k - rr5k*yr*yr;
      term2k = (numtyp)2.0*rr5k*yr;
      term3k = rr7k*yr*yr - rr5k;
      term4k = (numtyp)2.0*rr5k;
      term5k = (numtyp)5.0*rr7k*yr;
      term6k = rr9k*yr*yr;
      tiyy = vali*term1i + corei*term1core + diy*term2i - dir*term3i -
        qiyy*term4i + qiy*term5i - qir*term6i + (qix*xr+qiz*zr)*rr7i;
      tkyy = valk*term1k + corek*term1core - dky*term2k + dkr*term3k -
        qkyy*term4k + qky*term5k - qkr*term6k + (qkx*xr+qkz*zr)*rr7k;

      term1i = rr3i - rr5i*zr*zr;
      term1core = rr3core - rr5core*zr*zr;
      term2i = (numtyp)2.0*rr5i*zr;
      term3i = rr7i*zr*zr - rr5i;
      term4i = (numtyp)2.0*rr5i;
      term5i = (numtyp)5.0*rr7i*zr;
      term6i = rr9i*zr*zr;
      term1k = rr3k - rr5k*zr*zr;
      term2k = (numtyp)2.0*rr5k*zr;
      term3k = rr7k*zr*zr - rr5k;
      term4k = (numtyp)2.0*rr5k;
      term5k = (numtyp)5.0*rr7k*zr;
      term6k = rr9k*zr*zr;
      tizz = vali*term1i + corei*term1core + diz*term2i - dir*term3i -
        qizz*term4i + qiz*term5i - qir*term6i + (qix*xr+qiy*yr)*rr7i;
      tkzz = valk*term1k + corek*term1core - dkz*term2k + dkr*term3k -
        qkzz*term4k + qkz*term5k - qkr*term6k + (qkx*xr+qky*yr)*rr7k;

      term2i = rr5i*xr ;
      term1i = yr * term2i;
      term1core = rr5core*xr*yr;
      term3i = rr5i*yr;
      term4i = yr * (rr7i*xr);
      term5i = (numtyp)2.0*rr5i;
      term6i = (numtyp)2.0*rr7i*xr;
      term7i = (numtyp)2.0*rr7i*yr;
      term8i = yr*rr9i*xr;
      term2k = rr5k*xr;
      term1k = yr * term2k;
      term3k = rr5k*yr;
      term4k = yr * (rr7k*xr);
      term5k = (numtyp)2.0*rr5k;
      term6k = (numtyp)2.0*rr7k*xr;
      term7k = (numtyp)2.0*rr7k*yr;
      term8k = yr*rr9k*xr;
      tixy = -vali*term1i - corei*term1core + diy*term2i + dix*term3i -
        dir*term4i - qixy*term5i + qiy*term6i + qix*term7i - qir*term8i;
      tkxy = -valk*term1k - corek*term1core - dky*term2k - dkx*term3k +
        dkr*term4k - qkxy*term5k + qky*term6k + qkx*term7k - qkr*term8k;

      term2i = rr5i*xr;
      term1i = zr * term2i;
      term1core = rr5core*xr*zr;
      term3i = rr5i*zr;
      term4i = zr * (rr7i*xr);
      term5i = (numtyp)2.0*rr5i;
      term6i = (numtyp)2.0*rr7i*xr;
      term7i = (numtyp)2.0*rr7i*zr;
      term8i = zr*rr9i*xr;
      term2k = rr5k*xr;
      term1k = zr * term2k;
      term3k = rr5k*zr;
      term4k = zr * (rr7k*xr);
      term5k = (numtyp)2.0*rr5k;
      term6k = (numtyp)2.0*rr7k*xr;
      term7k = (numtyp)2.0*rr7k*zr;
      term8k = zr*rr9k*xr;
      tixz = -vali*term1i - corei*term1core + diz*term2i + dix*term3i -
        dir*term4i - qixz*term5i + qiz*term6i + qix*term7i - qir*term8i;
      tkxz = -valk*term1k - corek*term1core - dkz*term2k - dkx*term3k +
        dkr*term4k - qkxz*term5k + qkz*term6k + qkx*term7k - qkr*term8k;

      term2i = rr5i*yr;
      term1i = zr * term2i;
      term1core = rr5core*yr*zr;
      term3i = rr5i*zr;
      term4i = zr * (rr7i*yr);
      term5i = (numtyp)2.0*rr5i;
      term6i = (numtyp)2.0*rr7i*yr;
      term7i = (numtyp)2.0*rr7i*zr;
      term8i = zr*rr9i*yr;
      term2k = rr5k*yr;
      term1k = zr * term2k;
      term3k = rr5k*zr;
      term4k = zr * (rr7k*yr);
      term5k = (numtyp)2.0*rr5k;
      term6k = (numtyp)2.0*rr7k*yr;
      term7k = (numtyp)2.0*rr7k*zr;
      term8k = zr*rr9k*yr;
      tiyz = -vali*term1i - corei*term1core + diz*term2i + diy*term3i -
        dir*term4i - qiyz*term5i + qiz*term6i + qiy*term7i - qir*term8i;
      tkyz = -valk*term1k - corek*term1core - dkz*term2k - dky*term3k +
        dkr*term4k - qkyz*term5k + qkz*term6k + qky*term7k - qkr*term8k;

      numtyp depx = tixx*ukx + tixy*uky + tixz*ukz - tkxx*uix - tkxy*uiy - tkxz*uiz;
      numtyp depy = tixy*ukx + tiyy*uky + tiyz*ukz - tkxy*uix - tkyy*uiy - tkyz*uiz;
      numtyp depz = tixz*ukx + tiyz*uky + tizz*ukz - tkxz*uix - tkyz*uiy - tkzz*uiz;

      numtyp frcx = (numtyp)-2.0 * depx;
      numtyp frcy = (numtyp)-2.0 * depy;
      numtyp frcz = (numtyp)-2.0 * depz;

      numtyp term1,term2,term3;

      // get the dEp/dR terms used for direct polarization force
      // poltyp == MUTUAL && hippo
      // tixx and tkxx
      term1 = (numtyp)2.0 * rr5ik;
      term2 = term1*xr;
      term3 = rr5ik - rr7ik*xr*xr;
      tixx = uix*term2 + uir*term3;
      tkxx = ukx*term2 + ukr*term3;

      // tiyy and tkyy
      term2 = term1*yr;
      term3 = rr5ik - rr7ik*yr*yr;
      tiyy = uiy*term2 + uir*term3;
      tkyy = uky*term2 + ukr*term3;

      // tiz and tkzz
      term2 = term1*zr;
      term3 = rr5ik - rr7ik*zr*zr;
      tizz = uiz*term2 + uir*term3;
      tkzz = ukz*term2 + ukr*term3;

      // tixy and tkxy
      term1 = rr5ik*yr;
      term2 = rr5ik*xr;
      term3 = yr * (rr7ik*xr);
      tixy = uix*term1 + uiy*term2 - uir*term3;
      tkxy = ukx*term1 + uky*term2 - ukr*term3;

      // tixx and tkxx
      term1 = rr5ik * zr;
      term3 = zr * (rr7ik*xr);
      tixz = uix*term1 + uiz*term2 - uir*term3;
      tkxz = ukx*term1 + ukz*term2 - ukr*term3;

      // tiyz and tkyz
      term2 = rr5ik*yr;
      term3 = zr * (rr7ik*yr);
      tiyz = uiy*term1 + uiz*term2 - uir*term3;
      tkyz = uky*term1 + ukz*term2 - ukr*term3;

      depx = tixx*ukxp + tixy*ukyp + tixz*ukzp + tkxx*uixp + tkxy*uiyp + tkxz*uizp;
      depy = tixy*ukxp + tiyy*ukyp + tiyz*ukzp + tkxy*uixp + tkyy*uiyp + tkyz*uizp;
      depz = tixz*ukxp + tiyz*ukyp + tizz*ukzp + tkxz*uixp + tkyz*uiyp + tkzz*uizp;

      frcx = frcx - depx;
      frcy = frcy - depy;
      frcz = frcz - depz;

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

__kernel void k_hippo_fphi_uind(const __global numtyp4 *restrict thetai1,
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

__kernel void k_hippo_fphi_mpole(const __global numtyp4 *restrict thetai1,
                          const __global numtyp4 *restrict thetai2,
                          const __global numtyp4 *restrict thetai3,
                          const __global int *restrict igrid,
                          const __global numtyp2 *restrict grid,
                          __global acctyp *restrict fphi,
                          const int bsorder, const int inum,
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
      fphi[idx] = buf[m];
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

__kernel void k_hippo_special15(__global int * dev_nbor,
                          const __global int * dev_packed,
                          const __global tagint *restrict tag,
                          const __global int *restrict nspecial15,
                          const __global tagint *restrict special15,
                          const int inum, const int nall, const int nbor_pitch,
                          const int t_per_atom) {
  int tid, ii, offset, n_stride, j;
  atom_info(t_per_atom,ii,tid,offset);

  if (ii<inum) {

    int numj, nbor, nbor_end;
    nbor_info(dev_nbor,dev_packed,nbor_pitch,t_per_atom,ii,offset,j,numj,
              n_stride,nbor_end,nbor);

    int n15 = nspecial15[ii];

    for ( ; nbor<nbor_end; nbor+=n_stride) {

      int sj=dev_packed[nbor];
      int which = sj >> SBBITS & 3;
      j = sj & NEIGHMASK;
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

__kernel void k_hippo_short_nbor(const __global numtyp4 *restrict x_,
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
