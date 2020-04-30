// **************************************************************************
//                                 aux_fun1.h
//                             -------------------
//                           W. Michael Brown (ORNL)
//
//  Device code for pair style auxiliary functions
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                : Sat Oct 22 2011
//    email                : brownw@ornl.gov
// ***************************************************************************/

#if defined(NV_KERNEL) || defined(USE_HIP)
#include "lal_preprocessor.h"
#endif

#define atom_info(t_per_atom, ii, tid, offset)                               \
  tid=THREAD_ID_X;                                                           \
  offset=tid & (t_per_atom-1);                                               \
  ii=fast_mul((int)BLOCK_ID_X,(int)(BLOCK_SIZE_X)/t_per_atom)+tid/t_per_atom;

#define nbor_info(dev_nbor, dev_packed, nbor_pitch, t_per_atom, ii, offset,  \
                  i, numj, n_stride, nbor_end, nbor_begin)                   \
  i=dev_nbor[ii];                                                            \
  nbor_begin=ii+nbor_pitch;                                                  \
  numj=dev_nbor[nbor_begin];                                                 \
  if (dev_nbor==dev_packed) {                                                \
    nbor_begin+=nbor_pitch+fast_mul(ii,t_per_atom-1);                        \
    n_stride=fast_mul(t_per_atom,nbor_pitch);                                \
    nbor_end=nbor_begin+fast_mul(numj/t_per_atom,n_stride)+(numj & (t_per_atom-1)); \
    nbor_begin+=offset;                                                      \
  } else {                                                                   \
    nbor_begin+=nbor_pitch;                                                  \
    nbor_begin=dev_nbor[nbor_begin];                                         \
    nbor_end=nbor_begin+numj;                                                \
    n_stride=t_per_atom;                                                     \
    nbor_begin+=offset;                                                      \
  }

#if (ARCH < 300)

#define store_answers(f, energy, virial, ii, inum, tid, t_per_atom, offset, \
                      eflag, vflag, ans, engv)                              \
  if (t_per_atom>1) {                                                       \
    __local acctyp red_acc[6][BLOCK_PAIR];                                  \
    red_acc[0][tid]=f.x;                                                    \
    red_acc[1][tid]=f.y;                                                    \
    red_acc[2][tid]=f.z;                                                    \
    red_acc[3][tid]=energy;                                                 \
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                         \
      if (offset < s) {                                                     \
        for (int r=0; r<4; r++)                                             \
          red_acc[r][tid] += red_acc[r][tid+s];                             \
      }                                                                     \
    }                                                                       \
    f.x=red_acc[0][tid];                                                    \
    f.y=red_acc[1][tid];                                                    \
    f.z=red_acc[2][tid];                                                    \
    energy=red_acc[3][tid];                                                 \
    if (vflag>0) {                                                          \
      for (int r=0; r<6; r++)                                               \
        red_acc[r][tid]=virial[r];                                          \
      for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                       \
        if (offset < s) {                                                   \
          for (int r=0; r<6; r++)                                           \
            red_acc[r][tid] += red_acc[r][tid+s];                           \
        }                                                                   \
      }                                                                     \
      for (int r=0; r<6; r++)                                               \
        virial[r]=red_acc[r][tid];                                          \
    }                                                                       \
  }                                                                         \
  if (offset==0) {                                                          \
    int ei=ii;                                                              \
    if (eflag>0) {                                                          \
      engv[ei]=energy*(acctyp)0.5;                                          \
      ei+=inum;                                                             \
    }                                                                       \
    if (vflag>0) {                                                          \
      for (int i=0; i<6; i++) {                                             \
        engv[ei]=virial[i]*(acctyp)0.5;                                     \
        ei+=inum;                                                           \
      }                                                                     \
    }                                                                       \
    ans[ii]=f;                                                              \
  }

#define store_answers_q(f, energy, e_coul, virial, ii, inum, tid,           \
                        t_per_atom, offset, eflag, vflag, ans, engv)        \
  if (t_per_atom>1) {                                                       \
    __local acctyp red_acc[6][BLOCK_PAIR];                                  \
    red_acc[0][tid]=f.x;                                                    \
    red_acc[1][tid]=f.y;                                                    \
    red_acc[2][tid]=f.z;                                                    \
    red_acc[3][tid]=energy;                                                 \
    red_acc[4][tid]=e_coul;                                                 \
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                         \
      if (offset < s) {                                                     \
        for (int r=0; r<5; r++)                                             \
          red_acc[r][tid] += red_acc[r][tid+s];                             \
      }                                                                     \
    }                                                                       \
    f.x=red_acc[0][tid];                                                    \
    f.y=red_acc[1][tid];                                                    \
    f.z=red_acc[2][tid];                                                    \
    energy=red_acc[3][tid];                                                 \
    e_coul=red_acc[4][tid];                                                 \
    if (vflag>0) {                                                          \
      for (int r=0; r<6; r++)                                               \
        red_acc[r][tid]=virial[r];                                          \
      for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                       \
        if (offset < s) {                                                   \
          for (int r=0; r<6; r++)                                           \
            red_acc[r][tid] += red_acc[r][tid+s];                           \
        }                                                                   \
      }                                                                     \
      for (int r=0; r<6; r++)                                               \
        virial[r]=red_acc[r][tid];                                          \
    }                                                                       \
  }                                                                         \
  if (offset==0) {                                                          \
    int ei=ii;                                                              \
    if (eflag>0) {                                                          \
      engv[ei]=energy*(acctyp)0.5;                                          \
      ei+=inum;                                                             \
      engv[ei]=e_coul*(acctyp)0.5;                                          \
      ei+=inum;                                                             \
    }                                                                       \
    if (vflag>0) {                                                          \
      for (int i=0; i<6; i++) {                                             \
        engv[ei]=virial[i]*(acctyp)0.5;                                     \
        ei+=inum;                                                           \
      }                                                                     \
    }                                                                       \
    ans[ii]=f;                                                              \
  }

#else

#define store_answers(f, energy, virial, ii, inum, tid, t_per_atom, offset, \
                      eflag, vflag, ans, engv)                              \
  if (t_per_atom>1) {                                                       \
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                         \
        f.x += shfl_xor(f.x, s, t_per_atom);                                \
        f.y += shfl_xor(f.y, s, t_per_atom);                                \
        f.z += shfl_xor(f.z, s, t_per_atom);                                \
        energy += shfl_xor(energy, s, t_per_atom);                          \
    }                                                                       \
    if (vflag>0) {                                                          \
      for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                       \
          for (int r=0; r<6; r++)                                           \
            virial[r] += shfl_xor(virial[r], s, t_per_atom);                \
      }                                                                     \
    }                                                                       \
  }                                                                         \
  if (offset==0) {                                                          \
    int ei=ii;                                                              \
    if (eflag>0) {                                                          \
      engv[ei]=energy*(acctyp)0.5;                                          \
      ei+=inum;                                                             \
    }                                                                       \
    if (vflag>0) {                                                          \
      for (int i=0; i<6; i++) {                                             \
        engv[ei]=virial[i]*(acctyp)0.5;                                     \
        ei+=inum;                                                           \
      }                                                                     \
    }                                                                       \
    ans[ii]=f;                                                              \
  }

#define store_answers_q(f, energy, e_coul, virial, ii, inum, tid,           \
                        t_per_atom, offset, eflag, vflag, ans, engv)        \
  if (t_per_atom>1) {                                                       \
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                         \
      f.x += shfl_xor(f.x, s, t_per_atom);                                  \
      f.y += shfl_xor(f.y, s, t_per_atom);                                  \
      f.z += shfl_xor(f.z, s, t_per_atom);                                  \
      energy += shfl_xor(energy, s, t_per_atom);                            \
      e_coul += shfl_xor(e_coul, s, t_per_atom);                            \
    }                                                                       \
    if (vflag>0) {                                                          \
      for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                       \
          for (int r=0; r<6; r++)                                           \
            virial[r] += shfl_xor(virial[r], s, t_per_atom);                \
      }                                                                     \
    }                                                                       \
  }                                                                         \
  if (offset==0) {                                                          \
    int ei=ii;                                                              \
    if (eflag>0) {                                                          \
      engv[ei]=energy*(acctyp)0.5;                                          \
      ei+=inum;                                                             \
      engv[ei]=e_coul*(acctyp)0.5;                                          \
      ei+=inum;                                                             \
    }                                                                       \
    if (vflag>0) {                                                          \
      for (int i=0; i<6; i++) {                                             \
        engv[ei]=virial[i]*(acctyp)0.5;                                     \
        ei+=inum;                                                           \
      }                                                                     \
    }                                                                       \
    ans[ii]=f;                                                              \
  }

#endif

