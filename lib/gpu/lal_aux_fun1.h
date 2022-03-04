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

#define nbor_info_p(nbor_mem, nbor_stride, t_per_atom, ii, offset,           \
                    i, numj, stride, nbor_end, nbor_begin)                   \
    i=nbor_mem[ii];                                                          \
    nbor_begin=ii+nbor_stride;                                               \
    numj=nbor_mem[nbor_begin];                                               \
    nbor_begin+=nbor_stride+ii*(t_per_atom-1);                               \
    stride=fast_mul(t_per_atom,nbor_stride);                                 \
    nbor_end=nbor_begin+fast_mul(numj/t_per_atom,stride)+(numj &             \
                                                          (t_per_atom-1));   \
    nbor_begin+=offset;

#if (SHUFFLE_AVAIL == 0)

#define simd_reduce_add1(width, local, offset, tid, one)                    \
  local[0][tid]=one;                                                        \
  for (unsigned int s=width/2; s>0; s>>=1) {                                \
    simdsync();                                                             \
    if (offset < s) local[0][tid] += local[0][tid+s];                       \
  }                                                                         \
  if (offset==0) one=local[0][tid];

#define simd_reduce_add2(width, local, offset, tid, one, two)               \
  local[0][tid]=one;                                                        \
  local[1][tid]=two;                                                        \
  for (unsigned int s=width/2; s>0; s>>=1) {                                \
    simdsync();                                                             \
    if (offset < s) {                                                       \
      local[0][tid] += local[0][tid+s];                                     \
      local[1][tid] += local[1][tid+s];                                     \
    }                                                                       \
  }                                                                         \
  if (offset==0) {                                                          \
    one=local[0][tid];                                                      \
    two=local[1][tid];                                                      \
  }

#define simd_reduce_add3(width, local, offset, tid, one, two, three)        \
  local[0][tid]=one;                                                        \
  local[1][tid]=two;                                                        \
  local[2][tid]=three;                                                      \
  for (unsigned int s=width/2; s>0; s>>=1) {                                \
    simdsync();                                                             \
    if (offset < s) {                                                       \
      local[0][tid] += local[0][tid+s];                                     \
      local[1][tid] += local[1][tid+s];                                     \
      local[2][tid] += local[2][tid+s];                                     \
    }                                                                       \
  }                                                                         \
  if (offset==0) {                                                          \
    one=local[0][tid];                                                      \
    two=local[1][tid];                                                      \
    three=local[2][tid];                                                    \
  }

#define simd_reduce_add6(width, local, offset, tid, one, two, three,        \
                         four, five, six)                                   \
  local[0][tid]=one;                                                        \
  local[1][tid]=two;                                                        \
  local[2][tid]=three;                                                      \
  local[3][tid]=four;                                                       \
  local[4][tid]=five;                                                       \
  local[5][tid]=six;                                                        \
  for (unsigned int s=width/2; s>0; s>>=1) {                                \
    simdsync();                                                             \
    if (offset < s) {                                                       \
      local[0][tid] += local[0][tid+s];                                     \
      local[1][tid] += local[1][tid+s];                                     \
      local[2][tid] += local[2][tid+s];                                     \
      local[3][tid] += local[3][tid+s];                                     \
      local[4][tid] += local[4][tid+s];                                     \
      local[5][tid] += local[5][tid+s];                                     \
    }                                                                       \
  }                                                                         \
  if (offset==0) {                                                          \
    one=local[0][tid];                                                      \
    two=local[1][tid];                                                      \
    three=local[2][tid];                                                    \
    four=local[3][tid];                                                     \
    five=local[4][tid];                                                     \
    six=local[5][tid];                                                      \
  }

#define simd_reduce_arr(trip, width, local, offset, tid, arr)               \
  for (int r=0; r<trip; r++)                                                \
    local[r][tid]=arr[r];                                                   \
  for (unsigned int s=width/2; s>0; s>>=1) {                                \
    simdsync();                                                             \
    if (offset < s) {                                                       \
      for (int r=0; r<trip; r++)                                            \
        local[r][tid] += local[r][tid+s];                                   \
    }                                                                       \
  }                                                                         \
  if (offset==0) {                                                          \
    for (int r=0; r<trip; r++)                                              \
      arr[r]=local[r][tid];                                                 \
  }

#define block_reduce_add1(width, local, tid, one)                           \
  local[0][tid]=one;                                                        \
  for (unsigned int s=BLOCK_SIZE_X/2; s>width/2; s>>=1) {                   \
    __syncthreads();                                                        \
    if (tid < s) local[0][tid] += local[0][tid+s];                          \
  }                                                                         \
  if (tid<width) {                                                          \
    for (unsigned int s=width/2; s>0; s>>=1) {                              \
      simdsync();                                                           \
      if (tid < s) local[0][tid] += local[0][tid+s];                        \
    }                                                                       \
    if (tid==0) one=local[0][tid];                                          \
  }

#define block_reduce_add2(width, local, tid, one, two)                      \
  local[0][tid]=one;                                                        \
  local[1][tid]=two;                                                        \
  for (unsigned int s=BLOCK_SIZE_X/2; s>width/2; s>>=1) {                   \
    __syncthreads();                                                        \
    if (tid < s) {                                                          \
      local[0][tid] += local[0][tid+s];                                     \
      local[1][tid] += local[1][tid+s];                                     \
    }                                                                       \
  }                                                                         \
  if (tid<width) {                                                          \
    for (unsigned int s=width/2; s>0; s>>=1) {                              \
      simdsync();                                                           \
      if (tid < s) {                                                        \
        local[0][tid] += local[0][tid+s];                                   \
        local[1][tid] += local[1][tid+s];                                   \
      }                                                                     \
    }                                                                       \
    if (tid==0) {                                                           \
      one=local[0][tid];                                                    \
      two=local[1][tid];                                                    \
    }                                                                       \
  }

#define block_reduce_arr(trip, width, local, tid, arr)                      \
  for (int r=0; r<trip; r++)                                                \
    local[r][tid]=arr[r];                                                   \
  for (unsigned int s=BLOCK_SIZE_X/2; s>width/2; s>>=1) {                   \
    __syncthreads();                                                        \
    if (tid < s) {                                                          \
      for (int r=0; r<trip; r++)                                            \
        local[r][tid] += local[r][tid+s];                                   \
    }                                                                       \
  }                                                                         \
  if (tid<width) {                                                          \
    for (unsigned int s=width/2; s>0; s>>=1) {                              \
      simdsync();                                                           \
      if (tid < s) {                                                        \
        for (int r=0; r<trip; r++)                                          \
          local[r][tid] += local[r][tid+s];                                 \
      }                                                                     \
    }                                                                       \
    if (tid==0) {                                                           \
      for (int r=0; r<trip; r++)                                            \
        arr[r]=local[r][tid];                                               \
    }                                                                       \
  }

#define local_allocate_store_pair()                                         \
    __local acctyp red_acc[6][BLOCK_PAIR];
#define local_allocate_store_charge()                                       \
    __local acctyp red_acc[6][BLOCK_PAIR];
#define local_allocate_store_bio()                                          \
    __local acctyp red_acc[6][BLOCK_BIO_PAIR];
#define local_allocate_store_ellipse()                                      \
    __local acctyp red_acc[6][BLOCK_ELLIPSE];
#define local_allocate_store_three()                                        \
    __local acctyp red_acc[6][BLOCK_ELLIPSE];

#define store_answers(f, energy, virial, ii, inum, tid,                     \
                        t_per_atom, offset, eflag, vflag, ans, engv)        \
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
  if (offset==0 && ii<inum) ans[ii]=f;                                      \
  if (EVFLAG && (eflag || vflag)) {                                         \
    int ei=BLOCK_ID_X;                                                      \
    if (eflag!=2 && vflag!=2) {                                             \
      const int ev_stride=NUM_BLOCKS_X;                                     \
      if (eflag) {                                                          \
        simdsync();                                                         \
        block_reduce_add1(simd_size(), red_acc, tid, energy);               \
        if (vflag) __syncthreads();                                         \
        if (tid==0) {                                                       \
          engv[ei]=energy*(acctyp)0.5;                                      \
          ei+=ev_stride;                                                    \
        }                                                                   \
      }                                                                     \
      if (vflag) {                                                          \
        simdsync();                                                         \
        block_reduce_arr(6, simd_size(), red_acc, tid, virial);             \
        if (tid==0) {                                                       \
          for (int r=0; r<6; r++) {                                         \
            engv[ei]=virial[r]*(acctyp)0.5;                                 \
            ei+=ev_stride;                                                  \
          }                                                                 \
        }                                                                   \
      }                                                                     \
    } else if (offset==0 && ii<inum) {                                      \
      int ei=ii;                                                            \
      if (EVFLAG && eflag) {                                                \
        engv[ei]=energy*(acctyp)0.5;                                        \
        ei+=inum;                                                           \
      }                                                                     \
      if (EVFLAG && vflag) {                                                \
        for (int i=0; i<6; i++) {                                           \
          engv[ei]=virial[i]*(acctyp)0.5;                                   \
          ei+=inum;                                                         \
        }                                                                   \
      }                                                                     \
    }                                                                       \
  }

#define store_answers_q(f, energy, e_coul, virial, ii, inum, tid,           \
                        t_per_atom, offset, eflag, vflag, ans, engv)        \
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
  if (offset==0 && ii<inum) ans[ii]=f;                                      \
  if (EVFLAG && (eflag || vflag)) {                                         \
    int ei=BLOCK_ID_X;                                                      \
    const int ev_stride=NUM_BLOCKS_X;                                       \
    if (eflag!=2 && vflag!=2) {                                             \
      if (eflag) {                                                          \
        simdsync();                                                         \
        block_reduce_add2(simd_size(), red_acc, tid, energy, e_coul);       \
        if (vflag) __syncthreads();                                         \
        if (tid==0) {                                                       \
          engv[ei]=energy*(acctyp)0.5;                                      \
          ei+=ev_stride;                                                    \
          engv[ei]=e_coul*(acctyp)0.5;                                      \
          ei+=ev_stride;                                                    \
        }                                                                   \
      }                                                                     \
      if (vflag) {                                                          \
        simdsync();                                                         \
        block_reduce_arr(6, simd_size(), red_acc, tid, virial);             \
        if (tid==0) {                                                       \
          for (int r=0; r<6; r++) {                                         \
            engv[ei]=virial[r]*(acctyp)0.5;                                 \
            ei+=ev_stride;                                                  \
          }                                                                 \
        }                                                                   \
      }                                                                     \
    } else if (offset==0 && ii<inum) {                                      \
      int ei=ii;                                                            \
      if (EVFLAG && eflag) {                                                \
        engv[ei]=energy*(acctyp)0.5;                                        \
        ei+=inum;                                                           \
        engv[ei]=e_coul*(acctyp)0.5;                                        \
        ei+=inum;                                                           \
      }                                                                     \
      if (EVFLAG && vflag) {                                                \
        for (int i=0; i<6; i++) {                                           \
          engv[ei]=virial[i]*(acctyp)0.5;                                   \
          ei+=inum;                                                         \
        }                                                                   \
      }                                                                     \
    }                                                                       \
  }

#else

#define simd_reduce_add1(width, one)                                        \
  for (unsigned int s=width/2; s>0; s>>=1) one += shfl_down(one, s, width);

#define simd_reduce_add2(width, one, two)                                   \
  for (unsigned int s=width/2; s>0; s>>=1) {                                \
    one += shfl_down(one, s, width);                                        \
    two += shfl_down(two, s, width);                                        \
  }

#define simd_reduce_add3(width, one, two, three)                            \
  for (unsigned int s=width/2; s>0; s>>=1) {                                \
    one += shfl_down(one, s, width);                                        \
    two += shfl_down(two, s, width);                                        \
    three += shfl_down(three, s, width);                                    \
  }

#define simd_reduce_add6(width, one, two, three, four, five, six)           \
  for (unsigned int s=width/2; s>0; s>>=1) {                                \
    one += shfl_down(one, s, width);                                        \
    two += shfl_down(two, s, width);                                        \
    three += shfl_down(three, s, width);                                    \
    four += shfl_down(four, s, width);                                      \
    five += shfl_down(five, s, width);                                      \
    six += shfl_down(six, s, width);                                        \
  }

#define simd_reduce_arr(trip, width, arr)                                   \
  for (unsigned int s=width/2; s>0; s>>=1) {                                \
    for (int r=0; r<trip; r++)                                              \
      arr[r] += shfl_down(arr[r], s, width);                                \
  }

#if (EVFLAG == 1)

#define local_allocate_store_pair()                                         \
    __local acctyp red_acc[7][BLOCK_PAIR / SIMD_SIZE];
#define local_allocate_store_charge()                                       \
    __local acctyp red_acc[8][BLOCK_PAIR / SIMD_SIZE];
#define local_allocate_store_bio()                                          \
    __local acctyp red_acc[8][BLOCK_BIO_PAIR / SIMD_SIZE];
#define local_allocate_store_ellipse()
#define local_allocate_store_three()                                        \
    __local acctyp red_acc[7][BLOCK_ELLIPSE / SIMD_SIZE];

#define store_answers(f, energy, virial, ii, inum, tid,                     \
                      t_per_atom, offset, eflag, vflag, ans, engv)          \
  if (t_per_atom>1) {                                                       \
    simd_reduce_add3(t_per_atom, f.x, f.y, f.z);                            \
    if (vflag==2 || eflag==2) {                                             \
      if (eflag)                                                            \
        simd_reduce_add1(t_per_atom,energy);                                \
      if (vflag)                                                            \
        simd_reduce_arr(6, t_per_atom,virial);                              \
    }                                                                       \
  }                                                                         \
  if (offset==0 && ii<inum) ans[ii]=f;                                      \
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
        const int ev_stride=NUM_BLOCKS_X;                                   \
        if (eflag) {                                                        \
          simd_reduce_add1(vwidth, energy);                                 \
          if (tid==0) {                                                     \
            engv[ei]=energy*(acctyp)0.5;                                    \
            ei+=ev_stride;                                                  \
          }                                                                 \
        }                                                                   \
        if (vflag) {                                                        \
          simd_reduce_arr(6, vwidth, virial);                               \
          if (tid==0) {                                                     \
            for (int r=0; r<6; r++) {                                       \
              engv[ei]=virial[r]*(acctyp)0.5;                               \
              ei+=ev_stride;                                                \
            }                                                               \
          }                                                                 \
        }                                                                   \
      }                                                                     \
    } else if (offset==0 && ii<inum) {                                      \
      int ei=ii;                                                            \
      if (eflag) {                                                          \
        engv[ei]=energy*(acctyp)0.5;                                        \
        ei+=inum;                                                           \
      }                                                                     \
      if (vflag) {                                                          \
        for (int i=0; i<6; i++) {                                           \
          engv[ei]=virial[i]*(acctyp)0.5;                                   \
          ei+=inum;                                                         \
        }                                                                   \
      }                                                                     \
    }                                                                       \
  }

#define store_answers_q(f, energy, e_coul, virial, ii, inum, tid,           \
                        t_per_atom, offset, eflag, vflag, ans, engv)        \
  if (t_per_atom>1) {                                                       \
    simd_reduce_add3(t_per_atom, f.x, f.y, f.z);                            \
    if (vflag==2 || eflag==2) {                                             \
      if (eflag)                                                            \
        simd_reduce_add2(t_per_atom,energy,e_coul);                         \
      if (vflag)                                                            \
        simd_reduce_arr(6, t_per_atom,virial);                              \
    }                                                                       \
  }                                                                         \
  if (offset==0 && ii<inum) ans[ii]=f;                                      \
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
        const int ev_stride=NUM_BLOCKS_X;                                   \
        if (eflag) {                                                        \
          simd_reduce_add2(vwidth, energy, e_coul);                         \
          if (tid==0) {                                                     \
            engv[ei]=energy*(acctyp)0.5;                                    \
            ei+=ev_stride;                                                  \
            engv[ei]=e_coul*(acctyp)0.5;                                    \
            ei+=ev_stride;                                                  \
          }                                                                 \
        }                                                                   \
        if (vflag) {                                                        \
          simd_reduce_arr(6, vwidth, virial);                               \
          if (tid==0) {                                                     \
            for (int r=0; r<6; r++) {                                       \
              engv[ei]=virial[r]*(acctyp)0.5;                               \
              ei+=ev_stride;                                                \
            }                                                               \
          }                                                                 \
        }                                                                   \
      }                                                                     \
    } else if (offset==0 && ii<inum) {                                      \
      int ei=ii;                                                            \
      if (eflag) {                                                          \
        engv[ei]=energy*(acctyp)0.5;                                        \
        ei+=inum;                                                           \
        engv[ei]=e_coul*(acctyp)0.5;                                        \
        ei+=inum;                                                           \
      }                                                                     \
      if (vflag) {                                                          \
        for (int i=0; i<6; i++) {                                           \
          engv[ei]=virial[i]*(acctyp)0.5;                                   \
          ei+=inum;                                                         \
        }                                                                   \
      }                                                                     \
    }                                                                       \
  }

#else

#define local_allocate_store_pair()
#define local_allocate_store_charge()
#define local_allocate_store_bio()
#define local_allocate_store_ellipse()
#define local_allocate_store_three()

#define store_answers(f, energy, virial, ii, inum, tid,                     \
                      t_per_atom, offset, eflag, vflag, ans, engv)          \
  if (t_per_atom>1)                                                         \
    simd_reduce_add3(t_per_atom, f.x, f.y, f.z);                            \
  if (offset==0 && ii<inum) ans[ii]=f;

#define store_answers_q(f, energy, e_coul, virial, ii, inum, tid,           \
                        t_per_atom, offset, eflag, vflag, ans, engv)        \
  if (t_per_atom>1)                                                         \
    simd_reduce_add3(t_per_atom, f.x, f.y, f.z);                            \
  if (offset==0 && ii<inum) ans[ii]=f;

#endif

#endif

