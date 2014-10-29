// **************************************************************************
//                              ellipsoid_extra.h
//                             -------------------
//                           W. Michael Brown (ORNL)
//
//  Device code for Ellipsoid math routines
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                :
//    email                : brownw@ornl.gov
// ***************************************************************************/

#ifndef LAL_ELLIPSOID_EXTRA_H
#define LAL_ELLIPSOID_EXTRA_H

enum{SPHERE_SPHERE,SPHERE_ELLIPSE,ELLIPSE_SPHERE,ELLIPSE_ELLIPSE};

#ifdef NV_KERNEL
#include "lal_aux_fun1.h"
#ifndef _DOUBLE_DOUBLE
texture<float4> pos_tex, quat_tex;
#else
texture<int4,1> pos_tex, quat_tex;
#endif
#else
#define pos_tex x_
#define quat_tex qif
#endif

#define nbor_info_e(nbor_mem, nbor_stride, t_per_atom, ii, offset,           \
                    i, numj, stride, nbor_end, nbor_begin)                   \
    i=nbor_mem[ii];                                                          \
    nbor_begin=ii+nbor_stride;                                               \
    numj=nbor_mem[nbor_begin];                                               \
    nbor_begin+=nbor_stride;                                                 \
    nbor_end=nbor_begin+fast_mul(nbor_stride,numj);                          \
    nbor_begin+=fast_mul(offset,nbor_stride);                                \
    stride=fast_mul(t_per_atom,nbor_stride);

#if (ARCH < 300)

#define store_answers_t(f, tor, energy, virial, ii, astride, tid,           \
                        t_per_atom, offset, eflag, vflag, ans, engv)        \
  if (t_per_atom>1) {                                                       \
    __local acctyp red_acc[7][BLOCK_PAIR];                                  \
    red_acc[0][tid]=f.x;                                                    \
    red_acc[1][tid]=f.y;                                                    \
    red_acc[2][tid]=f.z;                                                    \
    red_acc[3][tid]=tor.x;                                                  \
    red_acc[4][tid]=tor.y;                                                  \
    red_acc[5][tid]=tor.z;                                                  \
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                         \
      if (offset < s) {                                                     \
        for (int r=0; r<6; r++)                                             \
          red_acc[r][tid] += red_acc[r][tid+s];                             \
      }                                                                     \
    }                                                                       \
    f.x=red_acc[0][tid];                                                    \
    f.y=red_acc[1][tid];                                                    \
    f.z=red_acc[2][tid];                                                    \
    tor.x=red_acc[3][tid];                                                  \
    tor.y=red_acc[4][tid];                                                  \
    tor.z=red_acc[5][tid];                                                  \
    if (eflag>0 || vflag>0) {                                               \
      for (int r=0; r<6; r++)                                               \
        red_acc[r][tid]=virial[r];                                          \
      red_acc[6][tid]=energy;                                               \
      for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                       \
        if (offset < s) {                                                   \
          for (int r=0; r<7; r++)                                           \
            red_acc[r][tid] += red_acc[r][tid+s];                           \
        }                                                                   \
      }                                                                     \
      for (int r=0; r<6; r++)                                               \
        virial[r]=red_acc[r][tid];                                          \
      energy=red_acc[6][tid];                                               \
    }                                                                       \
  }                                                                         \
  if (offset==0) {                                                          \
    __global acctyp *ap1=engv+ii;                                           \
    if (eflag>0) {                                                          \
      *ap1=energy*(acctyp)0.5;                                              \
      ap1+=astride;                                                         \
    }                                                                       \
    if (vflag>0) {                                                          \
      for (int i=0; i<6; i++) {                                             \
        *ap1=virial[i]*(acctyp)0.5;                                         \
        ap1+=astride;                                                       \
      }                                                                     \
    }                                                                       \
    ans[ii]=f;                                                              \
    ans[ii+astride]=tor;                                                    \
  }

#define acc_answers(f, energy, virial, ii, inum, tid, t_per_atom, offset,   \
                    eflag, vflag, ans, engv)                                \
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
    engv+=ii;                                                               \
    if (eflag>0) {                                                          \
      *engv+=energy*(acctyp)0.5;                                            \
      engv+=inum;                                                           \
    }                                                                       \
    if (vflag>0) {                                                          \
      for (int i=0; i<6; i++) {                                             \
        *engv+=virial[i]*(acctyp)0.5;                                       \
        engv+=inum;                                                         \
      }                                                                     \
    }                                                                       \
    acctyp4 old=ans[ii];                                                    \
    old.x+=f.x;                                                             \
    old.y+=f.y;                                                             \
    old.z+=f.z;                                                             \
    ans[ii]=old;                                                            \
  }

#else

#define store_answers_t(f, tor, energy, virial, ii, astride, tid,           \
                        t_per_atom, offset, eflag, vflag, ans, engv)        \
  if (t_per_atom>1) {                                                       \
    for (unsigned int s=t_per_atom/2; s>0; s>>=1) {                         \
        f.x += shfl_xor(f.x, s, t_per_atom);                                \
        f.y += shfl_xor(f.y, s, t_per_atom);                                \
        f.z += shfl_xor(f.z, s, t_per_atom);                                \
        tor.x += shfl_xor(tor.x, s, t_per_atom);                            \
        tor.y += shfl_xor(tor.y, s, t_per_atom);                            \
        tor.z += shfl_xor(tor.z, s, t_per_atom);                            \
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
    __global acctyp *ap1=engv+ii;                                           \
    if (eflag>0) {                                                          \
      *ap1=energy*(acctyp)0.5;                                              \
      ap1+=astride;                                                         \
    }                                                                       \
    if (vflag>0) {                                                          \
      for (int i=0; i<6; i++) {                                             \
        *ap1=virial[i]*(acctyp)0.5;                                         \
        ap1+=astride;                                                       \
      }                                                                     \
    }                                                                       \
    ans[ii]=f;                                                              \
    ans[ii+astride]=tor;                                                    \
  }

#define acc_answers(f, energy, virial, ii, inum, tid, t_per_atom, offset,   \
                    eflag, vflag, ans, engv)                                \
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
    engv+=ii;                                                               \
    if (eflag>0) {                                                          \
      *engv+=energy*(acctyp)0.5;                                            \
      engv+=inum;                                                           \
    }                                                                       \
    if (vflag>0) {                                                          \
      for (int i=0; i<6; i++) {                                             \
        *engv+=virial[i]*(acctyp)0.5;                                       \
        engv+=inum;                                                         \
      }                                                                     \
    }                                                                       \
    acctyp4 old=ans[ii];                                                    \
    old.x+=f.x;                                                             \
    old.y+=f.y;                                                             \
    old.z+=f.z;                                                             \
    ans[ii]=old;                                                            \
  }

#endif

/* ----------------------------------------------------------------------
   dot product of 2 vectors
------------------------------------------------------------------------- */

ucl_inline numtyp gpu_dot3(const numtyp *v1, const numtyp *v2)
{
  return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
};

/* ----------------------------------------------------------------------
   cross product of 2 vectors
------------------------------------------------------------------------- */

ucl_inline void gpu_cross3(const numtyp *v1, const numtyp *v2, numtyp *ans)
{
  ans[0] = v1[1]*v2[2]-v1[2]*v2[1];
  ans[1] = v1[2]*v2[0]-v1[0]*v2[2];
  ans[2] = v1[0]*v2[1]-v1[1]*v2[0];
};

/* ----------------------------------------------------------------------
   determinant of a matrix
------------------------------------------------------------------------- */

ucl_inline numtyp gpu_det3(const numtyp m[9])
{
  numtyp ans = m[0]*m[4]*m[8] - m[0]*m[5]*m[7] - 
    m[3]*m[1]*m[8] + m[3]*m[2]*m[7] + 
    m[6]*m[1]*m[5] - m[6]*m[2]*m[4];
  return ans;
};

/* ----------------------------------------------------------------------
   diagonal matrix times a full matrix
------------------------------------------------------------------------- */

ucl_inline void gpu_diag_times3(const numtyp4 shape, const numtyp m[9], 
                              numtyp ans[9])
{
  ans[0] = shape.x*m[0];
  ans[1] = shape.x*m[1];
  ans[2] = shape.x*m[2];
  ans[3] = shape.y*m[3];
  ans[4] = shape.y*m[4];
  ans[5] = shape.y*m[5];
  ans[6] = shape.z*m[6];
  ans[7] = shape.z*m[7];
  ans[8] = shape.z*m[8];
};

/* ----------------------------------------------------------------------
   add two matrices
------------------------------------------------------------------------- */

ucl_inline void gpu_plus3(const numtyp m[9], const numtyp m2[9], numtyp ans[9])
{
  ans[0] = m[0]+m2[0];
  ans[1] = m[1]+m2[1];
  ans[2] = m[2]+m2[2];
  ans[3] = m[3]+m2[3];
  ans[4] = m[4]+m2[4];
  ans[5] = m[5]+m2[5];
  ans[6] = m[6]+m2[6];
  ans[7] = m[7]+m2[7];
  ans[8] = m[8]+m2[8];
};

/* ----------------------------------------------------------------------
   multiply the transpose of mat1 times mat2
------------------------------------------------------------------------- */

ucl_inline void gpu_transpose_times3(const numtyp m[9], const numtyp m2[9],
                                   numtyp ans[9])
{
  ans[0] = m[0]*m2[0]+m[3]*m2[3]+m[6]*m2[6];
  ans[1] = m[0]*m2[1]+m[3]*m2[4]+m[6]*m2[7];
  ans[2] = m[0]*m2[2]+m[3]*m2[5]+m[6]*m2[8];
  ans[3] = m[1]*m2[0]+m[4]*m2[3]+m[7]*m2[6];
  ans[4] = m[1]*m2[1]+m[4]*m2[4]+m[7]*m2[7];
  ans[5] = m[1]*m2[2]+m[4]*m2[5]+m[7]*m2[8];
  ans[6] = m[2]*m2[0]+m[5]*m2[3]+m[8]*m2[6];
  ans[7] = m[2]*m2[1]+m[5]*m2[4]+m[8]*m2[7];
  ans[8] = m[2]*m2[2]+m[5]*m2[5]+m[8]*m2[8];
};

/* ----------------------------------------------------------------------
   row vector times matrix
------------------------------------------------------------------------- */

ucl_inline void gpu_row_times3(const numtyp *v, const numtyp m[9], numtyp *ans)
{
  ans[0] = m[0]*v[0]+v[1]*m[3]+v[2]*m[6];
  ans[1] = v[0]*m[1]+m[4]*v[1]+v[2]*m[7];
  ans[2] = v[0]*m[2]+v[1]*m[5]+m[8]*v[2];
};

/* ----------------------------------------------------------------------
   solve Ax = b or M ans = v
   use gaussian elimination & partial pivoting on matrix
   error_flag set to 2 if bad matrix inversion attempted
------------------------------------------------------------------------- */

ucl_inline void gpu_mldivide3(const numtyp m[9], const numtyp *v, numtyp *ans,
                            __global int *error_flag)
{
  // create augmented matrix for pivoting

  numtyp aug[12], t;

  aug[3] = v[0];
  aug[0] = m[0];
  aug[1] = m[1];
  aug[2] = m[2];
  aug[7] = v[1];
  aug[4] = m[3];
  aug[5] = m[4];
  aug[6] = m[5];
  aug[11] = v[2];
  aug[8] = m[6];
  aug[9] = m[7];
  aug[10] = m[8];

  if (ucl_abs(aug[4]) > ucl_abs(aug[0])) {
    numtyp swapt;
    swapt=aug[0]; aug[0]=aug[4]; aug[4]=swapt;
    swapt=aug[1]; aug[1]=aug[5]; aug[5]=swapt;
    swapt=aug[2]; aug[2]=aug[6]; aug[6]=swapt;
    swapt=aug[3]; aug[3]=aug[7]; aug[7]=swapt;
  }
  if (ucl_abs(aug[8]) > ucl_abs(aug[0])) {
    numtyp swapt;
    swapt=aug[0]; aug[0]=aug[8]; aug[8]=swapt;
    swapt=aug[1]; aug[1]=aug[9]; aug[9]=swapt;
    swapt=aug[2]; aug[2]=aug[10]; aug[10]=swapt;
    swapt=aug[3]; aug[3]=aug[11]; aug[11]=swapt;
  }

  if (aug[0] != (numtyp)0.0) {
    if (0!=0) {
      numtyp swapt;
      swapt=aug[0]; aug[0]=aug[0]; aug[0]=swapt;
      swapt=aug[1]; aug[1]=aug[1]; aug[1]=swapt;
      swapt=aug[2]; aug[2]=aug[2]; aug[2]=swapt;
      swapt=aug[3]; aug[3]=aug[3]; aug[3]=swapt;
    }
  } else if (aug[4] != (numtyp)0.0) {
    if (1!=0) {
      numtyp swapt;
      swapt=aug[0]; aug[0]=aug[4]; aug[4]=swapt;
      swapt=aug[1]; aug[1]=aug[5]; aug[5]=swapt;
      swapt=aug[2]; aug[2]=aug[6]; aug[6]=swapt;
      swapt=aug[3]; aug[3]=aug[7]; aug[7]=swapt;
    }
  } else if (aug[8] != (numtyp)0.0) {
    if (2!=0) {
      numtyp swapt;
      swapt=aug[0]; aug[0]=aug[8]; aug[8]=swapt;
      swapt=aug[1]; aug[1]=aug[9]; aug[9]=swapt;
      swapt=aug[2]; aug[2]=aug[10]; aug[10]=swapt;
      swapt=aug[3]; aug[3]=aug[11]; aug[11]=swapt;
    }
  } else
    *error_flag=2;

  t = aug[4]/aug[0];
  aug[5]-=t*aug[1];
  aug[6]-=t*aug[2];
  aug[7]-=t*aug[3];
  t = aug[8]/aug[0];
  aug[9]-=t*aug[1];
  aug[10]-=t*aug[2];
  aug[11]-=t*aug[3];

  if (ucl_abs(aug[9]) > ucl_abs(aug[5])) {
    numtyp swapt;
    swapt=aug[4]; aug[4]=aug[8]; aug[8]=swapt;
    swapt=aug[5]; aug[5]=aug[9]; aug[9]=swapt;
    swapt=aug[6]; aug[6]=aug[10]; aug[10]=swapt;
    swapt=aug[7]; aug[7]=aug[11]; aug[11]=swapt;
  }

  if (aug[5] != (numtyp)0.0) {
    if (1!=1) {
      numtyp swapt;
      swapt=aug[4]; aug[4]=aug[4]; aug[4]=swapt;
      swapt=aug[5]; aug[5]=aug[5]; aug[5]=swapt;
      swapt=aug[6]; aug[6]=aug[6]; aug[6]=swapt;
      swapt=aug[7]; aug[7]=aug[7]; aug[7]=swapt;
    }
  } else if (aug[9] != (numtyp)0.0) {
    if (2!=1) {
      numtyp swapt;
      swapt=aug[4]; aug[4]=aug[8]; aug[8]=swapt;
      swapt=aug[5]; aug[5]=aug[9]; aug[9]=swapt;
      swapt=aug[6]; aug[6]=aug[10]; aug[10]=swapt;
      swapt=aug[7]; aug[7]=aug[11]; aug[11]=swapt;
    }
  }

  t = aug[9]/aug[5];
  aug[10]-=t*aug[6];
  aug[11]-=t*aug[7];
  
  if (aug[10] == (numtyp)0.0)
    *error_flag=2;

  ans[2] = aug[11]/aug[10];
  t = (numtyp)0.0;
  t += aug[6]*ans[2];
  ans[1] = (aug[7]-t) / aug[5];
  t = (numtyp)0.0;
  t += aug[1]*ans[1];
  t += aug[2]*ans[2];
  ans[0] = (aug[3]-t) / aug[0];
};

/* ----------------------------------------------------------------------
   compute rotation matrix from quaternion conjugate
   quat = [w i j k]
------------------------------------------------------------------------- */

ucl_inline void gpu_quat_to_mat_trans(__global const numtyp4 *qif, const int qi, 
                                    numtyp mat[9])
{
  numtyp4 q; fetch4(q,qi,quat_tex);
  
  numtyp w2 = q.x*q.x;
  numtyp i2 = q.y*q.y;
  numtyp j2 = q.z*q.z;
  numtyp k2 = q.w*q.w;
  numtyp twoij = (numtyp)2.0*q.y*q.z;
  numtyp twoik = (numtyp)2.0*q.y*q.w;
  numtyp twojk = (numtyp)2.0*q.z*q.w;
  numtyp twoiw = (numtyp)2.0*q.y*q.x;
  numtyp twojw = (numtyp)2.0*q.z*q.x;
  numtyp twokw = (numtyp)2.0*q.w*q.x;

  mat[0] = w2+i2-j2-k2;
  mat[3] = twoij-twokw;
  mat[6] = twojw+twoik;

  mat[1] = twoij+twokw;
  mat[4] = w2-i2+j2-k2;
  mat[7] = twojk-twoiw;
	
  mat[2] = twoik-twojw;
  mat[5] = twojk+twoiw;
  mat[8] = w2-i2-j2+k2;
};

/* ----------------------------------------------------------------------
   transposed matrix times diagonal matrix
------------------------------------------------------------------------- */

ucl_inline void gpu_transpose_times_diag3(const numtyp m[9],
                                        const numtyp4 d, numtyp ans[9])
{
  ans[0] = m[0]*d.x;
  ans[1] = m[3]*d.y;
  ans[2] = m[6]*d.z;
  ans[3] = m[1]*d.x;
  ans[4] = m[4]*d.y;
  ans[5] = m[7]*d.z;
  ans[6] = m[2]*d.x;
  ans[7] = m[5]*d.y;
  ans[8] = m[8]*d.z;
};

/* ----------------------------------------------------------------------
   multiply mat1 times mat2
------------------------------------------------------------------------- */

ucl_inline void gpu_times3(const numtyp m[9], const numtyp m2[9],
                         numtyp ans[9])
{
  ans[0] = m[0]*m2[0] + m[1]*m2[3] + m[2]*m2[6];
  ans[1] = m[0]*m2[1] + m[1]*m2[4] + m[2]*m2[7];
  ans[2] = m[0]*m2[2] + m[1]*m2[5] + m[2]*m2[8];
  ans[3] = m[3]*m2[0] + m[4]*m2[3] + m[5]*m2[6];
  ans[4] = m[3]*m2[1] + m[4]*m2[4] + m[5]*m2[7];
  ans[5] = m[3]*m2[2] + m[4]*m2[5] + m[5]*m2[8];
  ans[6] = m[6]*m2[0] + m[7]*m2[3] + m[8]*m2[6];
  ans[7] = m[6]*m2[1] + m[7]*m2[4] + m[8]*m2[7];
  ans[8] = m[6]*m2[2] + m[7]*m2[5] + m[8]*m2[8];
};

/* ----------------------------------------------------------------------
   Apply principal rotation generator about x to rotation matrix m
------------------------------------------------------------------------- */

ucl_inline void gpu_rotation_generator_x(const numtyp m[9], numtyp ans[9])
{
  ans[0] = 0;
  ans[1] = -m[2];
  ans[2] = m[1];
  ans[3] = 0;
  ans[4] = -m[5];
  ans[5] = m[4];
  ans[6] = 0;
  ans[7] = -m[8];
  ans[8] = m[7];
};

/* ----------------------------------------------------------------------
   Apply principal rotation generator about y to rotation matrix m
------------------------------------------------------------------------- */

ucl_inline void gpu_rotation_generator_y(const numtyp m[9], numtyp ans[9])
{
  ans[0] = m[2];
  ans[1] = 0;
  ans[2] = -m[0];
  ans[3] = m[5];
  ans[4] = 0;
  ans[5] = -m[3];
  ans[6] = m[8];
  ans[7] = 0;
  ans[8] = -m[6];
};

/* ----------------------------------------------------------------------
   Apply principal rotation generator about z to rotation matrix m
------------------------------------------------------------------------- */

ucl_inline void gpu_rotation_generator_z(const numtyp m[9], numtyp ans[9])
{
  ans[0] = -m[1];
  ans[1] = m[0];
  ans[2] = 0;
  ans[3] = -m[4];
  ans[4] = m[3];
  ans[5] = 0;
  ans[6] = -m[7];
  ans[7] = m[6];
  ans[8] = 0;
};

/* ----------------------------------------------------------------------
   matrix times vector
------------------------------------------------------------------------- */

ucl_inline void gpu_times_column3(const numtyp m[9], const numtyp v[3],
                                numtyp ans[3]) 
{
  ans[0] = m[0]*v[0] + m[1]*v[1] + m[2]*v[2];
  ans[1] = m[3]*v[0] + m[4]*v[1] + m[5]*v[2];
  ans[2] = m[6]*v[0] + m[7]*v[1] + m[8]*v[2];
};

#endif
