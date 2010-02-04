/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors: Mike Brown (SNL), wmbrown@sandia.gov
                         Peng Wang (Nvidia), penwang@nvidia.com
                         Paul Crozier (SNL), pscrozi@sandia.gov
------------------------------------------------------------------------- */

#ifndef GB_GPU_EXTRA_H
#define GB_GPU_EXTRA_H

#include "math.h"
#include "stdio.h"
#include "string.h"

/* ----------------------------------------------------------------------
   Atomic update of global memory
------------------------------------------------------------------------- */
/*
template <class numtyp> __device__ 
inline void atomicAdd(numtyp *address, numtyp val);

template <>
__device__ inline void atomicAdd<float>(float *address, float val)
{
  int i_val = __float_as_int(val);
  int tmp0 = 0;
  int tmp1;

  while( (tmp1 = atomicCAS((int *)address, tmp0, i_val)) != tmp0) {
    tmp0 = tmp1;
    i_val = __float_as_int(val + __int_as_float(tmp1));
  }
}*/

/* ----------------------------------------------------------------------
   dot product of 2 vectors
------------------------------------------------------------------------- */

template <class numtyp>
static __inline__ __device__ numtyp gpu_dot3(const numtyp *v1, const numtyp *v2)
{
  return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

/* ----------------------------------------------------------------------
   cross product of 2 vectors
------------------------------------------------------------------------- */

template <class numtyp>
static __inline__ __device__ void gpu_cross3(const numtyp *v1, 
                                             const numtyp *v2, numtyp *ans)
{
  ans[0] = v1[1]*v2[2]-v1[2]*v2[1];
  ans[1] = v1[2]*v2[0]-v1[0]*v2[2];
  ans[2] = v1[0]*v2[1]-v1[1]*v2[0];
}

/* ----------------------------------------------------------------------
   determinant of a matrix
------------------------------------------------------------------------- */

template <class numtyp>
static __inline__ __device__ numtyp gpu_det3(const numtyp m[9])
{
  numtyp ans = m[0]*m[4]*m[8] - m[0]*m[5]*m[7] - 
    m[3]*m[1]*m[8] + m[3]*m[2]*m[7] + 
    m[6]*m[1]*m[5] - m[6]*m[2]*m[4];
  return ans;
}

/* ----------------------------------------------------------------------
   diagonal matrix times a full matrix
------------------------------------------------------------------------- */

template <class numtyp>
static __inline__ __device__ void gpu_well_times3(const int i, const numtyp m[9],
                                                  numtyp ans[9])
{
  ans[0] = _well_<numtyp>(i,0)*m[0];
  ans[1] = _well_<numtyp>(i,0)*m[1];
  ans[2] = _well_<numtyp>(i,0)*m[2];
  ans[3] = _well_<numtyp>(i,1)*m[3];
  ans[4] = _well_<numtyp>(i,1)*m[4];
  ans[5] = _well_<numtyp>(i,1)*m[5];
  ans[6] = _well_<numtyp>(i,2)*m[6];
  ans[7] = _well_<numtyp>(i,2)*m[7];
  ans[8] = _well_<numtyp>(i,2)*m[8];
}

/* ----------------------------------------------------------------------
   diagonal matrix times a full matrix
------------------------------------------------------------------------- */

template <class numtyp>
static __inline__ __device__ void gpu_shape_times3(const int i, const numtyp m[9],
                                                   numtyp ans[9])
{
  ans[0] = _shape_<numtyp>(i,0)*m[0];
  ans[1] = _shape_<numtyp>(i,0)*m[1];
  ans[2] = _shape_<numtyp>(i,0)*m[2];
  ans[3] = _shape_<numtyp>(i,1)*m[3];
  ans[4] = _shape_<numtyp>(i,1)*m[4];
  ans[5] = _shape_<numtyp>(i,1)*m[5];
  ans[6] = _shape_<numtyp>(i,2)*m[6];
  ans[7] = _shape_<numtyp>(i,2)*m[7];
  ans[8] = _shape_<numtyp>(i,2)*m[8];
}

/* ----------------------------------------------------------------------
   add two matrices
------------------------------------------------------------------------- */

template <class numtyp>
static __inline__ __device__ void gpu_plus3(const numtyp m[9], 
                                            const numtyp m2[9], numtyp ans[9])
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
}

/* ----------------------------------------------------------------------
   multiply the transpose of mat1 times mat2
------------------------------------------------------------------------- */

template <class numtyp>
static __inline__ __device__ void gpu_transpose_times3(const numtyp m[9], 
                                                       const numtyp m2[9],
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
}

/* ----------------------------------------------------------------------
   row vector times matrix
------------------------------------------------------------------------- */

template <class numtyp>
static __inline__ __device__ void gpu_row_times3(const numtyp *v, 
                                                 const numtyp m[9], numtyp *ans)
{
  ans[0] = m[0]*v[0]+v[1]*m[3]+v[2]*m[6];
  ans[1] = v[0]*m[1]+m[4]*v[1]+v[2]*m[7];
  ans[2] = v[0]*m[2]+v[1]*m[5]+m[8]*v[2];
}

/* ----------------------------------------------------------------------
   solve Ax = b or M ans = v
   use gaussian elimination & partial pivoting on matrix
   error_flag set to 2 if bad matrix inversion attempted
------------------------------------------------------------------------- */

template <class numtyp>
static __inline__ __device__ void gpu_mldivide3(const numtyp m[9], 
                                                const numtyp *v, numtyp *ans,
                                                int *error_flag)
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

  if (fabs(aug[4]) > fabs(aug[0])) {
    numtyp swapt;
    swapt=aug[0]; aug[0]=aug[4]; aug[4]=swapt;
    swapt=aug[1]; aug[1]=aug[5]; aug[5]=swapt;
    swapt=aug[2]; aug[2]=aug[6]; aug[6]=swapt;
    swapt=aug[3]; aug[3]=aug[7]; aug[7]=swapt;
  }
  if (fabs(aug[8]) > fabs(aug[0])) {
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

  if (fabs(aug[9]) > fabs(aug[5])) {
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
}

/* ----------------------------------------------------------------------
   compute rotation matrix from quaternion conjugate
   quat = [w i j k]
------------------------------------------------------------------------- */

template <class numtyp>
static __inline__ __device__ void gpu_quat_to_mat_trans(const int qi, 
                                                        numtyp mat[9])
{
  numtyp qi3=_x_<numtyp>(qi,3);
  numtyp qi4=_x_<numtyp>(qi,4);
  numtyp qi5=_x_<numtyp>(qi,5);
  numtyp qi6=_x_<numtyp>(qi,6);

  numtyp w2 = qi3*qi3;
  numtyp i2 = qi4*qi4;
  numtyp j2 = qi5*qi5;
  numtyp k2 = qi6*qi6;
  numtyp twoij = (numtyp)2.0*qi4*qi5;
  numtyp twoik = (numtyp)2.0*qi4*qi6;
  numtyp twojk = (numtyp)2.0*qi5*qi6;
  numtyp twoiw = (numtyp)2.0*qi4*qi3;
  numtyp twojw = (numtyp)2.0*qi5*qi3;
  numtyp twokw = (numtyp)2.0*qi6*qi3;

  mat[0] = w2+i2-j2-k2;
  mat[3] = twoij-twokw;
  mat[6] = twojw+twoik;

  mat[1] = twoij+twokw;
  mat[4] = w2-i2+j2-k2;
  mat[7] = twojk-twoiw;
	
  mat[2] = twoik-twojw;
  mat[5] = twojk+twoiw;
  mat[8] = w2-i2-j2+k2;
}

#endif
