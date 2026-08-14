// clang-format off
/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_MATH_EXTRA_KOKKOS_H
#define LMP_MATH_EXTRA_KOKKOS_H

#include <cmath>

#include "kokkos_type.h"

// NOTE: 'double' is still used in various quaternion related functions below.
// This is temporary to support current atom_vec_ellipsoid_kokkos bonus struct
// which still uses double for shape and quat and doesn't (yet) support KK_FLOAT.

namespace MathExtraKokkos {

  // 3 vector operations

  KOKKOS_INLINE_FUNCTION void norm3(KK_FLOAT *v);
  KOKKOS_INLINE_FUNCTION void normalize3(const KK_FLOAT *v, KK_FLOAT *ans);
  KOKKOS_INLINE_FUNCTION void snormalize3(const KK_FLOAT, const KK_FLOAT *v, KK_FLOAT *ans);
  KOKKOS_INLINE_FUNCTION void negate3(KK_FLOAT *v);
  KOKKOS_INLINE_FUNCTION void scale3(KK_FLOAT s, KK_FLOAT *v);
  KOKKOS_INLINE_FUNCTION void scale3(KK_FLOAT s, const KK_FLOAT *v, KK_FLOAT *ans);
  KOKKOS_INLINE_FUNCTION void add3(const KK_FLOAT *v1, const KK_FLOAT *v2, KK_FLOAT *ans);
  KOKKOS_INLINE_FUNCTION void sub3(const KK_FLOAT *v1, const KK_FLOAT *v2, KK_FLOAT *ans);
  KOKKOS_INLINE_FUNCTION KK_FLOAT len3(const KK_FLOAT *v);
  KOKKOS_INLINE_FUNCTION KK_FLOAT lensq3(const KK_FLOAT *v);
  KOKKOS_INLINE_FUNCTION KK_FLOAT dot3(const KK_FLOAT *v1, const KK_FLOAT *v2);
  KOKKOS_INLINE_FUNCTION void cross3(const KK_FLOAT *v1, const KK_FLOAT *v2, KK_FLOAT *ans);

  // 3x3 matrix operations

  KOKKOS_INLINE_FUNCTION KK_FLOAT det3(const KK_FLOAT mat[3][3]);
  KOKKOS_INLINE_FUNCTION void diag_times3(const KK_FLOAT *diagonal, const KK_FLOAT mat[3][3],
                          KK_FLOAT ans[3][3]);
  KOKKOS_INLINE_FUNCTION void plus3(const KK_FLOAT m[3][3], const KK_FLOAT m2[3][3],
                    KK_FLOAT ans[3][3]);
  KOKKOS_INLINE_FUNCTION void times3(const KK_FLOAT m[3][3], const KK_FLOAT m2[3][3],
                     KK_FLOAT ans[3][3]);
  KOKKOS_INLINE_FUNCTION void transpose_times3(const KK_FLOAT mat1[3][3],
                               const KK_FLOAT mat2[3][3],
                               KK_FLOAT ans[3][3]);
  KOKKOS_INLINE_FUNCTION void times3_transpose(const KK_FLOAT mat1[3][3],
                               const KK_FLOAT mat2[3][3],
                               KK_FLOAT ans[3][3]);
  KOKKOS_INLINE_FUNCTION void invert3(const KK_FLOAT mat[3][3], KK_FLOAT ans[3][3]);
  KOKKOS_INLINE_FUNCTION void matvec(const KK_FLOAT mat[3][3], const KK_FLOAT*vec, KK_FLOAT *ans);
  KOKKOS_INLINE_FUNCTION void matvec(const KK_FLOAT *ex, const KK_FLOAT *ey, const KK_FLOAT *ez,
                     const KK_FLOAT *vec, KK_FLOAT *ans);
  KOKKOS_INLINE_FUNCTION void transpose_matvec(const KK_FLOAT mat[3][3], const KK_FLOAT*vec,
                               KK_FLOAT *ans);
  KOKKOS_INLINE_FUNCTION void transpose_matvec(const KK_FLOAT *ex, const KK_FLOAT *ey,
                               const KK_FLOAT *ez, const KK_FLOAT *v,
                               KK_FLOAT *ans);
  KOKKOS_INLINE_FUNCTION void transpose_diag3(const KK_FLOAT mat[3][3], const KK_FLOAT*vec,
                              KK_FLOAT ans[3][3]);
  KOKKOS_INLINE_FUNCTION void vecmat(const KK_FLOAT *v, const KK_FLOAT m[3][3], KK_FLOAT *ans);
  KOKKOS_INLINE_FUNCTION void scalar_times3(const KK_FLOAT f, KK_FLOAT m[3][3]);

  KOKKOS_INLINE_FUNCTION void richardson(double *q, KK_FLOAT *m, KK_FLOAT *w, KK_FLOAT *moments, KK_FLOAT dtq);

  // quaternion operations
  KOKKOS_INLINE_FUNCTION void qnormalize(double *q);
  KOKKOS_INLINE_FUNCTION void qconjugate(KK_FLOAT *q, KK_FLOAT *qc);
  KOKKOS_INLINE_FUNCTION void vecquat(KK_FLOAT *a, double *b, KK_FLOAT *c);
  KOKKOS_INLINE_FUNCTION void axisangle_to_quat(const KK_FLOAT *v, const KK_FLOAT angle,
                                KK_FLOAT *quat);

  KOKKOS_INLINE_FUNCTION void mq_to_omega(KK_FLOAT *m, double *q, KK_FLOAT *moments, KK_FLOAT *w);
  KOKKOS_INLINE_FUNCTION void quat_to_mat(const double *quat, KK_FLOAT mat[3][3]);
}

/* ----------------------------------------------------------------------
   normalize a vector in place
------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::norm3(KK_FLOAT *v)
{
  KK_FLOAT scale = 1.0/sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  v[0] *= scale;
  v[1] *= scale;
  v[2] *= scale;
}

/* ----------------------------------------------------------------------
   normalize a vector, return in ans
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::normalize3(const KK_FLOAT *v, KK_FLOAT *ans)
{
  KK_FLOAT scale = 1.0/sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  ans[0] = v[0]*scale;
  ans[1] = v[1]*scale;
  ans[2] = v[2]*scale;
}

/* ----------------------------------------------------------------------
   scale a vector to length
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::snormalize3(const KK_FLOAT length, const KK_FLOAT *v, KK_FLOAT *ans)
{
  KK_FLOAT scale = length/sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  ans[0] = v[0]*scale;
  ans[1] = v[1]*scale;
  ans[2] = v[2]*scale;
}

/* ----------------------------------------------------------------------
   negate vector v in place
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::negate3(KK_FLOAT *v)
{
  v[0] = -v[0];
  v[1] = -v[1];
  v[2] = -v[2];
}

/* ----------------------------------------------------------------------
   scale vector v by s in place
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::scale3(KK_FLOAT s, KK_FLOAT *v)
{
  v[0] *= s;
  v[1] *= s;
  v[2] *= s;
}

/* ----------------------------------------------------------------------
   scale vector v by s, return in ans
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::scale3(KK_FLOAT s, const KK_FLOAT *v, KK_FLOAT *ans)
{
  ans[0] = s*v[0];
  ans[1] = s*v[1];
  ans[2] = s*v[2];
}

/* ----------------------------------------------------------------------
   ans = v1 + v2
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::add3(const KK_FLOAT *v1, const KK_FLOAT *v2, KK_FLOAT *ans)
{
  ans[0] = v1[0] + v2[0];
  ans[1] = v1[1] + v2[1];
  ans[2] = v1[2] + v2[2];
}

/* ----------------------------------------------------------------------
   ans = v1 - v2
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::sub3(const KK_FLOAT *v1, const KK_FLOAT *v2, KK_FLOAT *ans)
{
  ans[0] = v1[0] - v2[0];
  ans[1] = v1[1] - v2[1];
  ans[2] = v1[2] - v2[2];
}

/* ----------------------------------------------------------------------
   length of vector v
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
KK_FLOAT MathExtraKokkos::len3(const KK_FLOAT *v)
{
  return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

/* ----------------------------------------------------------------------
   squared length of vector v, or dot product of v with itself
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
KK_FLOAT MathExtraKokkos::lensq3(const KK_FLOAT *v)
{
  return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}

/* ----------------------------------------------------------------------
   dot product of 2 vectors
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
KK_FLOAT MathExtraKokkos::dot3(const KK_FLOAT *v1, const KK_FLOAT *v2)
{
  return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

/* ----------------------------------------------------------------------
   cross product of 2 vectors
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::cross3(const KK_FLOAT *v1, const KK_FLOAT *v2, KK_FLOAT *ans)
{
  ans[0] = v1[1]*v2[2] - v1[2]*v2[1];
  ans[1] = v1[2]*v2[0] - v1[0]*v2[2];
  ans[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

/* ----------------------------------------------------------------------
   determinant of a matrix
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
KK_FLOAT MathExtraKokkos::det3(const KK_FLOAT m[3][3])
{
  KK_FLOAT ans = m[0][0]*m[1][1]*m[2][2] - m[0][0]*m[1][2]*m[2][1] -
    m[1][0]*m[0][1]*m[2][2] + m[1][0]*m[0][2]*m[2][1] +
    m[2][0]*m[0][1]*m[1][2] - m[2][0]*m[0][2]*m[1][1];
  return ans;
}

/* ----------------------------------------------------------------------
   diagonal matrix times a full matrix
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::diag_times3(const KK_FLOAT *d, const KK_FLOAT m[3][3],
                            KK_FLOAT ans[3][3])
{
  ans[0][0] = d[0]*m[0][0];
  ans[0][1] = d[0]*m[0][1];
  ans[0][2] = d[0]*m[0][2];
  ans[1][0] = d[1]*m[1][0];
  ans[1][1] = d[1]*m[1][1];
  ans[1][2] = d[1]*m[1][2];
  ans[2][0] = d[2]*m[2][0];
  ans[2][1] = d[2]*m[2][1];
  ans[2][2] = d[2]*m[2][2];
}

/* ----------------------------------------------------------------------
   add two matrices
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::plus3(const KK_FLOAT m[3][3], const KK_FLOAT m2[3][3],
                      KK_FLOAT ans[3][3])
{
  ans[0][0] = m[0][0]+m2[0][0];
  ans[0][1] = m[0][1]+m2[0][1];
  ans[0][2] = m[0][2]+m2[0][2];
  ans[1][0] = m[1][0]+m2[1][0];
  ans[1][1] = m[1][1]+m2[1][1];
  ans[1][2] = m[1][2]+m2[1][2];
  ans[2][0] = m[2][0]+m2[2][0];
  ans[2][1] = m[2][1]+m2[2][1];
  ans[2][2] = m[2][2]+m2[2][2];
}

/* ----------------------------------------------------------------------
   multiply mat1 times mat2
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::times3(const KK_FLOAT m[3][3], const KK_FLOAT m2[3][3],
                       KK_FLOAT ans[3][3])
{
  ans[0][0] = m[0][0]*m2[0][0] + m[0][1]*m2[1][0] + m[0][2]*m2[2][0];
  ans[0][1] = m[0][0]*m2[0][1] + m[0][1]*m2[1][1] + m[0][2]*m2[2][1];
  ans[0][2] = m[0][0]*m2[0][2] + m[0][1]*m2[1][2] + m[0][2]*m2[2][2];
  ans[1][0] = m[1][0]*m2[0][0] + m[1][1]*m2[1][0] + m[1][2]*m2[2][0];
  ans[1][1] = m[1][0]*m2[0][1] + m[1][1]*m2[1][1] + m[1][2]*m2[2][1];
  ans[1][2] = m[1][0]*m2[0][2] + m[1][1]*m2[1][2] + m[1][2]*m2[2][2];
  ans[2][0] = m[2][0]*m2[0][0] + m[2][1]*m2[1][0] + m[2][2]*m2[2][0];
  ans[2][1] = m[2][0]*m2[0][1] + m[2][1]*m2[1][1] + m[2][2]*m2[2][1];
  ans[2][2] = m[2][0]*m2[0][2] + m[2][1]*m2[1][2] + m[2][2]*m2[2][2];
}

/* ----------------------------------------------------------------------
   multiply the transpose of mat1 times mat2
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::transpose_times3(const KK_FLOAT m[3][3], const KK_FLOAT m2[3][3],
                                 KK_FLOAT ans[3][3])
{
  ans[0][0] = m[0][0]*m2[0][0] + m[1][0]*m2[1][0] + m[2][0]*m2[2][0];
  ans[0][1] = m[0][0]*m2[0][1] + m[1][0]*m2[1][1] + m[2][0]*m2[2][1];
  ans[0][2] = m[0][0]*m2[0][2] + m[1][0]*m2[1][2] + m[2][0]*m2[2][2];
  ans[1][0] = m[0][1]*m2[0][0] + m[1][1]*m2[1][0] + m[2][1]*m2[2][0];
  ans[1][1] = m[0][1]*m2[0][1] + m[1][1]*m2[1][1] + m[2][1]*m2[2][1];
  ans[1][2] = m[0][1]*m2[0][2] + m[1][1]*m2[1][2] + m[2][1]*m2[2][2];
  ans[2][0] = m[0][2]*m2[0][0] + m[1][2]*m2[1][0] + m[2][2]*m2[2][0];
  ans[2][1] = m[0][2]*m2[0][1] + m[1][2]*m2[1][1] + m[2][2]*m2[2][1];
  ans[2][2] = m[0][2]*m2[0][2] + m[1][2]*m2[1][2] + m[2][2]*m2[2][2];
}

/* ----------------------------------------------------------------------
   multiply mat1 times transpose of mat2
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::times3_transpose(const KK_FLOAT m[3][3], const KK_FLOAT m2[3][3],
                                 KK_FLOAT ans[3][3])
{
  ans[0][0] = m[0][0]*m2[0][0] + m[0][1]*m2[0][1] + m[0][2]*m2[0][2];
  ans[0][1] = m[0][0]*m2[1][0] + m[0][1]*m2[1][1] + m[0][2]*m2[1][2];
  ans[0][2] = m[0][0]*m2[2][0] + m[0][1]*m2[2][1] + m[0][2]*m2[2][2];
  ans[1][0] = m[1][0]*m2[0][0] + m[1][1]*m2[0][1] + m[1][2]*m2[0][2];
  ans[1][1] = m[1][0]*m2[1][0] + m[1][1]*m2[1][1] + m[1][2]*m2[1][2];
  ans[1][2] = m[1][0]*m2[2][0] + m[1][1]*m2[2][1] + m[1][2]*m2[2][2];
  ans[2][0] = m[2][0]*m2[0][0] + m[2][1]*m2[0][1] + m[2][2]*m2[0][2];
  ans[2][1] = m[2][0]*m2[1][0] + m[2][1]*m2[1][1] + m[2][2]*m2[1][2];
  ans[2][2] = m[2][0]*m2[2][0] + m[2][1]*m2[2][1] + m[2][2]*m2[2][2];
}

/* ----------------------------------------------------------------------
   invert a matrix
   does NOT checks for singular or badly scaled matrix
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::invert3(const KK_FLOAT m[3][3], KK_FLOAT ans[3][3])
{
  KK_FLOAT den = m[0][0]*m[1][1]*m[2][2]-m[0][0]*m[1][2]*m[2][1];
  den += -m[1][0]*m[0][1]*m[2][2]+m[1][0]*m[0][2]*m[2][1];
  den += m[2][0]*m[0][1]*m[1][2]-m[2][0]*m[0][2]*m[1][1];

  ans[0][0] = (m[1][1]*m[2][2]-m[1][2]*m[2][1]) / den;
  ans[0][1] = -(m[0][1]*m[2][2]-m[0][2]*m[2][1]) / den;
  ans[0][2] = (m[0][1]*m[1][2]-m[0][2]*m[1][1]) / den;
  ans[1][0] = -(m[1][0]*m[2][2]-m[1][2]*m[2][0]) / den;
  ans[1][1] = (m[0][0]*m[2][2]-m[0][2]*m[2][0]) / den;
  ans[1][2] = -(m[0][0]*m[1][2]-m[0][2]*m[1][0]) / den;
  ans[2][0] = (m[1][0]*m[2][1]-m[1][1]*m[2][0]) / den;
  ans[2][1] = -(m[0][0]*m[2][1]-m[0][1]*m[2][0]) / den;
  ans[2][2] = (m[0][0]*m[1][1]-m[0][1]*m[1][0]) / den;
}

/* ----------------------------------------------------------------------
   matrix times vector
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::matvec(const KK_FLOAT m[3][3], const KK_FLOAT *v, KK_FLOAT *ans)
{
  ans[0] = m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2];
  ans[1] = m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2];
  ans[2] = m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2];
}

/* ----------------------------------------------------------------------
   matrix times vector
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::matvec(const KK_FLOAT *ex, const KK_FLOAT *ey, const KK_FLOAT *ez,
                       const KK_FLOAT *v, KK_FLOAT *ans)
{
  ans[0] = ex[0]*v[0] + ey[0]*v[1] + ez[0]*v[2];
  ans[1] = ex[1]*v[0] + ey[1]*v[1] + ez[1]*v[2];
  ans[2] = ex[2]*v[0] + ey[2]*v[1] + ez[2]*v[2];
}

/* ----------------------------------------------------------------------
   transposed matrix times vector
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::transpose_matvec(const KK_FLOAT m[3][3], const KK_FLOAT *v,
                                 KK_FLOAT *ans)
{
  ans[0] = m[0][0]*v[0] + m[1][0]*v[1] + m[2][0]*v[2];
  ans[1] = m[0][1]*v[0] + m[1][1]*v[1] + m[2][1]*v[2];
  ans[2] = m[0][2]*v[0] + m[1][2]*v[1] + m[2][2]*v[2];
}

/* ----------------------------------------------------------------------
   transposed matrix times vector
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::transpose_matvec(const KK_FLOAT *ex, const KK_FLOAT *ey,
                                 const KK_FLOAT *ez, const KK_FLOAT *v,
                                 KK_FLOAT *ans)
{
  ans[0] = ex[0]*v[0] + ex[1]*v[1] + ex[2]*v[2];
  ans[1] = ey[0]*v[0] + ey[1]*v[1] + ey[2]*v[2];
  ans[2] = ez[0]*v[0] + ez[1]*v[1] + ez[2]*v[2];
}

/* ----------------------------------------------------------------------
   transposed matrix times diagonal matrix
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::transpose_diag3(const KK_FLOAT m[3][3], const KK_FLOAT *d,
                                KK_FLOAT ans[3][3])
{
  ans[0][0] = m[0][0]*d[0];
  ans[0][1] = m[1][0]*d[1];
  ans[0][2] = m[2][0]*d[2];
  ans[1][0] = m[0][1]*d[0];
  ans[1][1] = m[1][1]*d[1];
  ans[1][2] = m[2][1]*d[2];
  ans[2][0] = m[0][2]*d[0];
  ans[2][1] = m[1][2]*d[1];
  ans[2][2] = m[2][2]*d[2];
}

/* ----------------------------------------------------------------------
   row vector times matrix
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::vecmat(const KK_FLOAT *v, const KK_FLOAT m[3][3], KK_FLOAT *ans)
{
  ans[0] = v[0]*m[0][0] + v[1]*m[1][0] + v[2]*m[2][0];
  ans[1] = v[0]*m[0][1] + v[1]*m[1][1] + v[2]*m[2][1];
  ans[2] = v[0]*m[0][2] + v[1]*m[1][2] + v[2]*m[2][2];
}

/* ----------------------------------------------------------------------
   matrix times scalar, in place
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::scalar_times3(const KK_FLOAT f, KK_FLOAT m[3][3])
{
  m[0][0] *= f; m[0][1] *= f; m[0][2] *= f;
  m[1][0] *= f; m[1][1] *= f; m[1][2] *= f;
  m[2][0] *= f; m[2][1] *= f; m[2][2] *= f;
}

/* ----------------------------------------------------------------------
   Richardson iteration to update quaternion from angular momentum
   return new normalized quaternion q
   also returns updated omega at 1/2 step
------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::richardson(double *q, KK_FLOAT *m, KK_FLOAT *w, KK_FLOAT *moments, KK_FLOAT dtq)
{
  // full update from dq/dt = 1/2 w q

  KK_FLOAT wq[4];
  MathExtraKokkos::vecquat(w,q,wq);

  double qfull[4];
  qfull[0] = q[0] + dtq * wq[0];
  qfull[1] = q[1] + dtq * wq[1];
  qfull[2] = q[2] + dtq * wq[2];
  qfull[3] = q[3] + dtq * wq[3];
  MathExtraKokkos::qnormalize(qfull);

  // 1st half update from dq/dt = 1/2 w q

  double qhalf[4];
  qhalf[0] = q[0] + 0.5*dtq * wq[0];
  qhalf[1] = q[1] + 0.5*dtq * wq[1];
  qhalf[2] = q[2] + 0.5*dtq * wq[2];
  qhalf[3] = q[3] + 0.5*dtq * wq[3];
  MathExtraKokkos::qnormalize(qhalf);

  // re-compute omega at 1/2 step from m at 1/2 step and q at 1/2 step
  // recompute wq

  MathExtraKokkos::mq_to_omega(m,qhalf,moments,w);
  MathExtraKokkos::vecquat(w,qhalf,wq);

  // 2nd half update from dq/dt = 1/2 w q

  qhalf[0] += 0.5*dtq * wq[0];
  qhalf[1] += 0.5*dtq * wq[1];
  qhalf[2] += 0.5*dtq * wq[2];
  qhalf[3] += 0.5*dtq * wq[3];
  MathExtraKokkos::qnormalize(qhalf);

  // corrected Richardson update

  q[0] = 2.0*qhalf[0] - qfull[0];
  q[1] = 2.0*qhalf[1] - qfull[1];
  q[2] = 2.0*qhalf[2] - qfull[2];
  q[3] = 2.0*qhalf[3] - qfull[3];
  MathExtraKokkos::qnormalize(q);
}

/* ----------------------------------------------------------------------
   normalize a quaternion
------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::qnormalize(double *q)
{
  KK_FLOAT norm = 1.0 / sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
  q[0] *= norm;
  q[1] *= norm;
  q[2] *= norm;
  q[3] *= norm;
}

/* ----------------------------------------------------------------------
   conjugate of a quaternion: qc = conjugate of q
   assume q is of unit length
------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::qconjugate(KK_FLOAT *q, KK_FLOAT *qc)
{
  qc[0] = q[0];
  qc[1] = -q[1];
  qc[2] = -q[2];
  qc[3] = -q[3];
}

/* ----------------------------------------------------------------------
   vector-quaternion multiply: c = a*b, where a = (0,a)
------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::vecquat(KK_FLOAT *a, double *b, KK_FLOAT *c)
{
  c[0] = -a[0] * b[1] - a[1] * b[2] - a[2] * b[3];
  c[1] = b[0] * a[0] + a[1] * b[3] - a[2] * b[2];
  c[2] = b[0] * a[1] + a[2] * b[1] - a[0] * b[3];
  c[3] = b[0] * a[2] + a[0] * b[2] - a[1] * b[1];
}

/* ----------------------------------------------------------------------
   compute quaternion from axis-angle rotation
   v MUST be a unit vector
------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::axisangle_to_quat(const KK_FLOAT *v, const KK_FLOAT angle,
                                  KK_FLOAT *quat)
{
  KK_FLOAT halfa = 0.5*angle;
  KK_FLOAT sina = sin(halfa);
  quat[0] = cos(halfa);
  quat[1] = v[0]*sina;
  quat[2] = v[1]*sina;
  quat[3] = v[2]*sina;
}

/* ----------------------------------------------------------------------
   compute omega from angular momentum
   w = omega = angular velocity in space frame
   wbody = angular velocity in body frame
   project space-frame angular momentum onto body axes
     and divide by principal moments
------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::mq_to_omega(KK_FLOAT *m, double *q, KK_FLOAT *moments, KK_FLOAT *w)
{
  KK_FLOAT wbody[3];
  KK_FLOAT rot[3][3];

  MathExtraKokkos::quat_to_mat(q,rot);
  MathExtraKokkos::transpose_matvec(rot,m,wbody);
  if (moments[0] == 0.0) wbody[0] = 0.0;
  else wbody[0] /= moments[0];
  if (moments[1] == 0.0) wbody[1] = 0.0;
  else wbody[1] /= moments[1];
  if (moments[2] == 0.0) wbody[2] = 0.0;
  else wbody[2] /= moments[2];
  MathExtraKokkos::matvec(rot,wbody,w);
}

/* ----------------------------------------------------------------------
   compute rotation matrix from quaternion
   quat = [w i j k]
------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::quat_to_mat(const double *quat, KK_FLOAT mat[3][3])
{
  KK_FLOAT w2 = quat[0]*quat[0];
  KK_FLOAT i2 = quat[1]*quat[1];
  KK_FLOAT j2 = quat[2]*quat[2];
  KK_FLOAT k2 = quat[3]*quat[3];
  KK_FLOAT twoij = 2.0*quat[1]*quat[2];
  KK_FLOAT twoik = 2.0*quat[1]*quat[3];
  KK_FLOAT twojk = 2.0*quat[2]*quat[3];
  KK_FLOAT twoiw = 2.0*quat[1]*quat[0];
  KK_FLOAT twojw = 2.0*quat[2]*quat[0];
  KK_FLOAT twokw = 2.0*quat[3]*quat[0];

  mat[0][0] = w2+i2-j2-k2;
  mat[0][1] = twoij-twokw;
  mat[0][2] = twojw+twoik;

  mat[1][0] = twoij+twokw;
  mat[1][1] = w2-i2+j2-k2;
  mat[1][2] = twojk-twoiw;

  mat[2][0] = twoik-twojw;
  mat[2][1] = twojk+twoiw;
  mat[2][2] = w2-i2-j2+k2;
}

#endif
