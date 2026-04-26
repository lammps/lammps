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

#include "kokkos_type.h"

// NOTE: 'double' is still used in various quaternion related functions below.
// This is temporary to support current atom_vec_ellipsoid_kokkos bonus struct
// which still uses double for shape and quat and doesn't (yet) support KK_FLOAT.
//
// UPDATE [2026/04 alphataubio]:
// functions now templated for T = KK_FLOAT / float / double / ...

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

  template <typename T>
  KOKKOS_INLINE_FUNCTION void matvec(const T mat[3][3], const T*vec, T *ans);

  template <typename T>
  KOKKOS_INLINE_FUNCTION void matvec(const T *ex, const T *ey, const T *ez,
                     const T *vec, T *ans);

  template <typename T>
  KOKKOS_INLINE_FUNCTION void transpose_matvec(const T mat[3][3], const T*vec,
                               T *ans);

  template <typename T>
  KOKKOS_INLINE_FUNCTION void transpose_matvec(const T *ex, const T *ey,
                               const T *ez, const T *v, T *ans);

  KOKKOS_INLINE_FUNCTION void transpose_diag3(const KK_FLOAT mat[3][3], const KK_FLOAT*vec,
                              KK_FLOAT ans[3][3]);
  KOKKOS_INLINE_FUNCTION void vecmat(const KK_FLOAT *v, const KK_FLOAT m[3][3], KK_FLOAT *ans);
  KOKKOS_INLINE_FUNCTION void scalar_times3(const KK_FLOAT f, KK_FLOAT m[3][3]);

  template <typename T>
  KOKKOS_INLINE_FUNCTION void richardson(T *q, KK_FLOAT *m, KK_FLOAT *w, KK_FLOAT *moments, KK_FLOAT dtq);

  // quaternion operations

  template <typename T>
  KOKKOS_INLINE_FUNCTION void qnormalize(T *q);

  KOKKOS_INLINE_FUNCTION void qconjugate(KK_FLOAT *q, KK_FLOAT *qc);

  template <typename T>
  KOKKOS_INLINE_FUNCTION void vecquat(KK_FLOAT *a, T *b, KK_FLOAT *c);

  KOKKOS_INLINE_FUNCTION void axisangle_to_quat(const KK_FLOAT *v, const KK_FLOAT angle,
                                KK_FLOAT *quat);

  template <typename T>
  KOKKOS_INLINE_FUNCTION
  void mq_to_omega(KK_FLOAT *m, T *q, KK_FLOAT *moments, KK_FLOAT *w);

  template <typename T>
  KOKKOS_INLINE_FUNCTION void quat_to_mat(const T *quat, T mat[3][3]);

  template <typename T>
  KOKKOS_INLINE_FUNCTION void exyz_to_q(const T *q, const T *ex, const T *ey, T *ez);

  template <typename T>
  KOKKOS_INLINE_FUNCTION void q_to_exyz(const T *q, T *ex, T *ey, T *ez);

  template <typename T>
  KOKKOS_INLINE_FUNCTION void quatvec(const T *a, const T *b, T *c);

  template <typename T>
  KOKKOS_INLINE_FUNCTION void invquatvec(const T *a, const KK_FLOAT *b, KK_FLOAT *c);

  template <typename T>
  KOKKOS_INLINE_FUNCTION void quatquat(const T *a, const T *b, T *c);

  template <typename T>
  KOKKOS_INLINE_FUNCTION void no_squish_rotate(int k, T *p, T *q,
                                               const KK_FLOAT *inertia, KK_FLOAT dt);

  KOKKOS_INLINE_FUNCTION void angmom_to_omega(const KK_FLOAT *m, const KK_FLOAT *ex,
                                              const KK_FLOAT *ey, const KK_FLOAT *ez,
                                              const KK_FLOAT *idiag, KK_FLOAT *w);
  KOKKOS_INLINE_FUNCTION void omega_to_angmom(const KK_FLOAT *w, const KK_FLOAT *ex,
                                              const KK_FLOAT *ey, const KK_FLOAT *ez,
                                              const KK_FLOAT *idiag, KK_FLOAT *m);

  template <typename T>
  KOKKOS_INLINE_FUNCTION
  int sym3x3_eigen(const T A[3][3], T evals[3], T evecs[3][3], int sort = 1) noexcept;

  using Kokkos::fma;

}

/* ----------------------------------------------------------------------
   normalize a vector in place
------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::norm3(KK_FLOAT *v)
{
  const KK_FLOAT scale = KK_FLOAT(1.0) / MathExtraKokkos::len3(v);
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
  const KK_FLOAT scale = KK_FLOAT(1.0) / MathExtraKokkos::len3(v);
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
  const KK_FLOAT scale = length / MathExtraKokkos::len3(v);
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
  return Kokkos::sqrt(lensq3(v));
}

/* ----------------------------------------------------------------------
   squared length of vector v, or dot product of v with itself
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
KK_FLOAT MathExtraKokkos::lensq3(const KK_FLOAT *v)
{
  // v0^2 + v1^2 + v2^2
  return fma(v[0], v[0], fma(v[1], v[1], v[2]*v[2]));
}

/* ----------------------------------------------------------------------
   dot product of 2 vectors
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
KK_FLOAT MathExtraKokkos::dot3(const KK_FLOAT *v1, const KK_FLOAT *v2)
{
  // v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
  return fma(v1[0], v2[0], fma(v1[1], v2[1], v1[2]*v2[2]));
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

template <typename T>
KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::matvec(const T m[3][3], const T *v, T *ans)
{
  ans[0] = m[0][0]*v[0] + m[0][1]*v[1] + m[0][2]*v[2];
  ans[1] = m[1][0]*v[0] + m[1][1]*v[1] + m[1][2]*v[2];
  ans[2] = m[2][0]*v[0] + m[2][1]*v[1] + m[2][2]*v[2];
}

/* ----------------------------------------------------------------------
   matrix times vector
------------------------------------------------------------------------- */

template <typename T>
KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::matvec(const T *ex, const T *ey, const T *ez,
                       const T *v, T *ans)
{
  ans[0] = fma(ez[0], v[2], fma(ey[0], v[1], ex[0]*v[0]));
  ans[1] = fma(ez[1], v[2], fma(ey[1], v[1], ex[1]*v[0]));
  ans[2] = fma(ez[2], v[2], fma(ey[2], v[1], ex[2]*v[0]));
}

/* ----------------------------------------------------------------------
   transposed matrix times vector
------------------------------------------------------------------------- */

template <typename T>
KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::transpose_matvec(const T m[3][3], const T *v, T *ans)
{
  ans[0] = fma(m[2][0], v[2], fma(m[1][0], v[1], m[0][0]*v[0]));
  ans[1] = fma(m[2][1], v[2], fma(m[1][1], v[1], m[0][1]*v[0]));
  ans[2] = fma(m[2][2], v[2], fma(m[1][2], v[1], m[0][2]*v[0]));
}

/* ----------------------------------------------------------------------
   transposed matrix times vector
------------------------------------------------------------------------- */

template <typename T>
KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::transpose_matvec(const T *ex, const T *ey,
                                       const T *ez, const T *v, T *ans)
{
  ans[0] = fma(ex[2], v[2], fma(ex[1], v[1], ex[0]*v[0]));
  ans[1] = fma(ey[2], v[2], fma(ey[1], v[1], ey[0]*v[0]));
  ans[2] = fma(ez[2], v[2], fma(ez[1], v[1], ez[0]*v[0]));
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
template <typename T>
KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::richardson(T *q, KK_FLOAT *m, KK_FLOAT *w, KK_FLOAT *moments, KK_FLOAT dtq)
{
  // full update from dq/dt = 1/2 w q

  T wq[4];
  MathExtraKokkos::vecquat(w,q,wq);

  T qfull[4];
  qfull[0] = q[0] + dtq * wq[0];
  qfull[1] = q[1] + dtq * wq[1];
  qfull[2] = q[2] + dtq * wq[2];
  qfull[3] = q[3] + dtq * wq[3];
  MathExtraKokkos::qnormalize(qfull);

  // 1st half update from dq/dt = 1/2 w q

  T qhalf[4];
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
template <typename T>
KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::qnormalize(T *q)
{
  const T norm = 1.0 / Kokkos::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
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
template <typename T>
KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::vecquat(KK_FLOAT *a, T *b, KK_FLOAT *c)
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
  KK_FLOAT sina = Kokkos::sin(halfa);
  quat[0] = Kokkos::cos(halfa);
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
template <typename T>
KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::mq_to_omega(KK_FLOAT *m, T *q, KK_FLOAT *moments, KK_FLOAT *w)
{
  KK_FLOAT wbody[3], rot[3][3], kk_q[4] = {
    static_cast<KK_FLOAT>(q[0]),
    static_cast<KK_FLOAT>(q[1]),
    static_cast<KK_FLOAT>(q[2]),
    static_cast<KK_FLOAT>(q[3])
  };
  MathExtraKokkos::quat_to_mat(kk_q,rot);
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
   create unit quaternion from space-frame ex,ey,ez
   ex,ey,ez are columns of a rotation matrix
------------------------------------------------------------------------- */

template <typename T>
KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::exyz_to_q(const T *ex, const T *ey, const T *ez, T *q)
{
  // squares of quaternion components
  T q0sq = 0.25 * (ex[0] + ey[1] + ez[2] + 1.0);
  T q1sq = q0sq - 0.5 * (ey[1] + ez[2]);
  T q2sq = q0sq - 0.5 * (ex[0] + ez[2]);
  T q3sq = q0sq - 0.5 * (ex[0] + ey[1]);

  // some component must be greater than 1/4 since they sum to 1
  // compute other components from it
  if (q0sq >= 0.25) {
    q[0] = Kokkos::sqrt(q0sq);
    q[1] = (ey[2] - ez[1]) / (4.0*q[0]);
    q[2] = (ez[0] - ex[2]) / (4.0*q[0]);
    q[3] = (ex[1] - ey[0]) / (4.0*q[0]);
  } else if (q1sq >= 0.25) {
    q[1] = Kokkos::sqrt(q1sq);
    q[0] = (ey[2] - ez[1]) / (4.0*q[1]);
    q[2] = (ey[0] + ex[1]) / (4.0*q[1]);
    q[3] = (ex[2] + ez[0]) / (4.0*q[1]);
  } else if (q2sq >= 0.25) {
    q[2] = Kokkos::sqrt(q2sq);
    q[0] = (ez[0] - ex[2]) / (4.0*q[2]);
    q[1] = (ey[0] + ex[1]) / (4.0*q[2]);
    q[3] = (ez[1] + ey[2]) / (4.0*q[2]);
  } else if (q3sq >= 0.25) {
    q[3] = Kokkos::sqrt(q3sq);
    q[0] = (ex[1] - ey[0]) / (4.0*q[3]);
    q[1] = (ez[0] + ex[2]) / (4.0*q[3]);
    q[2] = (ez[1] + ey[2]) / (4.0*q[3]);
  }
  qnormalize(q);
}

/* ----------------------------------------------------------------------
   compute rotation matrix from quaternion
   quat = [w i j k]
------------------------------------------------------------------------- */
template <typename T>
KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::quat_to_mat(const T *quat, T mat[3][3])
{
  const T w2 = quat[0]*quat[0];
  const T i2 = quat[1]*quat[1];
  const T j2 = quat[2]*quat[2];
  const T k2 = quat[3]*quat[3];
  const T twoij = 2.0*quat[1]*quat[2];
  const T twoik = 2.0*quat[1]*quat[3];
  const T twojk = 2.0*quat[2]*quat[3];
  const T twoiw = 2.0*quat[1]*quat[0];
  const T twojw = 2.0*quat[2]*quat[0];
  const T twokw = 2.0*quat[3]*quat[0];

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

/* ----------------------------------------------------------------------
   quaternion to body-frame axes (column vectors of rotation matrix)
------------------------------------------------------------------------- */
template <typename T>
KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::q_to_exyz(const T *q, T *ex, T *ey, T *ez)
{
  ex[0] = q[0]*q[0] + q[1]*q[1] - q[2]*q[2] - q[3]*q[3];
  ex[1] = 2.0 * (q[1]*q[2] + q[0]*q[3]);
  ex[2] = 2.0 * (q[1]*q[3] - q[0]*q[2]);

  ey[0] = 2.0 * (q[1]*q[2] - q[0]*q[3]);
  ey[1] = q[0]*q[0] - q[1]*q[1] + q[2]*q[2] - q[3]*q[3];
  ey[2] = 2.0 * (q[2]*q[3] + q[0]*q[1]);

  ez[0] = 2.0 * (q[1]*q[3] + q[0]*q[2]);
  ez[1] = 2.0 * (q[2]*q[3] - q[0]*q[1]);
  ez[2] = q[0]*q[0] - q[1]*q[1] - q[2]*q[2] + q[3]*q[3];
}

/* ----------------------------------------------------------------------
   quaternion-vector multiply: c = a*b, where b = (0,b)
------------------------------------------------------------------------- */
template <typename T>
KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::quatvec(const T *a, const T *b, T *c)
{
  c[0] = -a[1]*b[0] - a[2]*b[1] - a[3]*b[2];
  c[1] =  a[0]*b[0] + a[2]*b[2] - a[3]*b[1];
  c[2] =  a[0]*b[1] + a[3]*b[0] - a[1]*b[2];
  c[3] =  a[0]*b[2] + a[1]*b[1] - a[2]*b[0];
}

/* ----------------------------------------------------------------------
   inverse-quaternion-vector multiply: c = inv(a)*b
   a is a quaternion, b is a four-component vector, c is three-component
------------------------------------------------------------------------- */
template <typename T>
KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::invquatvec(const T *a, const KK_FLOAT *b, KK_FLOAT *c)
{
  c[0] = -a[1]*b[0] + a[0]*b[1] + a[3]*b[2] - a[2]*b[3];
  c[1] = -a[2]*b[0] - a[3]*b[1] + a[0]*b[2] + a[1]*b[3];
  c[2] = -a[3]*b[0] + a[2]*b[1] - a[1]*b[2] + a[0]*b[3];
}

/* ----------------------------------------------------------------------
   quaternion-quaternion multiply: c = a*b
------------------------------------------------------------------------- */
template <typename T>
KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::quatquat(const T *a, const T *b, T *c)
{
  c[0] = a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3];
  c[1] = a[0]*b[1] + b[0]*a[1] + a[2]*b[3] - a[3]*b[2];
  c[2] = a[0]*b[2] + b[0]*a[2] + a[3]*b[1] - a[1]*b[3];
  c[3] = a[0]*b[3] + b[0]*a[3] + a[1]*b[2] - a[2]*b[1];
}

/* ----------------------------------------------------------------------
   no-squish rotation for NH quaternion integration
------------------------------------------------------------------------- */
template <typename T>
KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::no_squish_rotate(int k, T *p, T *q,
                                       const KK_FLOAT *inertia, KK_FLOAT dt)
{
  KK_FLOAT phi,c_phi,s_phi;
  KK_FLOAT kp[4],kq[4];

  if (k == 1) {
    kq[0] = -q[1];  kp[0] = -p[1];
    kq[1] =  q[0];  kp[1] =  p[0];
    kq[2] =  q[3];  kp[2] =  p[3];
    kq[3] = -q[2];  kp[3] = -p[2];
  } else if (k == 2) {
    kq[0] = -q[2];  kp[0] = -p[2];
    kq[1] = -q[3];  kp[1] = -p[3];
    kq[2] =  q[0];  kp[2] =  p[0];
    kq[3] =  q[1];  kp[3] =  p[1];
  } else if (k == 3) {
    kq[0] = -q[3];  kp[0] = -p[3];
    kq[1] =  q[2];  kp[1] =  p[2];
    kq[2] = -q[1];  kp[2] = -p[1];
    kq[3] =  q[0];  kp[3] =  p[0];
  }

  phi = p[0]*kq[0] + p[1]*kq[1] + p[2]*kq[2] + p[3]*kq[3];
  if (inertia[k-1] == 0.0) phi = 0.0;
  else phi /= 4.0 * inertia[k-1];
  c_phi = Kokkos::cos(dt * phi);
  s_phi = Kokkos::sin(dt * phi);

  p[0] = fma(c_phi, p[0], s_phi*kp[0]);
  p[1] = fma(c_phi, p[1], s_phi*kp[1]);
  p[2] = fma(c_phi, p[2], s_phi*kp[2]);
  p[3] = fma(c_phi, p[3], s_phi*kp[3]);

  q[0] = fma(c_phi, q[0], s_phi*kq[0]);
  q[1] = fma(c_phi, q[1], s_phi*kq[1]);
  q[2] = fma(c_phi, q[2], s_phi*kq[2]);
  q[3] = fma(c_phi, q[3], s_phi*kq[3]);
}

/* ----------------------------------------------------------------------
   angular momentum to omega via body-frame axes
------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::angmom_to_omega(const KK_FLOAT *m, const KK_FLOAT *ex,
                                      const KK_FLOAT *ey, const KK_FLOAT *ez,
                                      const KK_FLOAT *idiag, KK_FLOAT *w)
{
  KK_FLOAT wbody[3];

  if (idiag[0] == 0.0) wbody[0] = 0.0;
  else wbody[0] = (m[0]*ex[0] + m[1]*ex[1] + m[2]*ex[2]) / idiag[0];
  if (idiag[1] == 0.0) wbody[1] = 0.0;
  else wbody[1] = (m[0]*ey[0] + m[1]*ey[1] + m[2]*ey[2]) / idiag[1];
  if (idiag[2] == 0.0) wbody[2] = 0.0;
  else wbody[2] = (m[0]*ez[0] + m[1]*ez[1] + m[2]*ez[2]) / idiag[2];

  w[0] = wbody[0]*ex[0] + wbody[1]*ey[0] + wbody[2]*ez[0];
  w[1] = wbody[0]*ex[1] + wbody[1]*ey[1] + wbody[2]*ez[1];
  w[2] = wbody[0]*ex[2] + wbody[1]*ey[2] + wbody[2]*ez[2];
}

/* ----------------------------------------------------------------------
   omega to angular momentum via body-frame axes
------------------------------------------------------------------------- */
KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::omega_to_angmom(const KK_FLOAT *w, const KK_FLOAT *ex,
                                      const KK_FLOAT *ey, const KK_FLOAT *ez,
                                      const KK_FLOAT *idiag, KK_FLOAT *m)
{
  KK_FLOAT mbody[3];

  mbody[0] = (w[0]*ex[0] + w[1]*ex[1] + w[2]*ex[2]) * idiag[0];
  mbody[1] = (w[0]*ey[0] + w[1]*ey[1] + w[2]*ey[2]) * idiag[1];
  mbody[2] = (w[0]*ez[0] + w[1]*ez[1] + w[2]*ez[2]) * idiag[2];

  m[0] = mbody[0]*ex[0] + mbody[1]*ey[0] + mbody[2]*ez[0];
  m[1] = mbody[0]*ex[1] + mbody[1]*ey[1] + mbody[2]*ez[1];
  m[2] = mbody[0]*ex[2] + mbody[1]*ey[2] + mbody[2]*ez[2];
}




//
// Device-compatible 3×3 real symmetric eigensolver, templated on scalar type T.
//
// Mixed-precision strategy:
//   T            — storage precision (float or double); input tensor + output
//   KK_ACC_FLOAT — accumulation precision; always double in LAMMPS Kokkos builds
//                  regardless of T, matching the convention used in pair styles.
//                  All intermediate Cardano arithmetic, norms, and trig run at
//                  this precision; results are cast back to T on write-out.
//
// Algorithm:
//   Eigenvalues  — Cardano trigonometric (Kopp 2008, arXiv:physics/0610206)
//   Eigenvectors — cross-product of (A-λI) rows; degenerate cases handled via
//                  analytic orthonormal complement (Hughes & Möller 1999)
//
// API is drop-in compatible with MathEigen::jacobi3:
//   sym3x3_eigen(tensor, evals, evecs, 1)   ascending sort
//   Always returns 0 (kept for parity).
//
// Convention (same as jacobi3):
//   evecs[i][0..2] = components of i-th eigenvector
//   evals[i]       = eigenvalue paired with evecs[i]

// ── Accumulation-precision helpers ──────────────────────────────────────────
// All geometry is done in KK_ACC_FLOAT regardless of T so that cross products
// and norms don't lose significance when T == float.

namespace impl {

// Two orthonormal vectors spanning the plane perpendicular to unit vector u.
// Uses the "most orthogonal axis" trick to avoid cancellation (Hughes & Möller).
KOKKOS_INLINE_FUNCTION
void orthonormal_complement(const KK_ACC_FLOAT u[3],
                                  KK_ACC_FLOAT v[3],
                                  KK_ACC_FLOAT w[3]) noexcept
{
  if (Kokkos::fabs(u[0]) > Kokkos::fabs(u[2])) {
    const KK_ACC_FLOAT inv =
        static_cast<KK_ACC_FLOAT>(1.0) /
        Kokkos::sqrt(u[0]*u[0] + u[1]*u[1]);
    v[0] = -u[1]*inv;  v[1] = u[0]*inv;  v[2] = static_cast<KK_ACC_FLOAT>(0);
  } else {
    const KK_ACC_FLOAT inv =
        static_cast<KK_ACC_FLOAT>(1.0) /
        Kokkos::sqrt(u[1]*u[1] + u[2]*u[2]);
    v[0] = static_cast<KK_ACC_FLOAT>(0);
    v[1] = -u[2]*inv;
    v[2] =  u[1]*inv;
  }
  MathExtraKokkos::cross3(u, v, w);
}


// Compute the eigenvector for eigenvalue lam as the cross product of the pair
// of rows of (A - lam·I) with the largest combined norm.
// Returns squared norm of the winning cross product.
// All six matrix entries are passed as KK_ACC_FLOAT (already upcast at call site).
KOKKOS_INLINE_FUNCTION
KK_ACC_FLOAT evec_from_eval(const KK_ACC_FLOAT a00,
                             const KK_ACC_FLOAT a11,
                             const KK_ACC_FLOAT a22,
                             const KK_ACC_FLOAT a01,
                             const KK_ACC_FLOAT a02,
                             const KK_ACC_FLOAT a12,
                             const KK_ACC_FLOAT lam,
                                   KK_ACC_FLOAT ev[3]) noexcept
{
  const KK_ACC_FLOAT r0[3] = { a00-lam, a01,     a02     };
  const KK_ACC_FLOAT r1[3] = { a01,     a11-lam, a12     };
  const KK_ACC_FLOAT r2[3] = { a02,     a12,     a22-lam };

  KK_ACC_FLOAT c01[3], c02[3], c12[3];
  MathExtraKokkos::cross3(r0, r1, c01);
  MathExtraKokkos::cross3(r0, r2, c02);
  MathExtraKokkos::cross3(r1, r2, c12);

  const KK_ACC_FLOAT n01 = MathExtraKokkos::dot3(c01, c01);
  const KK_ACC_FLOAT n02 = MathExtraKokkos::dot3(c02, c02);
  const KK_ACC_FLOAT n12 = MathExtraKokkos::dot3(c12, c12);

  KK_ACC_FLOAT best_n;
  if (n01 >= n02 && n01 >= n12) {
    ev[0]=c01[0]; ev[1]=c01[1]; ev[2]=c01[2]; best_n=n01;
  } else if (n02 >= n12) {
    ev[0]=c02[0]; ev[1]=c02[1]; ev[2]=c02[2]; best_n=n02;
  } else {
    ev[0]=c12[0]; ev[1]=c12[1]; ev[2]=c12[2]; best_n=n12;
  }
  return best_n;
}

} // !namespace impl

// ─────────────────────────────────────────────────────────────────────────────
//  sym3x3_eigen<T>
//
//  T       storage type  — float or double; controls input/output arrays
//  A       [in]   3×3 symmetric matrix (upper = lower, as the caller sets up)
//  evals   [out]  3 eigenvalues in type T
//  evecs   [out]  evecs[i][j] = j-th component of i-th eigenvector, in type T
//  sort    [in]   1 → ascending by eigenvalue (matches jacobi3 behaviour)
//
//  Returns 0 always (parity with MathEigen::jacobi3 error code).
// ─────────────────────────────────────────────────────────────────────────────
template <typename T>
KOKKOS_INLINE_FUNCTION
int MathExtraKokkos::sym3x3_eigen(const T A[3][3], T evals[3], T evecs[3][3], int sort) noexcept
{
  using acc_t = KK_ACC_FLOAT;   // accumulation type alias for readability

  // ── Upcast matrix entries to accumulation precision once ──────────────────
  const acc_t a00 = static_cast<acc_t>(A[0][0]);
  const acc_t a11 = static_cast<acc_t>(A[1][1]);
  const acc_t a22 = static_cast<acc_t>(A[2][2]);
  const acc_t a01 = static_cast<acc_t>(A[0][1]);
  const acc_t a02 = static_cast<acc_t>(A[0][2]);
  const acc_t a12 = static_cast<acc_t>(A[1][2]);

  // Working eigenvalues and eigenvector rows in accumulation precision.
  // Cast to T only on final write-out.
  acc_t ev[3];
  acc_t ew[3][3];

  // ── 1. Eigenvalues (Cardano) ──────────────────────────────────────────────

  // Off-diagonal Frobenius norm squared; if zero, A is already diagonal.
  const acc_t p1 = a01*a01 + a02*a02 + a12*a12;

  constexpr acc_t ZERO  = static_cast<acc_t>(0);
  constexpr acc_t ONE   = static_cast<acc_t>(1);
  constexpr acc_t TWO   = static_cast<acc_t>(2);
  constexpr acc_t THREE = static_cast<acc_t>(3);
  constexpr acc_t SIXTH = ONE / static_cast<acc_t>(6);

  // 2π/3 in accumulation precision
  constexpr acc_t TWO_PI_3 = static_cast<acc_t>(2.09439510239319550457);

  if (p1 < static_cast<acc_t>(1.0e-280)) {
    // Diagonal: trivial
    ev[0] = a00;  ev[1] = a11;  ev[2] = a22;
    ew[0][0]=ONE;  ew[0][1]=ZERO; ew[0][2]=ZERO;
    ew[1][0]=ZERO; ew[1][1]=ONE;  ew[1][2]=ZERO;
    ew[2][0]=ZERO; ew[2][1]=ZERO; ew[2][2]=ONE;

  } else {
    // Shift to trace-zero: B = A - q·I
    const acc_t q  = (a00 + a11 + a22) * (ONE/THREE);
    const acc_t b0 = a00 - q,  b1 = a11 - q,  b2 = a22 - q;

    // p = ||B||_F / sqrt(6)
    const acc_t p2 = (b0*b0 + b1*b1 + b2*b2 + TWO*p1) * SIXTH;
    const acc_t p  = Kokkos::sqrt(p2);
    const acc_t ip = ONE / p;

    // Scaled deviatoric matrix C = B/p
    const acc_t C00=b0*ip, C11=b1*ip, C22=b2*ip;
    const acc_t C01=a01*ip, C02=a02*ip, C12=a12*ip;

    // r = det(C)/2, clamped to [-1,1]
    acc_t r = static_cast<acc_t>(0.5) *
              (  C00*(C11*C22 - C12*C12)
               - C01*(C01*C22 - C12*C02)
               + C02*(C01*C12 - C11*C02) );
    r = (r < -ONE) ? -ONE : ((r > ONE) ? ONE : r);

    // Three roots via Cardano trig formula
    //   k=0 → middle, k=1 → largest, k=2 → smallest
    const acc_t phi = Kokkos::acos(r) * (ONE/THREE);
    ev[2] = q + TWO*p*Kokkos::cos(phi);
    ev[0] = q + TWO*p*Kokkos::cos(phi + TWO_PI_3);
    ev[1] = THREE*q - ev[0] - ev[2];  // trace residual: more stable than 3rd cos

    // ── 2. Eigenvectors ───────────────────────────────────────────────────
    const acc_t tol = static_cast<acc_t>(1.0e-8) *
                      (Kokkos::fabs(ev[2]) + Kokkos::fabs(ev[0]));
    const bool deg01 = (ev[1] - ev[0]) < tol;
    const bool deg12 = (ev[2] - ev[1]) < tol;

    if (deg01 && deg12) {
      // Triple degeneracy: any orthonormal frame
      ew[0][0]=ONE;  ew[0][1]=ZERO; ew[0][2]=ZERO;
      ew[1][0]=ZERO; ew[1][1]=ONE;  ew[1][2]=ZERO;
      ew[2][0]=ZERO; ew[2][1]=ZERO; ew[2][2]=ONE;

    } else if (deg01) {
      // ev[0] and ev[1] coalesce: compute distinct ev[2] then complement
      impl::evec_from_eval(a00,a11,a22,a01,a02,a12, ev[2], ew[2]);
      MathExtraKokkos::norm3(ew[2]);
      impl::orthonormal_complement(ew[2], ew[0], ew[1]);

    } else if (deg12) {
      // ev[1] and ev[2] coalesce: compute distinct ev[0] then complement
      impl::evec_from_eval(a00,a11,a22,a01,a02,a12, ev[0], ew[0]);
      MathExtraKokkos::norm3(ew[0]);
      impl::orthonormal_complement(ew[0], ew[1], ew[2]);

    } else {
      // All distinct: cross-product method for each
      for (int k = 0; k < 3; ++k) {
        impl::evec_from_eval(a00,a11,a22,a01,a02,a12, ev[k], ew[k]);
        MathExtraKokkos::norm3(ew[k]);
      }
    }
  }

  // ── 3. DESCENDING sort by eigenvalue (largest to smallest) ───────────────
  // Standard LAMMPS jacobi() sorts descending. Cardano outputs ascending, 
  // so we reverse the insertion sort check.
  if (sort) {
    for (int i = 1; i < 3; ++i) {
      const acc_t key_e  = ev[i];
      const acc_t key_v0 = ew[i][0];
      const acc_t key_v1 = ew[i][1];
      const acc_t key_v2 = ew[i][2];
      int j = i - 1;
      // FIX: Changed > to < to sort in DESCENDING order
      while (j >= 0 && ev[j] < key_e) { 
        ev[j+1]    = ev[j];
        ew[j+1][0] = ew[j][0];
        ew[j+1][1] = ew[j][1];
        ew[j+1][2] = ew[j][2];
        --j;
      }
      ev[j+1]    = key_e;
      ew[j+1][0] = key_v0;
      ew[j+1][1] = key_v1;
      ew[j+1][2] = key_v2;
    }
  }

  // ── 3.5 Enforce Right-Handed Coordinate System ───────────────────────────
  // jacobi() guarantees a right-handed system because it uses continuous rotations.
  // Analytic roots do not. We must ensure ew[0] x ew[1] \cdot ew[2] > 0
  const acc_t cross0 = ew[0][1]*ew[1][2] - ew[0][2]*ew[1][1];
  const acc_t cross1 = ew[0][2]*ew[1][0] - ew[0][0]*ew[1][2];
  const acc_t cross2 = ew[0][0]*ew[1][1] - ew[0][1]*ew[1][0];
  const acc_t dot = cross0*ew[2][0] + cross1*ew[2][1] + cross2*ew[2][2];
  
  if (dot < ZERO) {
    ew[2][0] = -ew[2][0];
    ew[2][1] = -ew[2][1];
    ew[2][2] = -ew[2][2];
  }

  // ── 4. Downcast to storage precision T and write out as COLUMNS ──────────
  // FIX: jacobi() stores the i-th eigenvector in the i-th COLUMN, not row.
  for (int i = 0; i < 3; ++i) {
    evals[i] = static_cast<T>(ev[i]);
    evecs[0][i] = static_cast<T>(ew[i][0]);
    evecs[1][i] = static_cast<T>(ew[i][1]);
    evecs[2][i] = static_cast<T>(ew[i][2]);
  }

  return 0;

}

#endif // !LMP_MATH_EXTRA_KOKKOS_H
