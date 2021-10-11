/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Mike Brown (SNL)
------------------------------------------------------------------------- */

#ifndef LMP_MATH_EXTRA_KOKKOS_H
#define LMP_MATH_EXTRA_KOKKOS_H

#include "kokkos_type.h"

#include <cmath>

namespace MathExtraKokkos {

// 3 vector operations

KOKKOS_INLINE_FUNCTION
void copy3(const double *v, double *ans);
KOKKOS_INLINE_FUNCTION
void zero3(double *v);
KOKKOS_INLINE_FUNCTION
void norm3(double *v);
KOKKOS_INLINE_FUNCTION
void normalize3(const double *v, double *ans);
KOKKOS_INLINE_FUNCTION
void snormalize3(const double, const double *v, double *ans);
KOKKOS_INLINE_FUNCTION
void negate3(double *v);
KOKKOS_INLINE_FUNCTION
void scale3(const double s, double *v);
KOKKOS_INLINE_FUNCTION
void scale3(const double s, const double *v, double *ans);
KOKKOS_INLINE_FUNCTION
void add3(const double *v1, const double *v2, double *ans);
KOKKOS_INLINE_FUNCTION
void scaleadd3(const double s, const double *v1, const double *v2, double *ans);
KOKKOS_INLINE_FUNCTION
void scaleadd3(const double s1, const double *v1, const double s2, const double *v2,
                      double *ans);
KOKKOS_INLINE_FUNCTION
void sub3(const double *v1, const double *v2, double *ans);
KOKKOS_INLINE_FUNCTION
double len3(const double *v);
KOKKOS_INLINE_FUNCTION
double lensq3(const double *v);
KOKKOS_INLINE_FUNCTION
double distsq3(const double *v1, const double *v2);
KOKKOS_INLINE_FUNCTION
double dot3(const double *v1, const double *v2);
KOKKOS_INLINE_FUNCTION
void cross3(const double *v1, const double *v2, double *ans);

// 3x3 matrix operations

KOKKOS_INLINE_FUNCTION
void zeromat3(double m[3][3]);
KOKKOS_INLINE_FUNCTION
void zeromat3(double **m);

KOKKOS_INLINE_FUNCTION
void col2mat(const double *ex, const double *ey, const double *ez, double m[3][3]);
KOKKOS_INLINE_FUNCTION
double det3(const double mat[3][3]);
KOKKOS_INLINE_FUNCTION
void diag_times3(const double *d, const double m[3][3], double ans[3][3]);
KOKKOS_INLINE_FUNCTION
void times3_diag(const double m[3][3], const double *d, double ans[3][3]);
KOKKOS_INLINE_FUNCTION
void plus3(const double m[3][3], const double m2[3][3], double ans[3][3]);
KOKKOS_INLINE_FUNCTION
void plus3(const double m[3][3], double **m2, double **ans);
KOKKOS_INLINE_FUNCTION
void minus3(const double m[3][3], const double m2[3][3], double ans[3][3]);
KOKKOS_INLINE_FUNCTION
void minus3(double **m, const double m2[3][3], double ans[3][3]);

KOKKOS_INLINE_FUNCTION
void times3(const double m[3][3], const double m2[3][3], double ans[3][3]);
KOKKOS_INLINE_FUNCTION
void transpose_times3(const double m[3][3], const double m2[3][3], double ans[3][3]);
KOKKOS_INLINE_FUNCTION
void times3_transpose(const double m[3][3], const double m2[3][3], double ans[3][3]);
KOKKOS_INLINE_FUNCTION
void invert3(const double mat[3][3], double ans[3][3]);
KOKKOS_INLINE_FUNCTION
void matvec(const double mat[3][3], const double *vec, double *ans);
KOKKOS_INLINE_FUNCTION
void matvec(const double *ex, const double *ey, const double *ez, const double *vec,
                   double *ans);
KOKKOS_INLINE_FUNCTION
void transpose_matvec(const double mat[3][3], const double *vec, double *ans);
KOKKOS_INLINE_FUNCTION
void transpose_matvec(const double *ex, const double *ey, const double *ez, const double *v,
                             double *ans);
KOKKOS_INLINE_FUNCTION
void transpose_diag3(const double m[3][3], const double *d, double ans[3][3]);
KOKKOS_INLINE_FUNCTION
void vecmat(const double *v, const double m[3][3], double *ans);
KOKKOS_INLINE_FUNCTION
void scalar_times3(const double f, double m[3][3]);
KOKKOS_INLINE_FUNCTION
void outer3(const double *v1, const double *v2, double ans[3][3]);

KOKKOS_INLINE_FUNCTION
void write3(const double mat[3][3]);
KOKKOS_INLINE_FUNCTION
int mldivide3(const double mat[3][3], const double *vec, double *ans);
KOKKOS_INLINE_FUNCTION
void rotate(double matrix[3][3], int i, int j, int k, int l, double s, double tau);
KOKKOS_INLINE_FUNCTION
void richardson(double *q, double *m, double *w, double *moments, double dtq);
KOKKOS_INLINE_FUNCTION
void no_squish_rotate(int k, double *p, double *q, double *inertia, double dt);

// shape matrix operations
// upper-triangular 3x3 matrix stored in Voigt ordering as 6-vector

KOKKOS_INLINE_FUNCTION
void multiply_shape_shape(const double *one, const double *two, double *ans);

// quaternion operations

KOKKOS_INLINE_FUNCTION
void qnormalize(double *q);
KOKKOS_INLINE_FUNCTION
void qconjugate(double *q, double *qc);
KOKKOS_INLINE_FUNCTION
void vecquat(double *a, double *b, double *c);
KOKKOS_INLINE_FUNCTION
void quatvec(double *a, double *b, double *c);
KOKKOS_INLINE_FUNCTION
void quatquat(double *a, double *b, double *c);
KOKKOS_INLINE_FUNCTION
void invquatvec(double *a, double *b, double *c);
KOKKOS_INLINE_FUNCTION
void axisangle_to_quat(const double *v, const double angle, double *quat);

void angmom_to_omega(double *m, double *ex, double *ey, double *ez, double *idiag, double *w);
void omega_to_angmom(double *w, double *ex, double *ey, double *ez, double *idiag, double *m);
void mq_to_omega(double *m, double *q, double *moments, double *w);
void exyz_to_q(double *ex, double *ey, double *ez, double *q);
void q_to_exyz(double *q, double *ex, double *ey, double *ez);
void quat_to_mat(const double *quat, double mat[3][3]);
void quat_to_mat_trans(const double *quat, double mat[3][3]);

// rotation operations

KOKKOS_INLINE_FUNCTION
void rotation_generator_x(const double m[3][3], double ans[3][3]);
KOKKOS_INLINE_FUNCTION
void rotation_generator_y(const double m[3][3], double ans[3][3]);
KOKKOS_INLINE_FUNCTION
void rotation_generator_z(const double m[3][3], double ans[3][3]);

KOKKOS_INLINE_FUNCTION
void BuildRxMatrix(double R[3][3], const double angle);
KOKKOS_INLINE_FUNCTION
void BuildRyMatrix(double R[3][3], const double angle);
KOKKOS_INLINE_FUNCTION
void BuildRzMatrix(double R[3][3], const double angle);

// moment of inertia operations

KOKKOS_INLINE_FUNCTION
void inertia_ellipsoid(double *shape, double *quat, double mass, double *inertia);
KOKKOS_INLINE_FUNCTION
void inertia_line(double length, double theta, double mass, double *inertia);
KOKKOS_INLINE_FUNCTION
void inertia_triangle(double *v0, double *v1, double *v2, double mass, double *inertia);
KOKKOS_INLINE_FUNCTION
void inertia_triangle(double *idiag, double *quat, double mass, double *inertia);
}    // namespace MathExtraKokkos

/* ----------------------------------------------------------------------
   copy a vector, return in ans
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::copy3(const double *v, double *ans)
{
  ans[0] = v[0];
  ans[1] = v[1];
  ans[2] = v[2];
}

/* ----------------------------------------------------------------------
   set vector equal to zero
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::zero3(double *v)
{
  v[0] = 0.0;
  v[1] = 0.0;
  v[2] = 0.0;
}

/* ----------------------------------------------------------------------
   normalize a vector in place
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::norm3(double *v)
{
  const double val = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
  if (val > 0.0) {
    const double scale = 1.0 / sqrt(val);
    v[0] *= scale;
    v[1] *= scale;
    v[2] *= scale;
  }
}

/* ----------------------------------------------------------------------
   normalize a vector, return in ans
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::normalize3(const double *v, double *ans)
{
  const double val = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
  if (val > 0.0) {
    double scale = 1.0 / sqrt(val);
    ans[0] = v[0] * scale;
    ans[1] = v[1] * scale;
    ans[2] = v[2] * scale;
  }
}

/* ----------------------------------------------------------------------
   scale a vector to length
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::snormalize3(const double length, const double *v, double *ans)
{
  const double val = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
  if (val > 0.0) {
    double scale = length / sqrt(val);
    ans[0] = v[0] * scale;
    ans[1] = v[1] * scale;
    ans[2] = v[2] * scale;
  }
}

/* ----------------------------------------------------------------------
   negate vector v
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::negate3(double *v)
{
  v[0] = -v[0];
  v[1] = -v[1];
  v[2] = -v[2];
}

/* ----------------------------------------------------------------------
   scale vector v by s
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::scale3(const double s, double *v)
{
  v[0] *= s;
  v[1] *= s;
  v[2] *= s;
}

/* ----------------------------------------------------------------------
   scale vector v by s and store in ans
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::scale3(const double s, const double *v, double *ans)
{
  ans[0] = s * v[0];
  ans[1] = s * v[1];
  ans[2] = s * v[2];
}

/* ----------------------------------------------------------------------
   ans = v1 + v2
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::add3(const double *v1, const double *v2, double *ans)
{
  ans[0] = v1[0] + v2[0];
  ans[1] = v1[1] + v2[1];
  ans[2] = v1[2] + v2[2];
}

/* ----------------------------------------------------------------------
   ans = s*v1 + v2
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::scaleadd3(const double s, const double *v1, const double *v2, double *ans)
{
  ans[0] = s * v1[0] + v2[0];
  ans[1] = s * v1[1] + v2[1];
  ans[2] = s * v1[2] + v2[2];
}

/* ----------------------------------------------------------------------
   ans = s1*v1 + s2*v2
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::scaleadd3(const double s1, const double *v1, const double s2,
                                 const double *v2, double *ans)
{
  ans[0] = s1 * v1[0] + s2 * v2[0];
  ans[1] = s1 * v1[1] + s2 * v2[1];
  ans[2] = s1 * v1[2] + s2 * v2[2];
}

/* ----------------------------------------------------------------------
   ans = v1 - v2
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::sub3(const double *v1, const double *v2, double *ans)
{
  ans[0] = v1[0] - v2[0];
  ans[1] = v1[1] - v2[1];
  ans[2] = v1[2] - v2[2];
}

/* ----------------------------------------------------------------------
   length of vector v
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
double MathExtraKokkos::len3(const double *v)
{
  return sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

/* ----------------------------------------------------------------------
   squared length of vector v, or dot product of v with itself
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
double MathExtraKokkos::lensq3(const double *v)
{
  return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
}

/* ----------------------------------------------------------------------
   ans = distance squared between pts v1 and v2
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
double MathExtraKokkos::distsq3(const double *v1, const double *v2)
{
  double dx = v1[0] - v2[0];
  double dy = v1[1] - v2[1];
  double dz = v1[2] - v2[2];
  return dx * dx + dy * dy + dz * dz;
}

/* ----------------------------------------------------------------------
   dot product of 2 vectors
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
double MathExtraKokkos::dot3(const double *v1, const double *v2)
{
  return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
}

/* ----------------------------------------------------------------------
   cross product of 2 vectors
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::cross3(const double *v1, const double *v2, double *ans)
{
  ans[0] = v1[1] * v2[2] - v1[2] * v2[1];
  ans[1] = v1[2] * v2[0] - v1[0] * v2[2];
  ans[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

/* ----------------------------------------------------------------------
   construct matrix from 3 column vectors
------------------------------------------------------------------------- */

void MathExtraKokkos::col2mat(const double *ex, const double *ey, const double *ez, double m[3][3])
{
  m[0][0] = ex[0];
  m[1][0] = ex[1];
  m[2][0] = ex[2];
  m[0][1] = ey[0];
  m[1][1] = ey[1];
  m[2][1] = ey[2];
  m[0][2] = ez[0];
  m[1][2] = ez[1];
  m[2][2] = ez[2];
}

/* ----------------------------------------------------------------------
   determinant of a matrix
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
double MathExtraKokkos::det3(const double m[3][3])
{
  double ans = m[0][0] * m[1][1] * m[2][2] - m[0][0] * m[1][2] * m[2][1] -
      m[1][0] * m[0][1] * m[2][2] + m[1][0] * m[0][2] * m[2][1] + m[2][0] * m[0][1] * m[1][2] -
      m[2][0] * m[0][2] * m[1][1];
  return ans;
}

/* ----------------------------------------------------------------------
   diagonal matrix times a full matrix
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::diag_times3(const double *d, const double m[3][3], double ans[3][3])
{
  ans[0][0] = d[0] * m[0][0];
  ans[0][1] = d[0] * m[0][1];
  ans[0][2] = d[0] * m[0][2];
  ans[1][0] = d[1] * m[1][0];
  ans[1][1] = d[1] * m[1][1];
  ans[1][2] = d[1] * m[1][2];
  ans[2][0] = d[2] * m[2][0];
  ans[2][1] = d[2] * m[2][1];
  ans[2][2] = d[2] * m[2][2];
}

/* ----------------------------------------------------------------------
   full matrix times a diagonal matrix
------------------------------------------------------------------------- */

void MathExtraKokkos::times3_diag(const double m[3][3], const double *d, double ans[3][3])
{
  ans[0][0] = m[0][0] * d[0];
  ans[0][1] = m[0][1] * d[1];
  ans[0][2] = m[0][2] * d[2];
  ans[1][0] = m[1][0] * d[0];
  ans[1][1] = m[1][1] * d[1];
  ans[1][2] = m[1][2] * d[2];
  ans[2][0] = m[2][0] * d[0];
  ans[2][1] = m[2][1] * d[1];
  ans[2][2] = m[2][2] * d[2];
}

/* ----------------------------------------------------------------------
   add two matrices
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::plus3(const double m[3][3], const double m2[3][3], double ans[3][3])
{
  ans[0][0] = m[0][0] + m2[0][0];
  ans[0][1] = m[0][1] + m2[0][1];
  ans[0][2] = m[0][2] + m2[0][2];
  ans[1][0] = m[1][0] + m2[1][0];
  ans[1][1] = m[1][1] + m2[1][1];
  ans[1][2] = m[1][2] + m2[1][2];
  ans[2][0] = m[2][0] + m2[2][0];
  ans[2][1] = m[2][1] + m2[2][1];
  ans[2][2] = m[2][2] + m2[2][2];
}

/* ----------------------------------------------------------------------
   multiply mat1 times mat2
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::times3(const double m[3][3], const double m2[3][3], double ans[3][3])
{
  ans[0][0] = m[0][0] * m2[0][0] + m[0][1] * m2[1][0] + m[0][2] * m2[2][0];
  ans[0][1] = m[0][0] * m2[0][1] + m[0][1] * m2[1][1] + m[0][2] * m2[2][1];
  ans[0][2] = m[0][0] * m2[0][2] + m[0][1] * m2[1][2] + m[0][2] * m2[2][2];
  ans[1][0] = m[1][0] * m2[0][0] + m[1][1] * m2[1][0] + m[1][2] * m2[2][0];
  ans[1][1] = m[1][0] * m2[0][1] + m[1][1] * m2[1][1] + m[1][2] * m2[2][1];
  ans[1][2] = m[1][0] * m2[0][2] + m[1][1] * m2[1][2] + m[1][2] * m2[2][2];
  ans[2][0] = m[2][0] * m2[0][0] + m[2][1] * m2[1][0] + m[2][2] * m2[2][0];
  ans[2][1] = m[2][0] * m2[0][1] + m[2][1] * m2[1][1] + m[2][2] * m2[2][1];
  ans[2][2] = m[2][0] * m2[0][2] + m[2][1] * m2[1][2] + m[2][2] * m2[2][2];
}

/* ----------------------------------------------------------------------
   multiply the transpose of mat1 times mat2
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::transpose_times3(const double m[3][3], const double m2[3][3],
                                        double ans[3][3])
{
  ans[0][0] = m[0][0] * m2[0][0] + m[1][0] * m2[1][0] + m[2][0] * m2[2][0];
  ans[0][1] = m[0][0] * m2[0][1] + m[1][0] * m2[1][1] + m[2][0] * m2[2][1];
  ans[0][2] = m[0][0] * m2[0][2] + m[1][0] * m2[1][2] + m[2][0] * m2[2][2];
  ans[1][0] = m[0][1] * m2[0][0] + m[1][1] * m2[1][0] + m[2][1] * m2[2][0];
  ans[1][1] = m[0][1] * m2[0][1] + m[1][1] * m2[1][1] + m[2][1] * m2[2][1];
  ans[1][2] = m[0][1] * m2[0][2] + m[1][1] * m2[1][2] + m[2][1] * m2[2][2];
  ans[2][0] = m[0][2] * m2[0][0] + m[1][2] * m2[1][0] + m[2][2] * m2[2][0];
  ans[2][1] = m[0][2] * m2[0][1] + m[1][2] * m2[1][1] + m[2][2] * m2[2][1];
  ans[2][2] = m[0][2] * m2[0][2] + m[1][2] * m2[1][2] + m[2][2] * m2[2][2];
}

/* ----------------------------------------------------------------------
   multiply mat1 times transpose of mat2
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::times3_transpose(const double m[3][3], const double m2[3][3],
                                        double ans[3][3])
{
  ans[0][0] = m[0][0] * m2[0][0] + m[0][1] * m2[0][1] + m[0][2] * m2[0][2];
  ans[0][1] = m[0][0] * m2[1][0] + m[0][1] * m2[1][1] + m[0][2] * m2[1][2];
  ans[0][2] = m[0][0] * m2[2][0] + m[0][1] * m2[2][1] + m[0][2] * m2[2][2];
  ans[1][0] = m[1][0] * m2[0][0] + m[1][1] * m2[0][1] + m[1][2] * m2[0][2];
  ans[1][1] = m[1][0] * m2[1][0] + m[1][1] * m2[1][1] + m[1][2] * m2[1][2];
  ans[1][2] = m[1][0] * m2[2][0] + m[1][1] * m2[2][1] + m[1][2] * m2[2][2];
  ans[2][0] = m[2][0] * m2[0][0] + m[2][1] * m2[0][1] + m[2][2] * m2[0][2];
  ans[2][1] = m[2][0] * m2[1][0] + m[2][1] * m2[1][1] + m[2][2] * m2[1][2];
  ans[2][2] = m[2][0] * m2[2][0] + m[2][1] * m2[2][1] + m[2][2] * m2[2][2];
}

/* ----------------------------------------------------------------------
   invert a matrix
   does NOT check for singular or badly scaled matrix
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::invert3(const double m[3][3], double ans[3][3])
{
  double den = m[0][0] * m[1][1] * m[2][2] - m[0][0] * m[1][2] * m[2][1];
  den += -m[1][0] * m[0][1] * m[2][2] + m[1][0] * m[0][2] * m[2][1];
  den += m[2][0] * m[0][1] * m[1][2] - m[2][0] * m[0][2] * m[1][1];

  ans[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) / den;
  ans[0][1] = -(m[0][1] * m[2][2] - m[0][2] * m[2][1]) / den;
  ans[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) / den;
  ans[1][0] = -(m[1][0] * m[2][2] - m[1][2] * m[2][0]) / den;
  ans[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) / den;
  ans[1][2] = -(m[0][0] * m[1][2] - m[0][2] * m[1][0]) / den;
  ans[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) / den;
  ans[2][1] = -(m[0][0] * m[2][1] - m[0][1] * m[2][0]) / den;
  ans[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / den;
}

/* ----------------------------------------------------------------------
   matrix times vector
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::matvec(const double m[3][3], const double *v, double *ans)
{
  ans[0] = m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2];
  ans[1] = m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2];
  ans[2] = m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2];
}

/* ----------------------------------------------------------------------
   matrix times vector
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::matvec(const double *ex, const double *ey, const double *ez, const double *v,
                              double *ans)
{
  ans[0] = ex[0] * v[0] + ey[0] * v[1] + ez[0] * v[2];
  ans[1] = ex[1] * v[0] + ey[1] * v[1] + ez[1] * v[2];
  ans[2] = ex[2] * v[0] + ey[2] * v[1] + ez[2] * v[2];
}

/* ----------------------------------------------------------------------
   transposed matrix times vector
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::transpose_matvec(const double m[3][3], const double *v, double *ans)
{
  ans[0] = m[0][0] * v[0] + m[1][0] * v[1] + m[2][0] * v[2];
  ans[1] = m[0][1] * v[0] + m[1][1] * v[1] + m[2][1] * v[2];
  ans[2] = m[0][2] * v[0] + m[1][2] * v[1] + m[2][2] * v[2];
}

/* ----------------------------------------------------------------------
   transposed matrix times vector
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::transpose_matvec(const double *ex, const double *ey, const double *ez,
                                        const double *v, double *ans)
{
  ans[0] = ex[0] * v[0] + ex[1] * v[1] + ex[2] * v[2];
  ans[1] = ey[0] * v[0] + ey[1] * v[1] + ey[2] * v[2];
  ans[2] = ez[0] * v[0] + ez[1] * v[1] + ez[2] * v[2];
}

/* ----------------------------------------------------------------------
   transposed matrix times diagonal matrix
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::transpose_diag3(const double m[3][3], const double *d, double ans[3][3])
{
  ans[0][0] = m[0][0] * d[0];
  ans[0][1] = m[1][0] * d[1];
  ans[0][2] = m[2][0] * d[2];
  ans[1][0] = m[0][1] * d[0];
  ans[1][1] = m[1][1] * d[1];
  ans[1][2] = m[2][1] * d[2];
  ans[2][0] = m[0][2] * d[0];
  ans[2][1] = m[1][2] * d[1];
  ans[2][2] = m[2][2] * d[2];
}

/* ----------------------------------------------------------------------
   row vector times matrix
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::vecmat(const double *v, const double m[3][3], double *ans)
{
  ans[0] = v[0] * m[0][0] + v[1] * m[1][0] + v[2] * m[2][0];
  ans[1] = v[0] * m[0][1] + v[1] * m[1][1] + v[2] * m[2][1];
  ans[2] = v[0] * m[0][2] + v[1] * m[1][2] + v[2] * m[2][2];
}

/* ----------------------------------------------------------------------
   matrix times scalar, in place
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::scalar_times3(const double f, double m[3][3])
{
  m[0][0] *= f;
  m[0][1] *= f;
  m[0][2] *= f;
  m[1][0] *= f;
  m[1][1] *= f;
  m[1][2] *= f;
  m[2][0] *= f;
  m[2][1] *= f;
  m[2][2] *= f;
}

/* ----------------------------------------------------------------------
   multiply 2 shape matrices
   upper-triangular 3x3, stored as 6-vector in Voigt ordering
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::multiply_shape_shape(const double *one, const double *two, double *ans)
{
  ans[0] = one[0] * two[0];
  ans[1] = one[1] * two[1];
  ans[2] = one[2] * two[2];
  ans[3] = one[1] * two[3] + one[3] * two[2];
  ans[4] = one[0] * two[4] + one[5] * two[3] + one[4] * two[2];
  ans[5] = one[0] * two[5] + one[5] * two[1];
}

/* ----------------------------------------------------------------------
   normalize a quaternion
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::qnormalize(double *q)
{
  double norm = 1.0 / sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
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
void MathExtraKokkos::qconjugate(double *q, double *qc)
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
void MathExtraKokkos::vecquat(double *a, double *b, double *c)
{
  c[0] = -a[0] * b[1] - a[1] * b[2] - a[2] * b[3];
  c[1] = b[0] * a[0] + a[1] * b[3] - a[2] * b[2];
  c[2] = b[0] * a[1] + a[2] * b[1] - a[0] * b[3];
  c[3] = b[0] * a[2] + a[0] * b[2] - a[1] * b[1];
}

/* ----------------------------------------------------------------------
   quaternion-vector multiply: c = a*b, where b = (0,b)
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::quatvec(double *a, double *b, double *c)
{
  c[0] = -a[1] * b[0] - a[2] * b[1] - a[3] * b[2];
  c[1] = a[0] * b[0] + a[2] * b[2] - a[3] * b[1];
  c[2] = a[0] * b[1] + a[3] * b[0] - a[1] * b[2];
  c[3] = a[0] * b[2] + a[1] * b[1] - a[2] * b[0];
}

/* ----------------------------------------------------------------------
   quaternion-quaternion multiply: c = a*b
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::quatquat(double *a, double *b, double *c)
{
  c[0] = a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3];
  c[1] = a[0] * b[1] + b[0] * a[1] + a[2] * b[3] - a[3] * b[2];
  c[2] = a[0] * b[2] + b[0] * a[2] + a[3] * b[1] - a[1] * b[3];
  c[3] = a[0] * b[3] + b[0] * a[3] + a[1] * b[2] - a[2] * b[1];
}

/* ----------------------------------------------------------------------
   quaternion multiply: c = inv(a)*b
   a is a quaternion
   b is a four component vector
   c is a three component vector
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::invquatvec(double *a, double *b, double *c)
{
  c[0] = -a[1] * b[0] + a[0] * b[1] + a[3] * b[2] - a[2] * b[3];
  c[1] = -a[2] * b[0] - a[3] * b[1] + a[0] * b[2] + a[1] * b[3];
  c[2] = -a[3] * b[0] + a[2] * b[1] - a[1] * b[2] + a[0] * b[3];
}

/* ----------------------------------------------------------------------
   compute quaternion from axis-angle rotation
   v MUST be a unit vector
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::axisangle_to_quat(const double *v, const double angle, double *quat)
{
  double halfa = 0.5 * angle;
  double sina = sin(halfa);
  quat[0] = cos(halfa);
  quat[1] = v[0] * sina;
  quat[2] = v[1] * sina;
  quat[3] = v[2] * sina;
}

/* ----------------------------------------------------------------------
   Apply principal rotation generator about x to rotation matrix m
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::rotation_generator_x(const double m[3][3], double ans[3][3])
{
  ans[0][0] = 0;
  ans[0][1] = -m[0][2];
  ans[0][2] = m[0][1];
  ans[1][0] = 0;
  ans[1][1] = -m[1][2];
  ans[1][2] = m[1][1];
  ans[2][0] = 0;
  ans[2][1] = -m[2][2];
  ans[2][2] = m[2][1];
}

/* ----------------------------------------------------------------------
   Apply principal rotation generator about y to rotation matrix m
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::rotation_generator_y(const double m[3][3], double ans[3][3])
{
  ans[0][0] = m[0][2];
  ans[0][1] = 0;
  ans[0][2] = -m[0][0];
  ans[1][0] = m[1][2];
  ans[1][1] = 0;
  ans[1][2] = -m[1][0];
  ans[2][0] = m[2][2];
  ans[2][1] = 0;
  ans[2][2] = -m[2][0];
}

/* ----------------------------------------------------------------------
   Apply principal rotation generator about z to rotation matrix m
------------------------------------------------------------------------- */

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::rotation_generator_z(const double m[3][3], double ans[3][3])
{
  ans[0][0] = -m[0][1];
  ans[0][1] = m[0][0];
  ans[0][2] = 0;
  ans[1][0] = -m[1][1];
  ans[1][1] = m[1][0];
  ans[1][2] = 0;
  ans[2][0] = -m[2][1];
  ans[2][1] = m[2][0];
  ans[2][2] = 0;
}

// set matrix to zero

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::zeromat3(double m[3][3])
{
  m[0][0] = m[0][1] = m[0][2] = 0.0;
  m[1][0] = m[1][1] = m[1][2] = 0.0;
  m[2][0] = m[2][1] = m[2][2] = 0.0;
}

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::zeromat3(double **m)
{
  m[0][0] = m[0][1] = m[0][2] = 0.0;
  m[1][0] = m[1][1] = m[1][2] = 0.0;
  m[2][0] = m[2][1] = m[2][2] = 0.0;
}

// add two matrices

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::plus3(const double m[3][3], double **m2, double **ans)
{
  ans[0][0] = m[0][0] + m2[0][0];
  ans[0][1] = m[0][1] + m2[0][1];
  ans[0][2] = m[0][2] + m2[0][2];
  ans[1][0] = m[1][0] + m2[1][0];
  ans[1][1] = m[1][1] + m2[1][1];
  ans[1][2] = m[1][2] + m2[1][2];
  ans[2][0] = m[2][0] + m2[2][0];
  ans[2][1] = m[2][1] + m2[2][1];
  ans[2][2] = m[2][2] + m2[2][2];
}

// subtract two matrices

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::minus3(const double m[3][3], const double m2[3][3], double ans[3][3])
{
  ans[0][0] = m[0][0] - m2[0][0];
  ans[0][1] = m[0][1] - m2[0][1];
  ans[0][2] = m[0][2] - m2[0][2];
  ans[1][0] = m[1][0] - m2[1][0];
  ans[1][1] = m[1][1] - m2[1][1];
  ans[1][2] = m[1][2] - m2[1][2];
  ans[2][0] = m[2][0] - m2[2][0];
  ans[2][1] = m[2][1] - m2[2][1];
  ans[2][2] = m[2][2] - m2[2][2];
}

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::minus3(double **m, const double m2[3][3], double ans[3][3])
{
  ans[0][0] = m[0][0] - m2[0][0];
  ans[0][1] = m[0][1] - m2[0][1];
  ans[0][2] = m[0][2] - m2[0][2];
  ans[1][0] = m[1][0] - m2[1][0];
  ans[1][1] = m[1][1] - m2[1][1];
  ans[1][2] = m[1][2] - m2[1][2];
  ans[2][0] = m[2][0] - m2[2][0];
  ans[2][1] = m[2][1] - m2[2][1];
  ans[2][2] = m[2][2] - m2[2][2];
}

// compute outer product of two vectors

KOKKOS_INLINE_FUNCTION
void MathExtraKokkos::outer3(const double *v1, const double *v2, double ans[3][3])
{
  ans[0][0] = v1[0] * v2[0];
  ans[0][1] = v1[0] * v2[1];
  ans[0][2] = v1[0] * v2[2];
  ans[1][0] = v1[1] * v2[0];
  ans[1][1] = v1[1] * v2[1];
  ans[1][2] = v1[1] * v2[2];
  ans[2][0] = v1[2] * v2[0];
  ans[2][1] = v1[2] * v2[1];
  ans[2][2] = v1[2] * v2[2];
}

#endif
