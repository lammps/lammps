/* -*- c++ -*- ----------------------------------------------------------
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
   Contributing author: Pieter J. in 't Veld (SNL)
------------------------------------------------------------------------- */

#ifndef LMP_MATH_VECTOR_H
#define LMP_MATH_VECTOR_H

#include <cmath>
#include <cstring>

#define VECTOR_NULL        {0, 0, 0}
#define SHAPE_NULL        {0, 0, 0, 0, 0, 0}
#define FORM_NULL        {0, 0, 0, 0, 0, 0}
#define MATRIX_NULL        {VECTOR_NULL, VECTOR_NULL, VECTOR_NULL}
#define VECTOR4_NULL        {0, 0, 0, 0}
#define QUATERNION_NULL        {0, 0, 0, 0}
#define FORM4_NULL        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}

#define FZERO 1e-15
#define fzero(x) (((x)>-FZERO) && ((x)<FZERO))

namespace LAMMPS_NS {

typedef double vector[3];                // 0:x  1:y  2:z
typedef int lvector[3];
typedef double shape[6];                // 0:xx 1:yy 2:zz 3:zy 4:zx 5:yx
typedef int lshape[6];                        // xy=0  xz=0  yz=0;
typedef double form[6];                        // 0:xx 1:yy 2:zz 3:zy 4:zx 5:yx
typedef int lform[6];                        // xy=yx xz=zx yz=zy;
typedef vector matrix[3];                // 00:xx 11:yy 22:zz 21:zy 20:zx 10:yx
typedef lvector lmatrix[3];
typedef double vector4[4];                // 4D vector
typedef double form4[10];                // 0:00 1:11 2:22 3:33 4:32
                                        // 5:31 6:30 7:21 8:20 9:10
                                        // 01=10 02=20 03=30 etc
typedef double quaternion[4];                // quaternion

// vector operators

inline double vec_dot(vector &a, vector &b) {                        // a.b
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }

inline void vec_null(vector &dest) {
  memset(dest, 0, sizeof(vector)); }

inline void vec_neg(vector &dest) {                                // -a
  dest[0] = -dest[0];
  dest[1] = -dest[1];
  dest[2] = -dest[2]; }

inline void vec_norm(vector &dest) {                                 // a/|a|
  double f = sqrt(vec_dot(dest, dest));
  dest[0] /= f;
  dest[1] /= f;
  dest[2] /= f; }

inline void vec_add(vector &dest, vector &src) {                // a+b
  dest[0] += src[0];
  dest[1] += src[1];
  dest[2] += src[2]; }

inline void vec_subtr(vector &dest, vector &src) {                // a-b
  dest[0] -= src[0];
  dest[1] -= src[1];
  dest[2] -= src[2]; }

inline void vec_mult(vector &dest, vector &src) {                // a*b
  dest[0] *= src[0];
  dest[1] *= src[1];
  dest[2] *= src[2]; }

inline void vec_div(vector &dest, vector &src) {                // a/b
  dest[0] /= src[0];
  dest[1] /= src[1];
  dest[2] /= src[2]; }

inline void vec_cross(vector &dest, vector &src) {                // a x b
  vector t = {
    dest[1]*src[2]-dest[2]*src[1],
    dest[2]*src[0]-dest[0]*src[2],
    dest[0]*src[1]-dest[1]*src[0]};
  memcpy(dest, t, sizeof(vector)); }

inline void vec_scalar_mult(vector &dest, double f) {                // f*a
  dest[0] *= f;
  dest[1] *= f;
  dest[2] *= f; }

inline void vec_to_lvec(lvector &dest, vector &src) {
  dest[0] = (int) src[0];
  dest[1] = (int) src[1];
  dest[2] = (int) src[2]; }

inline void lvec_to_vec(vector &dest, lvector &src) {
  dest[0] = (double) src[0];
  dest[1] = (double) src[1];
  dest[2] = (double) src[2]; }

// shape operators

inline double shape_det(shape &s) {
  return s[0]*s[1]*s[2]; }

inline void shape_null(shape &dest) {
  memset(dest, 0, sizeof(shape)); }

inline void shape_unit(shape &dest) {
  memset(dest, 0, sizeof(shape));
  dest[0] = dest[1] = dest[2] = 1.0; }

inline void shape_add(shape &dest, shape &src) {                // h_a+h_b
  dest[0] += src[0]; dest[1] += src[1]; dest[2] += src[2];
  dest[3] += src[3]; dest[4] += src[4]; dest[5] += src[5]; }

inline void shape_subtr(shape &dest, shape &src) {                // h_a-h_b
  dest[0] -= src[0]; dest[1] -= src[1]; dest[2] -= src[2];
  dest[3] -= src[3]; dest[4] -= src[4]; dest[5] -= src[5]; }

inline void shape_inv(shape &h_inv, shape &h) {                        // h^-1
  h_inv[0] = 1.0/h[0]; h_inv[1] = 1.0/h[1]; h_inv[2] = 1.0/h[2];
  h_inv[3] = -h[3]/(h[1]*h[2]);
  h_inv[4] = (h[3]*h[5]-h[1]*h[4])/(h[0]*h[1]*h[2]);
  h_inv[5] = -h[5]/(h[0]*h[1]); }

inline void shape_dot(shape &dest, shape &src) {                // h_a.h_b
  dest[3] = dest[1]*src[3]+dest[3]*src[2];
  dest[4] = dest[0]*src[4]+dest[5]*src[3]+dest[4]*src[2];
  dest[5] = dest[0]*src[5]+dest[5]*src[1];
  dest[0] *= src[0]; dest[1] *= src[1]; dest[2] *= src[2]; }

inline void shape_scalar_mult(shape &dest, double f) {                // f*h
  dest[0] *= f; dest[1] *= f; dest[2] *= f;
  dest[3] *= f; dest[4] *= f; dest[5] *= f; }

inline void shape_vec_dot(vector &dest, vector &src, shape &h) {// h.a
  dest[0] = h[0]*src[0]+h[5]*src[1]+h[4]*src[2];
  dest[1] =             h[1]*src[1]+h[3]*src[2];
  dest[2] =                         h[2]*src[2]; }

inline void vec_shape_dot(vector &dest, shape &h, vector &src) {// a.h
  dest[2] = h[4]*src[0]+h[3]*src[1]+h[2]*src[2];
  dest[1] = h[5]*src[0]+h[1]*src[1];
  dest[0] = h[0]*src[0]; }

inline void shape_to_matrix(matrix &dest, shape &h) {                // m = h
  dest[0][0] = h[0]; dest[1][0] = h[5]; dest[2][0] = h[4];
  dest[0][1] =  0.0; dest[1][1] = h[1]; dest[2][1] = h[3];
  dest[0][2] =  0.0; dest[1][2] =  0.0; dest[2][2] = h[2]; }

inline void shape_to_lshape(lshape &dest, shape &src) {
  dest[0] = (long)src[0]; dest[1] = (long)src[1]; dest[2] = (long)src[2];
  dest[3] = (long)src[3]; dest[4] = (long)src[4]; dest[5] = (long)src[5]; }

inline void lshape_to_shape(shape &dest, lshape &src) {
  dest[0] = (double)src[0]; dest[1] = (double)src[1]; dest[2] = (double)src[2];
  dest[3] = (double)src[3]; dest[4] = (double)src[4]; dest[5] = (double)src[5];}

// form operators

inline double form_det(form &m) {                                // |m|
  return m[0]*(m[1]*m[2]-m[3]*m[3])+
    m[4]*(2.0*m[3]*m[5]-m[1]*m[4])-m[2]*m[5]*m[5]; }

inline void form_null(form &dest) {
  memset(dest, 0, sizeof(form)); }

inline void form_unit(form &dest) {
  memset(dest, 0, sizeof(form));
  dest[0] = dest[1] = dest[2] = 1.0; }

inline void form_add(form &dest, form &src) {                        // m_a+m_b
  dest[0] += src[0]; dest[1] += src[1]; dest[2] += src[2];
  dest[3] += src[3]; dest[4] += src[4]; dest[5] += src[5]; }

inline void form_subtr(shape &dest, form &src) {                // m_a-m_b
  dest[0] -= src[0]; dest[1] -= src[1]; dest[2] -= src[2];
  dest[3] -= src[3]; dest[4] -= src[4]; dest[5] -= src[5]; }

inline int form_inv(form &m_inv, form &m) {                        // m^-1
  double det = form_det(m);
  if (fzero(det)) return 0;
  m_inv[0] = (m[1]*m[2]-m[3]*m[3])/det;
  m_inv[1] = (m[0]*m[2]-m[4]*m[4])/det;
  m_inv[2] = (m[0]*m[1]-m[5]*m[5])/det;
  m_inv[3] = (m[4]*m[5]-m[0]*m[3])/det;
  m_inv[4] = (m[3]*m[5]-m[1]*m[4])/det;
  m_inv[5] = (m[3]*m[4]-m[2]*m[5])/det;
  return 1; }

inline void form_dot(form &dest, form &src) {                        // m_a.m_b
  form m;
  memcpy(m, dest, sizeof(form));
  dest[0] = m[0]*src[0]+m[4]*src[4]+m[5]*src[5];
  dest[1] = m[1]*src[1]+m[3]*src[3]+m[5]*src[5];
  dest[2] = m[2]*src[2]+m[3]*src[3]+m[4]*src[4];
  dest[3] = m[3]*src[2]+m[1]*src[3]+m[5]*src[4];
  dest[4] = m[4]*src[2]+m[5]*src[3]+m[0]*src[4];
  dest[5] = m[5]*src[1]+m[4]*src[3]+m[0]*src[5]; }

inline void form_vec_dot(vector &dest, form &m) {                // m.a
  vector a;
  memcpy(a, dest, sizeof(vector));
  dest[0] = m[0]*a[0]+m[5]*a[1]+m[4]*a[2];
  dest[1] = m[5]*a[0]+m[1]*a[1]+m[3]*a[2];
  dest[2] = m[4]*a[0]+m[3]*a[1]+m[2]*a[2]; }

inline void form_to_lform(lform &dest, form &src) {
  dest[0] = (long)src[0]; dest[1] = (long)src[1]; dest[2] = (long)src[2];
  dest[3] = (long)src[3]; dest[4] = (long)src[4]; dest[5] = (long)src[5]; }

inline void lform_to_form(form &dest, lform &src) {
  dest[0] = (double)src[0]; dest[1] = (double)src[1]; dest[2] = (double)src[2];
  dest[3] = (double)src[3]; dest[4] = (double)src[4]; dest[5] = (double)src[5];}

// matrix operators

inline double matrix_det(matrix &m) {                                // |m|
  vector axb;
  memcpy(&axb, m[0], sizeof(vector));
  vec_cross(axb, m[1]);
  return vec_dot(axb, m[2]); }

inline void matrix_null(matrix &dest) {
  memset(dest, 0, sizeof(dest)); }

inline void matrix_unit(matrix &dest) {
  memset(dest, 0, sizeof(dest));
  dest[0][0] = dest[1][1] = dest[2][2] = 1.0; }

inline void matrix_scalar_mult(matrix &dest, double f) {        // f*m
  vec_scalar_mult(dest[0], f);
  vec_scalar_mult(dest[1], f);
  vec_scalar_mult(dest[2], f); }

inline void matrix_trans(matrix &dest) {                        // m^t
  double f = dest[0][1]; dest[0][1] = dest[1][0]; dest[1][0] = f;
  f = dest[0][2]; dest[0][2] = dest[2][0]; dest[2][0] = f;
  f = dest[1][2]; dest[1][2] = dest[2][1]; dest[2][1] = f; }

inline int matrix_inv(matrix &dest, matrix &src) {                // m^-1
  double f = matrix_det(src);
  if (fzero(f)) return 0;                                // singular matrix
  memcpy(dest[0], src[1], sizeof(vector));
  memcpy(dest[1], src[2], sizeof(vector));
  memcpy(dest[2], src[0], sizeof(vector));
  vec_cross(dest[0], src[2]);
  vec_cross(dest[1], src[0]);
  vec_cross(dest[2], src[1]);
  matrix_scalar_mult(dest, 1.0/f);
  matrix_trans(dest);
  return 0; }

inline void matrix_vec_dot(vector &dest, vector &src, matrix &m) { // m.a
  dest[0] = m[0][0]*src[0]+m[1][0]*src[1]+m[2][0]*src[2];
  dest[1] = m[0][1]*src[0]+m[1][1]*src[1]+m[2][1]*src[2];
  dest[2] = m[0][2]*src[0]+m[1][2]*src[1]+m[2][2]*src[2]; }

inline void matrix_to_shape(shape &dest, matrix &src) {                // h = m
  dest[0] = src[0][0]; dest[1] = src[1][1]; dest[2] = src[2][2];
  dest[3] = src[2][1]; dest[4] = src[2][0]; dest[5] = src[1][0]; }

inline void matrix_to_lmatrix(lmatrix &dest, matrix &src) {
  vec_to_lvec(dest[0], src[0]);
  vec_to_lvec(dest[1], src[1]);
  vec_to_lvec(dest[2], src[2]); }

inline void lmatrix_to_matrix(matrix &dest, lmatrix &src) {
  lvec_to_vec(dest[0], src[0]);
  lvec_to_vec(dest[1], src[1]);
  lvec_to_vec(dest[2], src[2]); }

// quaternion operators

inline double quat_dot(quaternion &p, quaternion &q) {                // p.q
  return p[0]*q[0]+p[1]*q[1]+p[2]*q[2]+p[3]*q[3];
}

inline void quat_norm(quaternion &q) {                                // q = q/|q|
  double f = sqrt(q[0]*q[0]+q[1]*q[1]+q[2]*q[2]+q[3]*q[3]);
  q[0] /= f; q[1] /= f; q[2] /= f; q[3] /= f;
}

inline void quat_conj(quaternion &q) {                                // q = conj(q)
  q[1] = -q[1]; q[2] = -q[2]; q[3] = -q[3];
}

inline void quat_mult(quaternion &dest, quaternion &src) {        // dest *= src
  quaternion q;
  memcpy(q, dest, sizeof(quaternion));
  dest[0] =  src[0]*q[3]+src[1]*q[2]-src[2]*q[1]+src[3]*q[0];
  dest[1] = -src[0]*q[2]+src[1]*q[3]+src[2]*q[0]+src[3]*q[1];
  dest[2] =  src[0]*q[1]-src[1]*q[0]+src[2]*q[3]+src[3]*q[2];
  dest[3] = -src[0]*q[0]-src[1]*q[1]-src[2]*q[2]+src[3]*q[3];
}

inline void quat_div(quaternion &dest, quaternion &src) {        // dest /= src
  quaternion q;
  memcpy(q, dest, sizeof(quaternion));
  dest[0] =  src[0]*q[3]-src[1]*q[2]+src[2]*q[1]-src[3]*q[0];
  dest[1] = -src[0]*q[2]-src[1]*q[3]-src[2]*q[0]-src[3]*q[1];
  dest[2] =  src[0]*q[1]+src[1]*q[0]-src[2]*q[3]-src[3]*q[2];
  dest[3] = -src[0]*q[0]+src[1]*q[1]+src[2]*q[2]-src[3]*q[3];
}
                                                        // dest = q*src*conj(q)
inline void quat_vec_rot(vector &dest, vector &src, quaternion &q) {
  quaternion aa={q[0]*q[0], q[1]*q[1], q[2]*q[2], q[3]*q[3]};
  form ab={q[0]*q[1], q[0]*q[2], q[0]*q[3], q[1]*q[2], q[1]*q[3], q[2]*q[3]};
  dest[0] = (aa[0]+aa[1]-aa[2]-aa[3])*src[0]+
            ((ab[3]-ab[2])*src[1]+(ab[1]+ab[4])*src[2])*2.0;
  dest[1] = (aa[0]-aa[1]+aa[2]-aa[3])*src[1]+
            ((ab[2]+ab[3])*src[0]+(ab[5]-ab[0])*src[2])*2.0;
  dest[2] = (aa[0]-aa[1]-aa[2]+aa[3])*src[2]+
            ((ab[4]-ab[1])*src[0]+(ab[0]+ab[5])*src[1])*2.0;
}

// vector4 operators

inline void vec4_null(vector4 &dest) {
  memset(dest, 0, sizeof(vector4));
}

inline double vec4_dot(vector4 &a, vector4 &b) {
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]+a[3]*b[3]; }

// form operators

inline void form4_null(form4 &dest) {
  memset(dest, 0, sizeof(form4)); }

inline void form4_unit(form4 &dest) {
  memset(dest, 0, sizeof(form4));
  dest[0] = dest[1] = dest[2] = dest[3] = 1.0; }

inline double form4_det(form4 &m) {
  double f = m[6]*m[7]-m[5]*m[8];
  return m[0]*(
      m[1]*(m[2]*m[3]-m[4]*m[4])+
      m[5]*(2.0*m[4]*m[7]-m[2]*m[5])-m[3]*m[7]*m[7])+f*f+
    -m[1]*(m[2]*m[6]*m[6]+m[8]*(m[3]*m[8]-2.0*m[4]*m[6]))+
    m[9]*(
        2.0*(m[2]*m[5]*m[6]+m[3]*m[7]*m[8]-m[4]*(m[6]*m[7]+m[5]*m[8]))+
        m[9]*(m[4]*m[4]-m[2]*m[3])); }

inline int form4_inv(form4 &m_inv, form4 &m) {
  double det = form4_det(m);
  if (fzero(det)) return 0;
  m_inv[0] = (m[1]*(m[2]*m[3]-m[4]*m[4])+
      m[5]*(2.0*m[4]*m[7]-m[2]*m[5])-m[3]*m[7]*m[7])/det;
  m_inv[1] = (m[0]*(m[2]*m[3]-m[4]*m[4])+
      m[6]*(2.0*m[4]*m[8]-m[2]*m[6])-m[3]*m[8]*m[8])/det;
  m_inv[2] = (m[0]*(m[1]*m[3]-m[5]*m[5])+
      m[6]*(2.0*m[5]*m[9]-m[1]*m[6])-m[3]*m[9]*m[9])/det;
  m_inv[3] = (m[0]*(m[1]*m[2]-m[7]*m[7])+
      m[8]*(2.0*m[7]*m[9]-m[1]*m[8])-m[2]*m[9]*m[9])/det;
  m_inv[4] = (m[0]*(m[5]*m[7]-m[1]*m[4])+m[1]*m[6]*m[8]+
      m[9]*(m[4]*m[9]-m[6]*m[7]-m[5]*m[8]))/det;
  m_inv[5] = (m[0]*(m[4]*m[7]-m[2]*m[5])+m[2]*m[6]*m[9]+
      m[8]*(m[5]*m[8]-m[6]*m[7]-m[4]*m[9]))/det;
  m_inv[6] = (m[1]*(m[4]*m[8]-m[2]*m[6])+m[2]*m[5]*m[9]+
      m[7]*(m[6]*m[7]-m[5]*m[8]-m[4]*m[9]))/det;
  m_inv[7] = (m[0]*(m[4]*m[5]-m[3]*m[7])+m[3]*m[8]*m[9]+
      m[6]*(m[6]*m[7]-m[5]*m[8]-m[4]*m[9]))/det;
  m_inv[8] = (m[1]*(m[4]*m[6]-m[3]*m[8])+m[3]*m[7]*m[9]+
      m[5]*(m[5]*m[8]-m[6]*m[7]-m[4]*m[9]))/det;
  m_inv[9] = (m[2]*(m[5]*m[6]-m[3]*m[9])+m[3]*m[7]*m[8]+
      m[4]*(m[4]*m[9]-m[6]*m[7]-m[5]*m[8]))/det;
  return 1; }

inline void form4_vec4_dot(vector4 &dest, form4 &m) {
  vector4 a;
  memcpy(a, dest, sizeof(vector4));
  dest[0] = m[0]*a[0]+m[9]*a[1]+m[7]*a[2]+m[6]*a[3];
  dest[1] = m[9]*a[0]+m[1]*a[1]+m[6]*a[2]+m[5]*a[3];
  dest[2] = m[8]*a[0]+m[7]*a[1]+m[2]*a[2]+m[4]*a[3];
  dest[3] = m[6]*a[0]+m[5]*a[1]+m[4]*a[2]+m[3]*a[3]; }

// square regression: y = eqn[0] + eqn[1]*x + eqn[2]*x*x

inline int regress2(vector &eqn, int order, double *x, double *y, int n) {
  form xtx = FORM_NULL, xtx_inv;
  vector xty = VECTOR_NULL;
  double xn, xi, yi;
  int i;

  vec_null(eqn);
  xtx[0] = n;
  if ((order = order%2)<0) order = -order;                // max: quad regress
  if (order<1) xtx[1] = 1.0;
  if (order<2) xtx[2] = 1.0;
  for (i=0; i<n; ++i) {
    xty[0] += (yi = y[i]);
    if (order<1) continue;
    xty[1] += yi*(xi = xn = x[i]); xtx[5] += xn; xtx[1] += (xn *= xi);
    if (order<2) continue;
    xty[2] += yi*xn; xtx[3] += (xn *= xi); xtx[2] += xn*xi;
  }
  if (order>1) xtx[4] = xtx[2];
  if (!form_inv(xtx_inv, xtx)) return 0;
  memcpy(eqn, xty, sizeof(vector));
  form_vec_dot(eqn, xtx_inv);
  return 1; }

// cubic regression: y = eqn[0] + eqn[1]*x + eqn[2]*x*x + eqn[3]*x*x*x

inline int regress3(vector4 &eqn, int order, double *x, double *y, int n) {
  form4 xtx = FORM4_NULL, xtx_inv;
  vector4 xty = VECTOR4_NULL;
  double xn, xi, yi;
  int i;

  vec4_null(eqn);
  xtx[0] = n;
  if ((order = order%3)<0) order = -order;                // max: cubic regress
  if (order<1) xtx[1] = 1.0;
  if (order<2) xtx[2] = 1.0;
  if (order<3) xtx[3] = 1.0;
  for (i=0; i<n; ++i) {
    xty[0] += (yi = y[i]);
    if (order<1) continue;
    xty[1] += yi*(xi = xn = x[i]); xtx[9] += xn; xtx[1] += (xn *= xi);
    if (order<2) continue;
    xty[2] += yi*xn; xtx[7] += (xn *= xi); xtx[2] += xn*xi;
    if (order<3) continue;
    xty[3] += yi*xn; xtx[4] += (xn *= xi*xi); xtx[3] += xn*xi;
  }
  if (order>1) xtx[8] = xtx[1];
  if (order>2) { xtx[6] = xtx[7]; xtx[5] = xtx[2]; }
  if (!form4_inv(xtx_inv, xtx)) return 0;
  memcpy(eqn, xty, sizeof(vector4));
  form4_vec4_dot(eqn, xtx_inv);
  return 1; }

}

#endif
