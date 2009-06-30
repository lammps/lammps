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
   Contributing author: Mike Brown (SNL)
------------------------------------------------------------------------- */

#ifndef MATH_EXTRA_H
#define MATH_EXTRA_H

#include "math.h"
#include "stdio.h"
#include "string.h"
#include "error.h"

namespace MathExtra {

  // 3 vector operations

  inline void normalize3(const double *v, double *ans);
  inline void snormalize3(const double, const double *v, double *ans);
  inline double dot3(const double *v1, const double *v2);
  inline void cross3(const double *v1, const double *v2, double *ans);

  // 3x3 matrix operations

  inline double det3(const double mat[3][3]);
  inline void diag_times3(const double *diagonal, const double mat[3][3],
                          double ans[3][3]);
  inline void plus3(const double m[3][3], const double m2[3][3],
                    double ans[3][3]);
  inline void times3(const double m[3][3], const double m2[3][3],
                     double ans[3][3]);
  inline void transpose_times3(const double mat1[3][3], 
                               const double mat2[3][3],
                               double ans[3][3]);
  inline void invert3(const double mat[3][3], double ans[3][3]);
  inline void times_column3(const double mat[3][3], const double*vec,
                            double *ans);
  inline void transpose_times_column3(const double mat[3][3],const double*vec,
                                      double *ans);
  inline void transpose_times_diag3(const double mat[3][3],const double*vec,
                                    double ans[3][3]);
  inline void row_times3(const double *v, const double m[3][3],
                         double *ans);
  inline void mldivide3(const double mat[3][3], const double *vec,
                        double *ans, LAMMPS_NS::Error *error);
  inline void write3(const double mat[3][3]);
      
  // quaternion operations
  
  inline void normalize4(double *quat);
  inline void quat_to_mat(const double *quat, double mat[3][3]);
  inline void quat_to_mat_trans(const double *quat, double mat[3][3]);
  inline void axisangle_to_quat(const double *v, const double angle,
                                double *quat);
  inline void multiply_quat_quat(const double *one, const double *two,
                                 double *ans);
  inline void multiply_vec_quat(const double *one, const double *two,
                                double *ans);

  // rotation operations
  
  inline void rotation_generator_x(const double m[3][3], double ans[3][3]);
  inline void rotation_generator_y(const double m[3][3], double ans[3][3]);
  inline void rotation_generator_z(const double m[3][3], double ans[3][3]);
  
  // shape matrix operations

  inline void multiply_shape_shape(const double *one, const double *two,
                                   double *ans);
};

/* ----------------------------------------------------------------------
   normalize a vector
------------------------------------------------------------------------- */

void MathExtra::normalize3(const double *v, double *ans)
{
  double scale = 1.0/sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  ans[0] = v[0]*scale;
  ans[1] = v[1]*scale;
  ans[2] = v[2]*scale;
}

/* ----------------------------------------------------------------------
   scale a vector to length
------------------------------------------------------------------------- */

void MathExtra::snormalize3(const double length, const double *v, double *ans)
{
  double scale = length/sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  ans[0] = v[0]*scale;
  ans[1] = v[1]*scale;
  ans[2] = v[2]*scale;
}

/* ----------------------------------------------------------------------
   dot product of 2 vectors
------------------------------------------------------------------------- */

double MathExtra::dot3(const double *v1, const double *v2)
{
  return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

/* ----------------------------------------------------------------------
   cross product of 2 vectors
------------------------------------------------------------------------- */

void MathExtra::cross3(const double *v1, const double *v2, double *ans)
{
  ans[0] = v1[1]*v2[2]-v1[2]*v2[1];
  ans[1] = v1[2]*v2[0]-v1[0]*v2[2];
  ans[2] = v1[0]*v2[1]-v1[1]*v2[0];
}

/* ----------------------------------------------------------------------
   determinant of a matrix
------------------------------------------------------------------------- */

double MathExtra::det3(const double m[3][3])
{
  double ans = m[0][0]*m[1][1]*m[2][2] - m[0][0]*m[1][2]*m[2][1] - 
    m[1][0]*m[0][1]*m[2][2] + m[1][0]*m[0][2]*m[2][1] + 
    m[2][0]*m[0][1]*m[1][2] - m[2][0]*m[0][2]*m[1][1];
  return ans;
}

/* ----------------------------------------------------------------------
   diagonal matrix times a full matrix
------------------------------------------------------------------------- */

void MathExtra::diag_times3(const double *d, const double m[3][3],
			    double ans[3][3])
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

void MathExtra::plus3(const double m[3][3], const double m2[3][3],
		      double ans[3][3])
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

void MathExtra::times3(const double m[3][3], const double m2[3][3],
                       double ans[3][3])
{
  ans[0][0] = m[0][0]*m2[0][0]+m[0][1]*m2[1][0]+m[0][2]*m2[2][0];
  ans[0][1] = m[0][0]*m2[0][1]+m[0][1]*m2[1][1]+m[0][2]*m2[2][1];
  ans[0][2] = m[0][0]*m2[0][2]+m[0][1]*m2[1][2]+m[0][2]*m2[2][2];
  ans[1][0] = m[1][0]*m2[0][0]+m[1][1]*m2[1][0]+m[1][2]*m2[2][0];
  ans[1][1] = m[1][0]*m2[0][1]+m[1][1]*m2[1][1]+m[1][2]*m2[2][1];
  ans[1][2] = m[1][0]*m2[0][2]+m[1][1]*m2[1][2]+m[1][2]*m2[2][2];
  ans[2][0] = m[2][0]*m2[0][0]+m[2][1]*m2[1][0]+m[2][2]*m2[2][0];
  ans[2][1] = m[2][0]*m2[0][1]+m[2][1]*m2[1][1]+m[2][2]*m2[2][1];
  ans[2][2] = m[2][0]*m2[0][2]+m[2][1]*m2[1][2]+m[2][2]*m2[2][2];
}

/* ----------------------------------------------------------------------
   multiply the transpose of mat1 times mat2
------------------------------------------------------------------------- */

void MathExtra::transpose_times3(const double m[3][3], const double m2[3][3],
                                 double ans[3][3])
{
  ans[0][0] = m[0][0]*m2[0][0]+m[1][0]*m2[1][0]+m[2][0]*m2[2][0];
  ans[0][1] = m[0][0]*m2[0][1]+m[1][0]*m2[1][1]+m[2][0]*m2[2][1];
  ans[0][2] = m[0][0]*m2[0][2]+m[1][0]*m2[1][2]+m[2][0]*m2[2][2];
  ans[1][0] = m[0][1]*m2[0][0]+m[1][1]*m2[1][0]+m[2][1]*m2[2][0];
  ans[1][1] = m[0][1]*m2[0][1]+m[1][1]*m2[1][1]+m[2][1]*m2[2][1];
  ans[1][2] = m[0][1]*m2[0][2]+m[1][1]*m2[1][2]+m[2][1]*m2[2][2];
  ans[2][0] = m[0][2]*m2[0][0]+m[1][2]*m2[1][0]+m[2][2]*m2[2][0];
  ans[2][1] = m[0][2]*m2[0][1]+m[1][2]*m2[1][1]+m[2][2]*m2[2][1];
  ans[2][2] = m[0][2]*m2[0][2]+m[1][2]*m2[1][2]+m[2][2]*m2[2][2];
}

/* ----------------------------------------------------------------------
   invert a matrix
   does NOT checks for singular or badly scaled matrix
------------------------------------------------------------------------- */

void MathExtra::invert3(const double m[3][3], double ans[3][3])
{
  double den = m[0][0]*m[1][1]*m[2][2]-m[0][0]*m[1][2]*m[2][1];
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

void MathExtra::times_column3(const double m[3][3], const double *v,
			      double *ans) 
{
  ans[0] = m[0][0]*v[0]+m[0][1]*v[1]+m[0][2]*v[2];
  ans[1] = m[1][0]*v[0]+m[1][1]*v[1]+m[1][2]*v[2];
  ans[2] = m[2][0]*v[0]+m[2][1]*v[1]+m[2][2]*v[2];
}

/* ----------------------------------------------------------------------
   transposed matrix times vector
------------------------------------------------------------------------- */

void MathExtra::transpose_times_column3(const double m[3][3], const double *v,
					double *ans)
{
  ans[0] = m[0][0]*v[0]+m[1][0]*v[1]+m[2][0]*v[2];
  ans[1] = m[0][1]*v[0]+m[1][1]*v[1]+m[2][1]*v[2];
  ans[2] = m[0][2]*v[0]+m[1][2]*v[1]+m[2][2]*v[2];
}

/* ----------------------------------------------------------------------
   transposed matrix times diagonal matrix
------------------------------------------------------------------------- */

void MathExtra::transpose_times_diag3(const double m[3][3],
                                      const double *d, double ans[3][3])
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

void MathExtra::row_times3(const double *v, const double m[3][3],
			   double *ans)
{
  ans[0] = m[0][0]*v[0]+v[1]*m[1][0]+v[2]*m[2][0];
  ans[1] = v[0]*m[0][1]+m[1][1]*v[1]+v[2]*m[2][1];
  ans[2] = v[0]*m[0][2]+v[1]*m[1][2]+m[2][2]*v[2];
}

/* ----------------------------------------------------------------------
   solve Ax = b or M ans = v
   use gaussian elimination & partial pivoting on matrix
------------------------------------------------------------------------- */

void MathExtra::mldivide3(const double m[3][3], const double *v, double *ans,
			  LAMMPS_NS::Error *error)
{
  // create augmented matrix for pivoting

  double aug[3][4];
  for (unsigned i = 0; i < 3; i++) {
    aug[i][3] = v[i];
    for (unsigned j = 0; j < 3; j++) aug[i][j] = m[i][j];
  }

  for (unsigned i = 0; i < 2; i++) {
    unsigned p = i;
    for (unsigned j = i+1; j < 3; j++) {
      if (fabs(aug[j][i]) > fabs(aug[i][i])) {
        double tempv[4];
        memcpy(tempv,aug[i],4*sizeof(double));
        memcpy(aug[i],aug[j],4*sizeof(double));
        memcpy(aug[j],tempv,4*sizeof(double));
      }
    }

    while (aug[p][i] == 0.0 && p < 3) p++;

    if (p == 3) error->all("Bad matrix inversion in MathExtra::mldivide3");
    else
      if (p != i) {
        double tempv[4];
        memcpy(tempv,aug[i],4*sizeof(double));
        memcpy(aug[i],aug[p],4*sizeof(double));
        memcpy(aug[p],tempv,4*sizeof(double));
      }

    for (unsigned j = i+1; j < 3; j++) {
      double m = aug[j][i]/aug[i][i];
      for (unsigned k=i+1; k<4; k++) aug[j][k]-=m*aug[i][k];
    }
  }

  if (aug[2][2] == 0.0)
    error->all("Bad matrix inversion in MathExtra::mldivide3");
  
  // back substitution

  ans[2] = aug[2][3]/aug[2][2];
  for (int i = 1; i >= 0; i--) {
    double sumax = 0.0;
    for (unsigned j = i+1; j < 3; j++) sumax += aug[i][j]*ans[j];
    ans[i] = (aug[i][3]-sumax) / aug[i][i];
  }
}

/* ----------------------------------------------------------------------
   output a matrix
------------------------------------------------------------------------- */

void MathExtra::write3(const double mat[3][3])
{
  for (unsigned i = 0; i < 3; i++) {
    for (unsigned j = 0; j < 3; j++) printf("%g ",mat[i][j]);
    printf("\n");
  }
}

/* ----------------------------------------------------------------------
   normalize a quaternion
------------------------------------------------------------------------- */

void MathExtra::normalize4(double *quat)
{
  double scale = 1.0/sqrt(quat[0]*quat[0]+quat[1]*quat[1] +
			  quat[2]*quat[2]+quat[3]*quat[3]);
  quat[0] *= scale;
  quat[1] *= scale;
  quat[2] *= scale;
  quat[3] *= scale;
}

/* ----------------------------------------------------------------------
   compute rotation matrix from quaternion
   quat = [w i j k]
------------------------------------------------------------------------- */

void MathExtra::quat_to_mat(const double *quat, double mat[3][3])
{
  double w2 = quat[0]*quat[0];
  double i2 = quat[1]*quat[1];
  double j2 = quat[2]*quat[2];
  double k2 = quat[3]*quat[3];
  double twoij = 2.0*quat[1]*quat[2];
  double twoik = 2.0*quat[1]*quat[3];
  double twojk = 2.0*quat[2]*quat[3];
  double twoiw = 2.0*quat[1]*quat[0];
  double twojw = 2.0*quat[2]*quat[0];
  double twokw = 2.0*quat[3]*quat[0];

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
   compute rotation matrix from quaternion conjugate
   quat = [w i j k]
------------------------------------------------------------------------- */

void MathExtra::quat_to_mat_trans(const double *quat, double mat[3][3])
{
  double w2 = quat[0]*quat[0];
  double i2 = quat[1]*quat[1];
  double j2 = quat[2]*quat[2];
  double k2 = quat[3]*quat[3];
  double twoij = 2.0*quat[1]*quat[2];
  double twoik = 2.0*quat[1]*quat[3];
  double twojk = 2.0*quat[2]*quat[3];
  double twoiw = 2.0*quat[1]*quat[0];
  double twojw = 2.0*quat[2]*quat[0];
  double twokw = 2.0*quat[3]*quat[0];

  mat[0][0] = w2+i2-j2-k2;
  mat[1][0] = twoij-twokw;
  mat[2][0] = twojw+twoik;

  mat[0][1] = twoij+twokw;
  mat[1][1] = w2-i2+j2-k2;
  mat[2][1] = twojk-twoiw;
	
  mat[0][2] = twoik-twojw;
  mat[1][2] = twojk+twoiw;
  mat[2][2] = w2-i2-j2+k2;
}

/* ----------------------------------------------------------------------
   compute quaternion from axis-angle rotation
   v MUST be a unit vector
------------------------------------------------------------------------- */

void MathExtra::axisangle_to_quat(const double *v, const double angle,
				  double *quat)
{
  double halfa = 0.5*angle;
  double sina = sin(halfa);
  quat[0] = cos(halfa);
  quat[1] = v[0]*sina;
  quat[2] = v[1]*sina;
  quat[3] = v[2]*sina;
}

/* ----------------------------------------------------------------------
   multiply 2 quaternions
   effectively concatenates rotations
   NOT a commutative operation
------------------------------------------------------------------------- */

void MathExtra::multiply_quat_quat(const double *one, const double *two,
				   double *ans)
{
  ans[0] = one[0]*two[0]-one[1]*two[1]-one[2]*two[2]-one[3]*two[3];
  ans[1] = one[0]*two[1]+one[1]*two[0]+one[2]*two[3]-one[3]*two[2];
  ans[2] = one[0]*two[2]-one[1]*two[3]+one[2]*two[0]+one[3]*two[1];
  ans[3] = one[0]*two[3]+one[1]*two[2]-one[2]*two[1]+one[3]*two[0];
}

/* ----------------------------------------------------------------------
   multiply 3 vector times quaternion
   3 vector one is treated as quaternion (0,one)
------------------------------------------------------------------------- */

void MathExtra::multiply_vec_quat(const double *one, const double *two,
				  double *ans)
{
  ans[0] = -one[0]*two[1]-one[1]*two[2]-one[2]*two[3];
  ans[1] =  one[0]*two[0]+one[1]*two[3]-one[2]*two[2];
  ans[2] = -one[0]*two[3]+one[1]*two[0]+one[2]*two[1];
  ans[3] =  one[0]*two[2]-one[1]*two[1]+one[2]*two[0];
}

/* ----------------------------------------------------------------------
   multiply 2 shape matrices
   upper-triangular 3x3, stored as 6-vector in Voigt notation
------------------------------------------------------------------------- */

void MathExtra::multiply_shape_shape(const double *one, const double *two,
                                     double *ans)
{
  ans[0] = one[0]*two[0];
  ans[1] = one[1]*two[1];
  ans[2] = one[2]*two[2];
  ans[3] = one[1]*two[3] + one[3]*two[2];
  ans[4] = one[0]*two[4] + one[5]*two[3] + one[4]*two[2];
  ans[5] = one[0]*two[5] + one[5]*two[1];
}

/* ----------------------------------------------------------------------
   Apply principal rotation generator about x to rotation matrix m
------------------------------------------------------------------------- */

void MathExtra::rotation_generator_x(const double m[3][3], double ans[3][3])
{
  ans[0][0]=0;
  ans[0][1]=-m[0][2];
  ans[0][2]=m[0][1];
  ans[1][0]=0;
  ans[1][1]=-m[1][2];
  ans[1][2]=m[1][1];
  ans[2][0]=0;
  ans[2][1]=-m[2][2];
  ans[2][2]=m[2][1];
}

/* ----------------------------------------------------------------------
   Apply principal rotation generator about y to rotation matrix m
------------------------------------------------------------------------- */

void MathExtra::rotation_generator_y(const double m[3][3], double ans[3][3])
{
  ans[0][0]=m[0][2];
  ans[0][1]=0;
  ans[0][2]=-m[0][0];
  ans[1][0]=m[1][2];
  ans[1][1]=0;
  ans[1][2]=-m[1][0];
  ans[2][0]=m[2][2];
  ans[2][1]=0;
  ans[2][2]=-m[2][0];
}

/* ----------------------------------------------------------------------
   Apply principal rotation generator about z to rotation matrix m
------------------------------------------------------------------------- */

void MathExtra::rotation_generator_z(const double m[3][3], double ans[3][3])
{
  ans[0][0]=-m[0][1];
  ans[0][1]=m[0][0];
  ans[0][2]=0;
  ans[1][0]=-m[1][1];
  ans[1][1]=m[1][0];
  ans[1][2]=0;
  ans[2][0]=-m[2][1];
  ans[2][1]=m[2][0];
  ans[2][2]=0;
}
 
#endif
