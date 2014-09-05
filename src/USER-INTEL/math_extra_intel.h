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
   Contributing author: W. Michael Brown (Intel)
------------------------------------------------------------------------- */

#ifndef LMP_MATH_EXTRA_INTEL_H
#define LMP_MATH_EXTRA_INTEL_H

#define ME_quat_to_mat_trans(quat, mat)		\
{						\
  flt_t quat_w = quat.w;			\
  flt_t quat_i = quat.i;			\
  flt_t quat_j = quat.j;			\
  flt_t quat_k = quat.k;			\
  flt_t w2 = quat_w * quat_w;			\
  flt_t i2 = quat_i * quat_i;			\
  flt_t j2 = quat_j * quat_j;			\
  flt_t k2 = quat_k * quat_k;			\
  flt_t twoij = (flt_t)2.0 * quat_i * quat_j;	\
  flt_t twoik = (flt_t)2.0 * quat_i * quat_k;	\
  flt_t twojk = (flt_t)2.0 * quat_j * quat_k;	\
  flt_t twoiw = (flt_t)2.0 * quat_i * quat_w;	\
  flt_t twojw = (flt_t)2.0 * quat_j * quat_w;	\
  flt_t twokw = (flt_t)2.0 * quat_k * quat_w;	\
  						\
  mat##_0 = w2 + i2 - j2 - k2;			\
  mat##_3 = twoij - twokw;			\
  mat##_6 = twojw + twoik;			\
  						\
  mat##_1 = twoij + twokw;			\
  mat##_4 = w2 - i2 + j2 - k2;			\
  mat##_7 = twojk - twoiw;			\
  						\
  mat##_2 = twoik - twojw;			\
  mat##_5 = twojk + twoiw;			\
  mat##_8 = w2 - i2 - j2 + k2;			\
}

/* ----------------------------------------------------------------------
   diagonal matrix times a full matrix
------------------------------------------------------------------------- */

#define ME_diag_times3(d, m, ans)			\
  {							\
  ans##_0 = d[0] * m##_0;				\
  ans##_1 = d[0] * m##_1;				\
  ans##_2 = d[0] * m##_2;				\
  ans##_3 = d[1] * m##_3;				\
  ans##_4 = d[1] * m##_4;				\
  ans##_5 = d[1] * m##_5;				\
  ans##_6 = d[2] * m##_6;				\
  ans##_7 = d[2] * m##_7;				\
  ans##_8 = d[2] * m##_8;				\
}

#define ME_diag_times3a(d, m, ans)			\
  {							\
  ans##_0 = d##_0 * m##_0;				\
  ans##_1 = d##_0 * m##_1;				\
  ans##_2 = d##_0 * m##_2;				\
  ans##_3 = d##_1 * m##_3;				\
  ans##_4 = d##_1 * m##_4;				\
  ans##_5 = d##_1 * m##_5;				\
  ans##_6 = d##_2 * m##_6;				\
  ans##_7 = d##_2 * m##_7;				\
  ans##_8 = d##_2 * m##_8;				\
}

/* ----------------------------------------------------------------------
   multiply the transpose of mat1 times mat2
------------------------------------------------------------------------- */

#define ME_transpose_times3(m1, m2, ans)                	\
{								\
  ans##_0 = m1##_0*m2##_0 + m1##_3*m2##_3 + m1##_6*m2##_6;	\
  ans##_1 = m1##_0*m2##_1 + m1##_3*m2##_4 + m1##_6*m2##_7;	\
  ans##_2 = m1##_0*m2##_2 + m1##_3*m2##_5 + m1##_6*m2##_8;	\
  ans##_3 = m1##_1*m2##_0 + m1##_4*m2##_3 + m1##_7*m2##_6;	\
  ans##_4 = m1##_1*m2##_1 + m1##_4*m2##_4 + m1##_7*m2##_7;	\
  ans##_5 = m1##_1*m2##_2 + m1##_4*m2##_5 + m1##_7*m2##_8;	\
  ans##_6 = m1##_2*m2##_0 + m1##_5*m2##_3 + m1##_8*m2##_6;	\
  ans##_7 = m1##_2*m2##_1 + m1##_5*m2##_4 + m1##_8*m2##_7;	\
  ans##_8 = m1##_2*m2##_2 + m1##_5*m2##_5 + m1##_8*m2##_8;	\
}

/* ----------------------------------------------------------------------
   normalize a vector, return in ans
------------------------------------------------------------------------- */

#define ME_normalize3(v0, v1, v2, ans)	        	\
{							\
  flt_t scale = (flt_t)1.0 / sqrt(v0*v0+v1*v1+v2*v2);	\
  ans##_0 = v0 * scale;					\
  ans##_1 = v1 * scale;					\
  ans##_2 = v2 * scale;					\
}

/* ----------------------------------------------------------------------
   add two matrices
------------------------------------------------------------------------- */

#define ME_plus3(m1, m2, ans)			\
{						\
  ans##_0 = m1##_0 + m2##_0;			\
  ans##_1 = m1##_1 + m2##_1;			\
  ans##_2 = m1##_2 + m2##_2;			\
  ans##_3 = m1##_3 + m2##_3;			\
  ans##_4 = m1##_4 + m2##_4;			\
  ans##_5 = m1##_5 + m2##_5;			\
  ans##_6 = m1##_6 + m2##_6;			\
  ans##_7 = m1##_7 + m2##_7;			\
  ans##_8 = m1##_8 + m2##_8;			\
}

/* ----------------------------------------------------------------------
   dot product of 2 vectors
------------------------------------------------------------------------- */

#define ME_dot3(v1, v2) \
  (v1##_0*v2##_0 + v1##_1 * v2##_1 + v1##_2 * v2##_2)

/* ----------------------------------------------------------------------
   determinant of a matrix
------------------------------------------------------------------------- */

#define ME_det3(m)				    \
  ( m##_0 * m##_4 * m##_8 - m##_0 * m##_5 * m##_7 - \
    m##_3 * m##_1 * m##_8 + m##_3 * m##_2 * m##_7 + \
    m##_6 * m##_1 * m##_5 - m##_6 * m##_2 * m##_4 )

/* ----------------------------------------------------------------------
   row vector times matrix
------------------------------------------------------------------------- */

#define ME_vecmat(v, m, ans)				    \
{							    \
  ans##_0 = v##_0 * m##_0 + v##_1 * m##_3 + v##_2 * m##_6;  \
  ans##_1 = v##_0 * m##_1 + v##_1 * m##_4 + v##_2 * m##_7;  \
  ans##_2 = v##_0 * m##_2 + v##_1 * m##_5 + v##_2 * m##_8;  \
}

/* ----------------------------------------------------------------------
   cross product of 2 vectors
------------------------------------------------------------------------- */

#define ME_cross3(v1, v2, ans)			\
{						\
  ans##_0 = v1##_1 * v2##_2 - v1##_2 * v2##_1;	\
  ans##_1 = v1##_2 * v2##_0 - v1##_0 * v2##_2;	\
  ans##_2 = v1##_0 * v2##_1 - v1##_1 * v2##_0;	\
}

/* ----------------------------------------------------------------------
   cross product of 2 vectors
------------------------------------------------------------------------- */

#define ME_mv0_cross3(m1, v2, ans)		\
{						\
  ans##_0 = m1##_1 * v2##_2 - m1##_2 * v2##_1;	\
  ans##_1 = m1##_2 * v2##_0 - m1##_0 * v2##_2;	\
  ans##_2 = m1##_0 * v2##_1 - m1##_1 * v2##_0;	\
}

#define ME_mv1_cross3(m1, v2, ans)		\
{						\
  ans##_0 = m1##_4 * v2##_2 - m1##_5 * v2##_1;	\
  ans##_1 = m1##_5 * v2##_0 - m1##_3 * v2##_2;	\
  ans##_2 = m1##_3 * v2##_1 - m1##_4 * v2##_0;	\
}

#define ME_mv2_cross3(m1, v2, ans)		\
{						\
  ans##_0 = m1##_7 * v2##_2 - m1##_8 * v2##_1;	\
  ans##_1 = m1##_8 * v2##_0 - m1##_6 * v2##_2;	\
  ans##_2 = m1##_6 * v2##_1 - m1##_7 * v2##_0;	\
}


#define ME_compute_eta_torque(m1, m2, s1, ans)                              \
{								            \
  flt_t den = m1##_3*m1##_2*m1##_7-m1##_0*m1##_5*m1##_7-		    \
    m1##_2*m1##_6*m1##_4+m1##_1*m1##_6*m1##_5-				    \
    m1##_3*m1##_1*m1##_8+m1##_0*m1##_4*m1##_8;				    \
  den = (flt_t)1.0 / den;						    \
									    \
  ans##_0 = s1##_0*(m1##_5*m1##_1*m2##_2+(flt_t)2.0*m1##_4*m1##_8*m2##_0-   \
 		   m1##_4*m2##_2*m1##_2-(flt_t)2.0*m1##_5*m2##_0*m1##_7+    \
		   m2##_1*m1##_2*m1##_7-m2##_1*m1##_1*m1##_8-	 	    \
		   m1##_3*m1##_8*m2##_1+m1##_6*m1##_5*m2##_1+		    \
		   m1##_3*m2##_2*m1##_7-m2##_2*m1##_6*m1##_4)*den;	    \
									    \
  ans##_1 = s1##_0*(m1##_2*m2##_0*m1##_7-m1##_8*m2##_0*m1##_1+		    \
		   (flt_t)2.0*m1##_0*m1##_8*m2##_1-m1##_0*m2##_2*m1##_5-    \
		   (flt_t)2.0*m1##_6*m1##_2*m2##_1+m2##_2*m1##_3*m1##_2-    \
		   m1##_8*m1##_3*m2##_0+m1##_6*m2##_0*m1##_5+		    \
		   m1##_6*m2##_2*m1##_1-m2##_2*m1##_0*m1##_7)*den;	    \
  									    \
  ans##_2 = s1##_0*(m1##_1*m1##_5*m2##_0-m1##_2*m2##_0*m1##_4-              \
		   m1##_0*m1##_5*m2##_1+m1##_3*m1##_2*m2##_1-		    \
		   m2##_1*m1##_0*m1##_7-m1##_6*m1##_4*m2##_0+		    \
		   (flt_t)2.0*m1##_4*m1##_0*m2##_2-                         \
		   (flt_t)2.0*m1##_3*m2##_2*m1##_1+			    \
		   m1##_3*m1##_7*m2##_0+m1##_6*m2##_1*m1##_1)*den;	    \
									    \
  ans##_3 = s1##_1*(-m1##_4*m2##_5*m1##_2+(flt_t)2.0*m1##_4*m1##_8*m2##_3+  \
		   m1##_5*m1##_1*m2##_5-(flt_t)2.0*m1##_5*m2##_3*m1##_7+    \
		   m2##_4*m1##_2*m1##_7-m2##_4*m1##_1*m1##_8-		    \
		   m1##_3*m1##_8*m2##_4+m1##_6*m1##_5*m2##_4- 		    \
		   m2##_5*m1##_6*m1##_4+m1##_3*m2##_5*m1##_7)*den;	    \
									    \
  ans##_4 = s1##_1*(m1##_2*m2##_3*m1##_7-m1##_1*m1##_8*m2##_3+		    \
		   (flt_t)2.0*m1##_8*m1##_0*m2##_4-m2##_5*m1##_0*m1##_5-    \
		   (flt_t)2.0*m1##_6*m2##_4*m1##_2-m1##_3*m1##_8*m2##_3+    \
		   m1##_6*m1##_5*m2##_3+m1##_3*m2##_5*m1##_2-		    \
		   m1##_0*m2##_5*m1##_7+m2##_5*m1##_1*m1##_6)*den;	    \
									    \
  ans##_5 = s1##_1*(m1##_1*m1##_5*m2##_3-m1##_2*m2##_3*m1##_4-		    \
		   m1##_0*m1##_5*m2##_4+m1##_3*m1##_2*m2##_4+		    \
		   (flt_t)2.0*m1##_4*m1##_0*m2##_5-m1##_0*m2##_4*m1##_7+    \
		   m1##_1*m1##_6*m2##_4-m2##_3*m1##_6*m1##_4-		    \
		   (flt_t)2.0*m1##_3*m1##_1*m2##_5+m1##_3*m2##_3*m1##_7)*   \
    den;							   	    \
									    \
  ans##_6 = s1##_2*(-m1##_4*m1##_2*m2##_8+m1##_1*m1##_5*m2##_8+  	    \
		   (flt_t)2.0*m1##_4*m2##_6*m1##_8-m1##_1*m2##_7*m1##_8+    \
		   m1##_2*m1##_7*m2##_7-(flt_t)2.0*m2##_6*m1##_7*m1##_5-    \
		   m1##_3*m2##_7*m1##_8+m1##_5*m1##_6*m2##_7-		    \
		   m1##_4*m1##_6*m2##_8+m1##_7*m1##_3*m2##_8)*den;	    \
									    \
  ans##_7 = s1##_2*-(m1##_1*m1##_8*m2##_6-m1##_2*m2##_6*m1##_7-		    \
		    (flt_t)2.0*m2##_7*m1##_0*m1##_8+m1##_5*m2##_8*m1##_0+   \
		    (flt_t)2.0*m2##_7*m1##_2*m1##_6+m1##_3*m2##_6*m1##_8-   \
		    m1##_3*m1##_2*m2##_8-m1##_5*m1##_6*m2##_6+		    \
		    m1##_0*m2##_8*m1##_7-m2##_8*m1##_1*m1##_6)*den;	    \
									    \
  ans##_8 = s1##_2*(m1##_1*m1##_5*m2##_6-m1##_2*m2##_6*m1##_4-		    \
		   m1##_0*m1##_5*m2##_7+m1##_3*m1##_2*m2##_7-		    \
		   m1##_4*m1##_6*m2##_6-m1##_7*m2##_7*m1##_0+		    \
		   (flt_t)2.0*m1##_4*m2##_8*m1##_0+m1##_7*m1##_3*m2##_6+    \
                    m1##_6*m1##_1*m2##_7-(flt_t)2.0*m2##_8*m1##_3*m1##_1)*  \
    den;								    \
}

#define ME_vcopy4(dst,src)			\
  dst##_0 = src##_0;				\
  dst##_1 = src##_1;				\
  dst##_2 = src##_2;				\
  dst##_3 = src##_3;

#define ME_mldivide3(m1, v_0, v_1, v_2, ans, error)	\
{							\
  flt_t aug_0, aug_1, aug_2, aug_3, aug_4, aug_5;	\
  flt_t aug_6, aug_7, aug_8, aug_9, aug_10, aug_11, t;	\
							\
  aug_3 = v_0;						\
  aug_0 = m1##_0;					\
  aug_1 = m1##_1;					\
  aug_2 = m1##_2;					\
  aug_7 = v_1;						\
  aug_4 = m1##_3;					\
  aug_5 = m1##_4;					\
  aug_6 = m1##_5;					\
  aug_11 = v_2;						\
  aug_8 = m1##_6;					\
  aug_9 = m1##_7;					\
  aug_10 = m1##_8;					\
							\
  if (fabs(aug_4) > fabs(aug_0)) {			\
    flt_t swapt;					\
    swapt = aug_0; aug_0 = aug_4; aug_4 = swapt;	\
    swapt = aug_1; aug_1 = aug_5; aug_5 = swapt;	\
    swapt = aug_2; aug_2 = aug_6; aug_6 = swapt;	\
    swapt = aug_3; aug_3 = aug_7; aug_7 = swapt;	\
  }							\
  if (fabs(aug_8) > fabs(aug_0)) {			\
    flt_t swapt;					\
    swapt = aug_0; aug_0 = aug_8; aug_8 = swapt;	\
    swapt = aug_1; aug_1 = aug_9; aug_9 = swapt;        \
    swapt = aug_2; aug_2 = aug_10; aug_10 = swapt;      \
    swapt = aug_3; aug_3 = aug_11; aug_11 = swapt;      \
  }							\
							\
  if (aug_0 != (flt_t)0.0) {				\
  } else if (aug_4 != (flt_t)0.0) {			\
    flt_t swapt;					\
    swapt = aug_0; aug_0 = aug_4; aug_4 = swapt;	\
    swapt = aug_1; aug_1 = aug_5; aug_5 = swapt;	\
    swapt = aug_2; aug_2 = aug_6; aug_6 = swapt;	\
    swapt = aug_3; aug_3 = aug_7; aug_7 = swapt;	\
  } else if (aug_8 != (flt_t)0.0) {			\
    flt_t swapt;					\
    swapt = aug_0; aug_0 = aug_8; aug_8 = swapt;	\
    swapt = aug_1; aug_1 = aug_9; aug_9 = swapt;	\
    swapt = aug_2; aug_2 = aug_10; aug_10 = swapt;	\
    swapt = aug_3; aug_3 = aug_11; aug_11 = swapt;	\
  } else						\
    error = 1;						\
							\
  t = aug_4 / aug_0;					\
  aug_5 -= t * aug_1;					\
  aug_6 -= t * aug_2;					\
  aug_7 -= t * aug_3;					\
  t = aug_8 / aug_0;					\
  aug_9 -= t * aug_1;					\
  aug_10 -= t * aug_2;					\
  aug_11 -= t * aug_3;					\
							\
  if (fabs(aug_9) > fabs(aug_5)) {			\
    flt_t swapt;					\
    swapt = aug_4; aug_4 = aug_8; aug_8 = swapt;	\
    swapt = aug_5; aug_5 = aug_9; aug_9 = swapt;	\
    swapt = aug_6; aug_6 = aug_10; aug_10 = swapt;	\
    swapt = aug_7; aug_7 = aug_11; aug_11 = swapt;	\
  }							\
							\
  if (aug_5 != (flt_t)0.0) {				\
  } else if (aug_9 != (flt_t)0.0) {			\
    flt_t swapt;					\
    swapt = aug_4; aug_4 = aug_8; aug_8 = swapt;        \
    swapt = aug_5; aug_5 = aug_9; aug_9 = swapt;	\
    swapt = aug_6; aug_6 = aug_10; aug_10 = swapt;	\
    swapt = aug_7; aug_7 = aug_11; aug_11 = swapt;	\
  }							\
							\
  t = aug_9 / aug_5;					\
  aug_10 -= t * aug_6;					\
  aug_11 -= t * aug_7;					\
							\
  if (aug_10 == (flt_t)0.0)				\
    error = 1;						\
							\
  ans##_2 = aug_11/aug_10;				\
  t = (flt_t)0.0;					\
  t += aug_6 * ans##_2;					\
  ans##_1 = (aug_7-t) / aug_5;				\
  t = (flt_t)0.0;					\
  t += aug_1 * ans##_1;					\
  t += aug_2 * ans##_2;					\
  ans##_0 = (aug_3 - t) / aug_0;			\
}

#endif
