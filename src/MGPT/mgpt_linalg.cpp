// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   This file is part of the MGPT implementation. See further comments
   in pair_mgpt.cpp and pair_mgpt.h.
------------------------------------------------------------------------- */

#include "mgpt_linalg.h"

#include <cstdio>
#include <cstdlib>

#ifdef x86_SIMD
/* Double precision SSE2 matrix multiplication */
#include "mgpt_mmul3d_526.c.h"
#include "mgpt_mmul3d_744.c.h"

/* Single precision SSE2 matrix multiplication */
#include "mgpt_mmul3_538.c.h"
#include "mgpt_mmul3_748.c.h"

/* Double precision SSE3 product trace */
#define real double
#include "mgpt_ttr_5022.c.h"
#include "mgpt_ttr_7022.c.h"
#undef real

/* Single precision SSE3 product trace */
#define real float
#include "mgpt_ttr_5042.c.h"
#include "mgpt_ttr_7042.c.h"
#undef real

#endif

static void transprod_generic(const double * restrict A,
                              const double * restrict B,
                              double * restrict C) {
  const int lda = 8,n = mgpt_linalg::matrix_size;
  int i,j,k;
  double s;
  for (i = 0; i<n; i++)
    for (j = 0; j<n; j++) {
      s = 0.0;
      for (k = 1; k<=n; k++)
        s = s + A[i*lda+k]*B[j*lda+k];
      C[i*lda+(j+1)] = s;
    }
}

static void transtrace3_generic(const double * restrict A,
                                const double * restrict B0,double * restrict tout0,
                                const double * restrict B1,double * restrict tout1,
                                const double * restrict B2,double * restrict tout2) {
  const int lda = 8,n = mgpt_linalg::matrix_size;
  double t0 = 0.0,t1 = 0.0,t2 = 0.0;
  int i,j;

  for (i = 0; i<n; i++)
    for (j = 1; j<=n; j++) {
      int idx = i*lda + j;
      double atmp = A[idx];
      t0 = t0 + atmp*B0[idx];
      t1 = t1 + atmp*B1[idx];
      t2 = t2 + atmp*B2[idx];
    }
  *tout0 = t0;
  *tout1 = t1;
  *tout2 = t2;
}

static void transprod_error(const double * restrict /*A*/,
                            const double * restrict /*B*/,
                            double * restrict /*C*/) {
  printf("Linear algebra subroutines not initialized (transprod).\n");
  exit(1);
}
static void transtrace3_error(const double * restrict /*A*/,
                              const double * restrict /*B0*/,double * restrict /*tout0*/,
                              const double * restrict /*B1*/,double * restrict /*tout1*/,
                              const double * restrict /*B2*/,double * restrict /*tout2*/) {
  printf("Linear algebra subroutines not initialized (transtrace3).\n");
  exit(1);
}


int mgpt_linalg::matrix_size;

mgpt_linalg::mgpt_linalg() {
  mgpt_linalg::matrix_size = 0;

  tr_mul = transprod_error;
  tr_trace = transtrace3_error;
  single = 0;
}

mgpt_linalg::mgpt_linalg(int n,int single_precision) {

  mgpt_linalg::matrix_size = n;

  tr_mul = transprod_generic;
  tr_trace = transtrace3_generic;
  single = 0;
  msg = "@@@ Choosing generic (unoptimized) linear algebra routines.\n";

#if defined(x86_SIMD)
  if (single_precision) {
    msg = "@@@ Choosing Intel/AMD single precision linear algebra routines.\n";
    if (n == 5) {
      tr_mul = (trmul_fun) mmul3_5_8_3x8v4;
      tr_trace = (trtrace3_fun) ttr_5_8_3_v4r2;
      single = 1;
    } else if (n == 7) {
      tr_mul = (trmul_fun) mmul3_7_8_4x8v4;
      tr_trace = (trtrace3_fun) ttr_7_8_3_v4r2;
      single = 1;
    }
  } else {
    msg = "@@@ Choosing Intel/AMD double precision linear algebra routines.\n";
    if (n == 5) {
      tr_mul = mmul3_5_8_2x6v2;
      tr_trace = ttr_5_8_3_v2r2;
    } else if (n == 7) {
      tr_mul = mmul3_7_8_4x4v2;
      tr_trace = ttr_7_8_3_v2r2;
    }
  }
#endif
}
