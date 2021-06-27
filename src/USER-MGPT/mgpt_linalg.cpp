// clang-format off
/* ----------------------------------------------------------------------
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
   This file is part of the MGPT implementation. See further comments
   in pair_mgpt.cpp and pair_mgpt.h.
------------------------------------------------------------------------- */

#include "mgpt_linalg.h"

#include <cstdio>
#include <cstdlib>

#define restrict __restrict__

#ifdef IBM_BG_SIMD
#include <builtins.h>

/* Double precision 440d (double hummer) matrix multiplication */
#define const
#include "mgpt_mmul_bg_552.c.h"
#include "mgpt_mmul_bg_722.c.h"
#include "mgpt_bgmul_7.c.h"

/* Double precision 440d (double hummer) product trace */
#define real double
#include "mgpt_ttr_5123.c.h"
#include "mgpt_ttr_7123.c.h"
#undef real
#undef const
#endif


#ifdef IBM_BGQ_SIMD
/* Double precision QPX matrix multiplication */
#include "mgpt_mmul_bgq_n5_lda8_2x8.c.h"
#include "mgpt_mmul_bgq_n7_lda8_4x8.c.h"

/* Double precision QPX product trace */
#include "mgpt_ttr_5141.c.h"
#include "mgpt_ttr_7141.c.h"
#endif


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

#if defined(IBM_BG_SIMD) || defined(IBM_BGQ_SIMD)
#define const
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


#if defined(IBM_BG_SIMD) || defined(IBM_BGQ_SIMD)
#undef const
#endif

#undef restrict

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

#ifdef IBM_BG_SIMD
  msg = "@@@ Choosing BG/L optimized linear algebra routines.\n";
  if (n == 5) {
    tr_mul = mmul_bg_5_8_5x2v2;
    tr_trace = ttr_bg_5_8_3_v2r3;
  } else if (n == 7) {
    //tr_mul = mmul_bg_7_8_2x2v2;
    tr_mul = (trmul_fun) bgmul_7;
    tr_trace = ttr_bg_7_8_3_v2r3;
  }
#elif defined(IBM_BGQ_SIMD)
  msg = "@@@ Choosing BG/Q optimized linear algebra routines.\n";
  if (1) {
    if (n == 5) {
      tr_mul = mmul_bgq_n5_lda8_2x8;
      tr_trace = ttr_bg_5_8_3_v4r1;
    } else if (n == 7) {
      tr_mul = mmul_bgq_n7_lda8_4x8;
      tr_trace = ttr_bg_7_8_3_v4r1;
    }
  }
#elif defined(x86_SIMD)
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
