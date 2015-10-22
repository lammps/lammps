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
   This file is part of the MGPT implementation. See further comments
   in pair_mgpt.cpp and pair_mgpt.h.
------------------------------------------------------------------------- */

#ifndef MGPT_LINALG__
#define MGPT_LINALG__


#ifdef __bg__

  #ifdef __bgq__
    #ifdef __VECTOR4DOUBLE__
      #define IBM_BGQ_SIMD
    #endif
  #else
    #define IBM_BG_SIMD
  #endif

#elif defined(__SSE3__)
  #define x86_SIMD
#endif

#define restrict __restrict__

#if defined(IBM_BG_SIMD) || defined(IBM_BGQ_SIMD)
#define const
#endif
typedef void (*trmul_fun) (const double * restrict A,
			   const double * restrict B,
			   double * restrict C);

typedef void (*trtrace3_fun) (const double * restrict A,
			      const double * restrict B1,double * restrict t1,
			      const double * restrict B2,double * restrict t2,
			      const double * restrict B3,double * restrict t3);
#if defined(IBM_BG_SIMD) || defined(IBM_BGQ_SIMD)
#undef const
#endif

class mgpt_linalg {
 public:
  static int matrix_size;

  trmul_fun tr_mul;
  trtrace3_fun tr_trace;
  int single;
  const char *msg;

  mgpt_linalg();
  mgpt_linalg(int n,int single_precision);
};

#undef restrict

#endif
