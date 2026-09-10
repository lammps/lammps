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

/* ----------------------------------------------------------------------
   This file is part of the MGPT implementation. See further comments
   in pair_mgpt.cpp and pair_mgpt.h.
------------------------------------------------------------------------- */

#ifndef MGPT_LINALG__
#define MGPT_LINALG__


#if defined(__SSE3__)
#define x86_SIMD
#define restrict __restrict__
#else
#define restrict
#endif


typedef void (*trmul_fun) (const double * restrict A,
                           const double * restrict B,
                           double * restrict C);

typedef void (*trtrace3_fun) (const double * restrict A,
                              const double * restrict B1,double * restrict t1,
                              const double * restrict B2,double * restrict t2,
                              const double * restrict B3,double * restrict t3);

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

#endif
