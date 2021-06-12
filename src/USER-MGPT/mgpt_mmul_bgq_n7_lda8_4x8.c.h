// clang-format off
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
   This file is part of the MGPT implementation. See further comments
   in pair_mgpt.cpp and pair_mgpt.h.
------------------------------------------------------------------------- */

#define vector vector4double
#define const
#define vec_fma vec_madd

void mmul_bgq_n7_lda8_4x8(const double * restrict A,
                          const double * restrict B,
                                double * restrict C) {
  vector
    Creg00,Creg04,
    Creg10,Creg14,
    Creg20,Creg24,
    Creg30,Creg34;
  vector Areg0,Areg1,Areg2,Areg3;
  vector Breg0,Breg4;
  vector Atmp;


    /*  Computing C(0:3,0:7)  */

      Areg0 = vec_ld(8*0,A) ;
      Areg1 = vec_ld(8*8,A) ;
      Areg2 = vec_ld(8*16,A) ;
      Areg3 = vec_ld(8*24,A) ;

      Breg0 = vec_ld(8*0,B) ;
      Breg4 = vec_ld(8*4,B) ;

      Atmp = vec_splat(Areg0,1) ;
      Creg00 = vec_mul(Atmp,Breg0) ;
      Creg04 = vec_mul(Atmp,Breg4) ;

      Atmp = vec_splat(Areg1,1) ;
      Creg10 = vec_mul(Atmp,Breg0) ;
      Creg14 = vec_mul(Atmp,Breg4) ;

      Atmp = vec_splat(Areg2,1) ;
      Creg20 = vec_mul(Atmp,Breg0) ;
      Creg24 = vec_mul(Atmp,Breg4) ;

      Atmp = vec_splat(Areg3,1) ;
      Creg30 = vec_mul(Atmp,Breg0) ;
      Creg34 = vec_mul(Atmp,Breg4) ;


      Breg0 = vec_ld(8*8,B) ;
      Breg4 = vec_ld(8*12,B) ;

      Atmp = vec_splat(Areg0,2) ;
      Creg00 = vec_fma(Atmp,Breg0,Creg00) ;
      Creg04 = vec_fma(Atmp,Breg4,Creg04) ;

      Atmp = vec_splat(Areg1,2) ;
      Creg10 = vec_fma(Atmp,Breg0,Creg10) ;
      Creg14 = vec_fma(Atmp,Breg4,Creg14) ;

      Atmp = vec_splat(Areg2,2) ;
      Creg20 = vec_fma(Atmp,Breg0,Creg20) ;
      Creg24 = vec_fma(Atmp,Breg4,Creg24) ;

      Atmp = vec_splat(Areg3,2) ;
      Creg30 = vec_fma(Atmp,Breg0,Creg30) ;
      Creg34 = vec_fma(Atmp,Breg4,Creg34) ;


      Breg0 = vec_ld(8*16,B) ;
      Breg4 = vec_ld(8*20,B) ;

      Atmp = vec_splat(Areg0,3) ;
      Creg00 = vec_fma(Atmp,Breg0,Creg00) ;
      Creg04 = vec_fma(Atmp,Breg4,Creg04) ;

      Atmp = vec_splat(Areg1,3) ;
      Creg10 = vec_fma(Atmp,Breg0,Creg10) ;
      Creg14 = vec_fma(Atmp,Breg4,Creg14) ;

      Atmp = vec_splat(Areg2,3) ;
      Creg20 = vec_fma(Atmp,Breg0,Creg20) ;
      Creg24 = vec_fma(Atmp,Breg4,Creg24) ;

      Atmp = vec_splat(Areg3,3) ;
      Creg30 = vec_fma(Atmp,Breg0,Creg30) ;
      Creg34 = vec_fma(Atmp,Breg4,Creg34) ;


      Areg0 = vec_ld(8*4,A) ;
      Areg1 = vec_ld(8*12,A) ;
      Areg2 = vec_ld(8*20,A) ;
      Areg3 = vec_ld(8*28,A) ;

      Breg0 = vec_ld(8*24,B) ;
      Breg4 = vec_ld(8*28,B) ;

      Atmp = vec_splat(Areg0,0) ;
      Creg00 = vec_fma(Atmp,Breg0,Creg00) ;
      Creg04 = vec_fma(Atmp,Breg4,Creg04) ;

      Atmp = vec_splat(Areg1,0) ;
      Creg10 = vec_fma(Atmp,Breg0,Creg10) ;
      Creg14 = vec_fma(Atmp,Breg4,Creg14) ;

      Atmp = vec_splat(Areg2,0) ;
      Creg20 = vec_fma(Atmp,Breg0,Creg20) ;
      Creg24 = vec_fma(Atmp,Breg4,Creg24) ;

      Atmp = vec_splat(Areg3,0) ;
      Creg30 = vec_fma(Atmp,Breg0,Creg30) ;
      Creg34 = vec_fma(Atmp,Breg4,Creg34) ;


      Breg0 = vec_ld(8*32,B) ;
      Breg4 = vec_ld(8*36,B) ;

      Atmp = vec_splat(Areg0,1) ;
      Creg00 = vec_fma(Atmp,Breg0,Creg00) ;
      Creg04 = vec_fma(Atmp,Breg4,Creg04) ;

      Atmp = vec_splat(Areg1,1) ;
      Creg10 = vec_fma(Atmp,Breg0,Creg10) ;
      Creg14 = vec_fma(Atmp,Breg4,Creg14) ;

      Atmp = vec_splat(Areg2,1) ;
      Creg20 = vec_fma(Atmp,Breg0,Creg20) ;
      Creg24 = vec_fma(Atmp,Breg4,Creg24) ;

      Atmp = vec_splat(Areg3,1) ;
      Creg30 = vec_fma(Atmp,Breg0,Creg30) ;
      Creg34 = vec_fma(Atmp,Breg4,Creg34) ;


      Breg0 = vec_ld(8*40,B) ;
      Breg4 = vec_ld(8*44,B) ;

      Atmp = vec_splat(Areg0,2) ;
      Creg00 = vec_fma(Atmp,Breg0,Creg00) ;
      Creg04 = vec_fma(Atmp,Breg4,Creg04) ;

      Atmp = vec_splat(Areg1,2) ;
      Creg10 = vec_fma(Atmp,Breg0,Creg10) ;
      Creg14 = vec_fma(Atmp,Breg4,Creg14) ;

      Atmp = vec_splat(Areg2,2) ;
      Creg20 = vec_fma(Atmp,Breg0,Creg20) ;
      Creg24 = vec_fma(Atmp,Breg4,Creg24) ;

      Atmp = vec_splat(Areg3,2) ;
      Creg30 = vec_fma(Atmp,Breg0,Creg30) ;
      Creg34 = vec_fma(Atmp,Breg4,Creg34) ;


      Breg0 = vec_ld(8*48,B) ;
      Breg4 = vec_ld(8*52,B) ;

      Atmp = vec_splat(Areg0,3) ;
      Creg00 = vec_fma(Atmp,Breg0,Creg00) ;
      Creg04 = vec_fma(Atmp,Breg4,Creg04) ;

      Atmp = vec_splat(Areg1,3) ;
      Creg10 = vec_fma(Atmp,Breg0,Creg10) ;
      Creg14 = vec_fma(Atmp,Breg4,Creg14) ;

      Atmp = vec_splat(Areg2,3) ;
      Creg20 = vec_fma(Atmp,Breg0,Creg20) ;
      Creg24 = vec_fma(Atmp,Breg4,Creg24) ;

      Atmp = vec_splat(Areg3,3) ;
      Creg30 = vec_fma(Atmp,Breg0,Creg30) ;
      Creg34 = vec_fma(Atmp,Breg4,Creg34) ;


      vec_st(Creg00,8*0,C) ;
      vec_st(Creg04,8*4,C) ;
      vec_st(Creg10,8*8,C) ;
      vec_st(Creg14,8*12,C) ;
      vec_st(Creg20,8*16,C) ;
      vec_st(Creg24,8*20,C) ;
      vec_st(Creg30,8*24,C) ;
      vec_st(Creg34,8*28,C) ;


    /*  Computing C(4:6,0:7)  */

      Areg0 = vec_ld(8*32,A) ;
      Areg1 = vec_ld(8*40,A) ;
      Areg2 = vec_ld(8*48,A) ;

      Breg0 = vec_ld(8*0,B) ;
      Breg4 = vec_ld(8*4,B) ;

      Atmp = vec_splat(Areg0,1) ;
      Creg00 = vec_mul(Atmp,Breg0) ;
      Creg04 = vec_mul(Atmp,Breg4) ;

      Atmp = vec_splat(Areg1,1) ;
      Creg10 = vec_mul(Atmp,Breg0) ;
      Creg14 = vec_mul(Atmp,Breg4) ;

      Atmp = vec_splat(Areg2,1) ;
      Creg20 = vec_mul(Atmp,Breg0) ;
      Creg24 = vec_mul(Atmp,Breg4) ;


      Breg0 = vec_ld(8*8,B) ;
      Breg4 = vec_ld(8*12,B) ;

      Atmp = vec_splat(Areg0,2) ;
      Creg00 = vec_fma(Atmp,Breg0,Creg00) ;
      Creg04 = vec_fma(Atmp,Breg4,Creg04) ;

      Atmp = vec_splat(Areg1,2) ;
      Creg10 = vec_fma(Atmp,Breg0,Creg10) ;
      Creg14 = vec_fma(Atmp,Breg4,Creg14) ;

      Atmp = vec_splat(Areg2,2) ;
      Creg20 = vec_fma(Atmp,Breg0,Creg20) ;
      Creg24 = vec_fma(Atmp,Breg4,Creg24) ;


      Breg0 = vec_ld(8*16,B) ;
      Breg4 = vec_ld(8*20,B) ;

      Atmp = vec_splat(Areg0,3) ;
      Creg00 = vec_fma(Atmp,Breg0,Creg00) ;
      Creg04 = vec_fma(Atmp,Breg4,Creg04) ;

      Atmp = vec_splat(Areg1,3) ;
      Creg10 = vec_fma(Atmp,Breg0,Creg10) ;
      Creg14 = vec_fma(Atmp,Breg4,Creg14) ;

      Atmp = vec_splat(Areg2,3) ;
      Creg20 = vec_fma(Atmp,Breg0,Creg20) ;
      Creg24 = vec_fma(Atmp,Breg4,Creg24) ;


      Areg0 = vec_ld(8*36,A) ;
      Areg1 = vec_ld(8*44,A) ;
      Areg2 = vec_ld(8*52,A) ;

      Breg0 = vec_ld(8*24,B) ;
      Breg4 = vec_ld(8*28,B) ;

      Atmp = vec_splat(Areg0,0) ;
      Creg00 = vec_fma(Atmp,Breg0,Creg00) ;
      Creg04 = vec_fma(Atmp,Breg4,Creg04) ;

      Atmp = vec_splat(Areg1,0) ;
      Creg10 = vec_fma(Atmp,Breg0,Creg10) ;
      Creg14 = vec_fma(Atmp,Breg4,Creg14) ;

      Atmp = vec_splat(Areg2,0) ;
      Creg20 = vec_fma(Atmp,Breg0,Creg20) ;
      Creg24 = vec_fma(Atmp,Breg4,Creg24) ;


      Breg0 = vec_ld(8*32,B) ;
      Breg4 = vec_ld(8*36,B) ;

      Atmp = vec_splat(Areg0,1) ;
      Creg00 = vec_fma(Atmp,Breg0,Creg00) ;
      Creg04 = vec_fma(Atmp,Breg4,Creg04) ;

      Atmp = vec_splat(Areg1,1) ;
      Creg10 = vec_fma(Atmp,Breg0,Creg10) ;
      Creg14 = vec_fma(Atmp,Breg4,Creg14) ;

      Atmp = vec_splat(Areg2,1) ;
      Creg20 = vec_fma(Atmp,Breg0,Creg20) ;
      Creg24 = vec_fma(Atmp,Breg4,Creg24) ;


      Breg0 = vec_ld(8*40,B) ;
      Breg4 = vec_ld(8*44,B) ;

      Atmp = vec_splat(Areg0,2) ;
      Creg00 = vec_fma(Atmp,Breg0,Creg00) ;
      Creg04 = vec_fma(Atmp,Breg4,Creg04) ;

      Atmp = vec_splat(Areg1,2) ;
      Creg10 = vec_fma(Atmp,Breg0,Creg10) ;
      Creg14 = vec_fma(Atmp,Breg4,Creg14) ;

      Atmp = vec_splat(Areg2,2) ;
      Creg20 = vec_fma(Atmp,Breg0,Creg20) ;
      Creg24 = vec_fma(Atmp,Breg4,Creg24) ;


      Breg0 = vec_ld(8*48,B) ;
      Breg4 = vec_ld(8*52,B) ;

      Atmp = vec_splat(Areg0,3) ;
      Creg00 = vec_fma(Atmp,Breg0,Creg00) ;
      Creg04 = vec_fma(Atmp,Breg4,Creg04) ;

      Atmp = vec_splat(Areg1,3) ;
      Creg10 = vec_fma(Atmp,Breg0,Creg10) ;
      Creg14 = vec_fma(Atmp,Breg4,Creg14) ;

      Atmp = vec_splat(Areg2,3) ;
      Creg20 = vec_fma(Atmp,Breg0,Creg20) ;
      Creg24 = vec_fma(Atmp,Breg4,Creg24) ;


      vec_st(Creg00,8*32,C) ;
      vec_st(Creg04,8*36,C) ;
      vec_st(Creg10,8*40,C) ;
      vec_st(Creg14,8*44,C) ;
      vec_st(Creg20,8*48,C) ;
      vec_st(Creg24,8*52,C) ;


}

#undef vector
#undef const
#undef vec_fma
