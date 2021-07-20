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

#include <xmmintrin.h>

void mmul3_5_8_3x8v4(const float * restrict A,
                     const float * restrict B,
                           float * restrict C) {
  __m128
    Creg00,Creg04,
    Creg10,Creg14,
    Creg20,Creg24;
  __m128 Areg0,Areg1,Areg2;
  __m128 Breg0,Breg4;
  __m128 Atmp,Btmp;


    /*  Computing C(0:2,0:5)  */

      Areg0 = _mm_load_ps(&A[0]) ;
      Areg1 = _mm_load_ps(&A[8]) ;
      Areg2 = _mm_load_ps(&A[16]) ;

      Breg0 = _mm_load_ps(&B[0]) ;
      Breg4 = _mm_load_ps(&B[4]) ;

      Atmp = (__m128) _mm_shuffle_epi32((__m128i) Areg0,_MM_SHUFFLE(1,1,1,1)) ;
      Creg00 = Breg0 ;
      Creg00 = _mm_mul_ps(Creg00,Atmp) ;
      Creg04 = Breg4 ;
      Creg04 = _mm_mul_ps(Creg04,Atmp) ;

      Atmp = (__m128) _mm_shuffle_epi32((__m128i) Areg1,_MM_SHUFFLE(1,1,1,1)) ;
      Creg10 = Breg0 ;
      Creg10 = _mm_mul_ps(Creg10,Atmp) ;
      Creg14 = Breg4 ;
      Creg14 = _mm_mul_ps(Creg14,Atmp) ;

      Atmp = (__m128) _mm_shuffle_epi32((__m128i) Areg2,_MM_SHUFFLE(1,1,1,1)) ;
      Creg20 = Breg0 ;
      Creg20 = _mm_mul_ps(Creg20,Atmp) ;
      Creg24 = Breg4 ;
      Creg24 = _mm_mul_ps(Creg24,Atmp) ;


      Breg0 = _mm_load_ps(&B[8]) ;
      Breg4 = _mm_load_ps(&B[12]) ;

      Atmp = (__m128) _mm_shuffle_epi32((__m128i) Areg0,_MM_SHUFFLE(2,2,2,2)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg00 = _mm_add_ps(Creg00,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg04 = _mm_add_ps(Creg04,Btmp) ;

      Atmp = (__m128) _mm_shuffle_epi32((__m128i) Areg1,_MM_SHUFFLE(2,2,2,2)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg10 = _mm_add_ps(Creg10,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg14 = _mm_add_ps(Creg14,Btmp) ;

      Atmp = (__m128) _mm_shuffle_epi32((__m128i) Areg2,_MM_SHUFFLE(2,2,2,2)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg20 = _mm_add_ps(Creg20,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg24 = _mm_add_ps(Creg24,Btmp) ;


      Breg0 = _mm_load_ps(&B[16]) ;
      Breg4 = _mm_load_ps(&B[20]) ;

      Atmp = (__m128) _mm_shuffle_epi32((__m128i) Areg0,_MM_SHUFFLE(3,3,3,3)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg00 = _mm_add_ps(Creg00,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg04 = _mm_add_ps(Creg04,Btmp) ;

      Atmp = (__m128) _mm_shuffle_epi32((__m128i) Areg1,_MM_SHUFFLE(3,3,3,3)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg10 = _mm_add_ps(Creg10,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg14 = _mm_add_ps(Creg14,Btmp) ;

      Atmp = (__m128) _mm_shuffle_epi32((__m128i) Areg2,_MM_SHUFFLE(3,3,3,3)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg20 = _mm_add_ps(Creg20,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg24 = _mm_add_ps(Creg24,Btmp) ;


      Areg0 = _mm_load_ps(&A[4]) ;
      Areg1 = _mm_load_ps(&A[12]) ;
      Areg2 = _mm_load_ps(&A[20]) ;

      Breg0 = _mm_load_ps(&B[24]) ;
      Breg4 = _mm_load_ps(&B[28]) ;

      Atmp = (__m128) _mm_shuffle_epi32((__m128i) Areg0,_MM_SHUFFLE(0,0,0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg00 = _mm_add_ps(Creg00,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg04 = _mm_add_ps(Creg04,Btmp) ;

      Atmp = (__m128) _mm_shuffle_epi32((__m128i) Areg1,_MM_SHUFFLE(0,0,0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg10 = _mm_add_ps(Creg10,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg14 = _mm_add_ps(Creg14,Btmp) ;

      Atmp = (__m128) _mm_shuffle_epi32((__m128i) Areg2,_MM_SHUFFLE(0,0,0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg20 = _mm_add_ps(Creg20,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg24 = _mm_add_ps(Creg24,Btmp) ;


      Breg0 = _mm_load_ps(&B[32]) ;
      Breg4 = _mm_load_ps(&B[36]) ;

      Atmp = (__m128) _mm_shuffle_epi32((__m128i) Areg0,_MM_SHUFFLE(1,1,1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg00 = _mm_add_ps(Creg00,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg04 = _mm_add_ps(Creg04,Btmp) ;

      Atmp = (__m128) _mm_shuffle_epi32((__m128i) Areg1,_MM_SHUFFLE(1,1,1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg10 = _mm_add_ps(Creg10,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg14 = _mm_add_ps(Creg14,Btmp) ;

      Atmp = (__m128) _mm_shuffle_epi32((__m128i) Areg2,_MM_SHUFFLE(1,1,1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg20 = _mm_add_ps(Creg20,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg24 = _mm_add_ps(Creg24,Btmp) ;


      _mm_store_ps(&C[0],Creg00) ;
      _mm_store_ps(&C[4],Creg04) ;
      _mm_store_ps(&C[8],Creg10) ;
      _mm_store_ps(&C[12],Creg14) ;
      _mm_store_ps(&C[16],Creg20) ;
      _mm_store_ps(&C[20],Creg24) ;


    /*  Computing C(3:4,0:5)  */

      Areg0 = _mm_load_ps(&A[24]) ;
      Areg1 = _mm_load_ps(&A[32]) ;

      Breg0 = _mm_load_ps(&B[0]) ;
      Breg4 = _mm_load_ps(&B[4]) ;

      Atmp = (__m128) _mm_shuffle_epi32((__m128i) Areg0,_MM_SHUFFLE(1,1,1,1)) ;
      Creg00 = Breg0 ;
      Creg00 = _mm_mul_ps(Creg00,Atmp) ;
      Creg04 = Breg4 ;
      Creg04 = _mm_mul_ps(Creg04,Atmp) ;

      Atmp = (__m128) _mm_shuffle_epi32((__m128i) Areg1,_MM_SHUFFLE(1,1,1,1)) ;
      Creg10 = Breg0 ;
      Creg10 = _mm_mul_ps(Creg10,Atmp) ;
      Creg14 = Breg4 ;
      Creg14 = _mm_mul_ps(Creg14,Atmp) ;


      Breg0 = _mm_load_ps(&B[8]) ;
      Breg4 = _mm_load_ps(&B[12]) ;

      Atmp = (__m128) _mm_shuffle_epi32((__m128i) Areg0,_MM_SHUFFLE(2,2,2,2)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg00 = _mm_add_ps(Creg00,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg04 = _mm_add_ps(Creg04,Btmp) ;

      Atmp = (__m128) _mm_shuffle_epi32((__m128i) Areg1,_MM_SHUFFLE(2,2,2,2)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg10 = _mm_add_ps(Creg10,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg14 = _mm_add_ps(Creg14,Btmp) ;


      Breg0 = _mm_load_ps(&B[16]) ;
      Breg4 = _mm_load_ps(&B[20]) ;

      Atmp = (__m128) _mm_shuffle_epi32((__m128i) Areg0,_MM_SHUFFLE(3,3,3,3)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg00 = _mm_add_ps(Creg00,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg04 = _mm_add_ps(Creg04,Btmp) ;

      Atmp = (__m128) _mm_shuffle_epi32((__m128i) Areg1,_MM_SHUFFLE(3,3,3,3)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg10 = _mm_add_ps(Creg10,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg14 = _mm_add_ps(Creg14,Btmp) ;


      Areg0 = _mm_load_ps(&A[28]) ;
      Areg1 = _mm_load_ps(&A[36]) ;

      Breg0 = _mm_load_ps(&B[24]) ;
      Breg4 = _mm_load_ps(&B[28]) ;

      Atmp = (__m128) _mm_shuffle_epi32((__m128i) Areg0,_MM_SHUFFLE(0,0,0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg00 = _mm_add_ps(Creg00,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg04 = _mm_add_ps(Creg04,Btmp) ;

      Atmp = (__m128) _mm_shuffle_epi32((__m128i) Areg1,_MM_SHUFFLE(0,0,0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg10 = _mm_add_ps(Creg10,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg14 = _mm_add_ps(Creg14,Btmp) ;


      Breg0 = _mm_load_ps(&B[32]) ;
      Breg4 = _mm_load_ps(&B[36]) ;

      Atmp = (__m128) _mm_shuffle_epi32((__m128i) Areg0,_MM_SHUFFLE(1,1,1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg00 = _mm_add_ps(Creg00,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg04 = _mm_add_ps(Creg04,Btmp) ;

      Atmp = (__m128) _mm_shuffle_epi32((__m128i) Areg1,_MM_SHUFFLE(1,1,1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg10 = _mm_add_ps(Creg10,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_ps(Btmp,Atmp) ;
      Creg14 = _mm_add_ps(Creg14,Btmp) ;


      _mm_store_ps(&C[24],Creg00) ;
      _mm_store_ps(&C[28],Creg04) ;
      _mm_store_ps(&C[32],Creg10) ;
      _mm_store_ps(&C[36],Creg14) ;


}
