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

void mmul3_5_8_2x6v2(const double * restrict A,
                     const double * restrict B,
                           double * restrict C) {
  __m128d
    Creg00,Creg02,Creg04,
    Creg10,Creg12,Creg14;
  __m128d Areg0,Areg1;
  __m128d Breg0,Breg2,Breg4;
  __m128d Atmp,Btmp;


    /*  Computing C(0:1,0:5)  */

      Areg0 = _mm_load_pd(&A[0]) ;
      Areg1 = _mm_load_pd(&A[8]) ;

      Breg0 = _mm_load_pd(&B[0]) ;
      Breg2 = _mm_load_pd(&B[2]) ;
      Breg4 = _mm_load_pd(&B[4]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(1,1)) ;
      Creg00 = Breg0 ;
      Creg00 = _mm_mul_pd(Creg00,Atmp) ;
      Creg02 = Breg2 ;
      Creg02 = _mm_mul_pd(Creg02,Atmp) ;
      Creg04 = Breg4 ;
      Creg04 = _mm_mul_pd(Creg04,Atmp) ;

      Atmp = Areg1 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg1,_MM_SHUFFLE2(1,1)) ;
      Creg10 = Breg0 ;
      Creg10 = _mm_mul_pd(Creg10,Atmp) ;
      Creg12 = Breg2 ;
      Creg12 = _mm_mul_pd(Creg12,Atmp) ;
      Creg14 = Breg4 ;
      Creg14 = _mm_mul_pd(Creg14,Atmp) ;


      Areg0 = _mm_load_pd(&A[2]) ;
      Areg1 = _mm_load_pd(&A[10]) ;

      Breg0 = _mm_load_pd(&B[8]) ;
      Breg2 = _mm_load_pd(&B[10]) ;
      Breg4 = _mm_load_pd(&B[12]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg00 = _mm_add_pd(Creg00,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg02 = _mm_add_pd(Creg02,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg04 = _mm_add_pd(Creg04,Btmp) ;

      Atmp = Areg1 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg1,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg10 = _mm_add_pd(Creg10,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg12 = _mm_add_pd(Creg12,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg14 = _mm_add_pd(Creg14,Btmp) ;


      Breg0 = _mm_load_pd(&B[16]) ;
      Breg2 = _mm_load_pd(&B[18]) ;
      Breg4 = _mm_load_pd(&B[20]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg00 = _mm_add_pd(Creg00,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg02 = _mm_add_pd(Creg02,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg04 = _mm_add_pd(Creg04,Btmp) ;

      Atmp = Areg1 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg1,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg10 = _mm_add_pd(Creg10,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg12 = _mm_add_pd(Creg12,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg14 = _mm_add_pd(Creg14,Btmp) ;


      Areg0 = _mm_load_pd(&A[4]) ;
      Areg1 = _mm_load_pd(&A[12]) ;

      Breg0 = _mm_load_pd(&B[24]) ;
      Breg2 = _mm_load_pd(&B[26]) ;
      Breg4 = _mm_load_pd(&B[28]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg00 = _mm_add_pd(Creg00,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg02 = _mm_add_pd(Creg02,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg04 = _mm_add_pd(Creg04,Btmp) ;

      Atmp = Areg1 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg1,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg10 = _mm_add_pd(Creg10,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg12 = _mm_add_pd(Creg12,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg14 = _mm_add_pd(Creg14,Btmp) ;


      Breg0 = _mm_load_pd(&B[32]) ;
      Breg2 = _mm_load_pd(&B[34]) ;
      Breg4 = _mm_load_pd(&B[36]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg00 = _mm_add_pd(Creg00,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg02 = _mm_add_pd(Creg02,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg04 = _mm_add_pd(Creg04,Btmp) ;

      Atmp = Areg1 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg1,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg10 = _mm_add_pd(Creg10,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg12 = _mm_add_pd(Creg12,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg14 = _mm_add_pd(Creg14,Btmp) ;


      _mm_store_pd(&C[0],Creg00) ;
      _mm_store_pd(&C[2],Creg02) ;
      _mm_store_pd(&C[4],Creg04) ;
      _mm_store_pd(&C[8],Creg10) ;
      _mm_store_pd(&C[10],Creg12) ;
      _mm_store_pd(&C[12],Creg14) ;


    /*  Computing C(2:3,0:5)  */

      Areg0 = _mm_load_pd(&A[16]) ;
      Areg1 = _mm_load_pd(&A[24]) ;

      Breg0 = _mm_load_pd(&B[0]) ;
      Breg2 = _mm_load_pd(&B[2]) ;
      Breg4 = _mm_load_pd(&B[4]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(1,1)) ;
      Creg00 = Breg0 ;
      Creg00 = _mm_mul_pd(Creg00,Atmp) ;
      Creg02 = Breg2 ;
      Creg02 = _mm_mul_pd(Creg02,Atmp) ;
      Creg04 = Breg4 ;
      Creg04 = _mm_mul_pd(Creg04,Atmp) ;

      Atmp = Areg1 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg1,_MM_SHUFFLE2(1,1)) ;
      Creg10 = Breg0 ;
      Creg10 = _mm_mul_pd(Creg10,Atmp) ;
      Creg12 = Breg2 ;
      Creg12 = _mm_mul_pd(Creg12,Atmp) ;
      Creg14 = Breg4 ;
      Creg14 = _mm_mul_pd(Creg14,Atmp) ;


      Areg0 = _mm_load_pd(&A[18]) ;
      Areg1 = _mm_load_pd(&A[26]) ;

      Breg0 = _mm_load_pd(&B[8]) ;
      Breg2 = _mm_load_pd(&B[10]) ;
      Breg4 = _mm_load_pd(&B[12]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg00 = _mm_add_pd(Creg00,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg02 = _mm_add_pd(Creg02,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg04 = _mm_add_pd(Creg04,Btmp) ;

      Atmp = Areg1 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg1,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg10 = _mm_add_pd(Creg10,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg12 = _mm_add_pd(Creg12,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg14 = _mm_add_pd(Creg14,Btmp) ;


      Breg0 = _mm_load_pd(&B[16]) ;
      Breg2 = _mm_load_pd(&B[18]) ;
      Breg4 = _mm_load_pd(&B[20]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg00 = _mm_add_pd(Creg00,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg02 = _mm_add_pd(Creg02,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg04 = _mm_add_pd(Creg04,Btmp) ;

      Atmp = Areg1 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg1,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg10 = _mm_add_pd(Creg10,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg12 = _mm_add_pd(Creg12,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg14 = _mm_add_pd(Creg14,Btmp) ;


      Areg0 = _mm_load_pd(&A[20]) ;
      Areg1 = _mm_load_pd(&A[28]) ;

      Breg0 = _mm_load_pd(&B[24]) ;
      Breg2 = _mm_load_pd(&B[26]) ;
      Breg4 = _mm_load_pd(&B[28]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg00 = _mm_add_pd(Creg00,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg02 = _mm_add_pd(Creg02,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg04 = _mm_add_pd(Creg04,Btmp) ;

      Atmp = Areg1 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg1,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg10 = _mm_add_pd(Creg10,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg12 = _mm_add_pd(Creg12,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg14 = _mm_add_pd(Creg14,Btmp) ;


      Breg0 = _mm_load_pd(&B[32]) ;
      Breg2 = _mm_load_pd(&B[34]) ;
      Breg4 = _mm_load_pd(&B[36]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg00 = _mm_add_pd(Creg00,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg02 = _mm_add_pd(Creg02,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg04 = _mm_add_pd(Creg04,Btmp) ;

      Atmp = Areg1 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg1,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg10 = _mm_add_pd(Creg10,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg12 = _mm_add_pd(Creg12,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg14 = _mm_add_pd(Creg14,Btmp) ;


      _mm_store_pd(&C[16],Creg00) ;
      _mm_store_pd(&C[18],Creg02) ;
      _mm_store_pd(&C[20],Creg04) ;
      _mm_store_pd(&C[24],Creg10) ;
      _mm_store_pd(&C[26],Creg12) ;
      _mm_store_pd(&C[28],Creg14) ;


    /*  Computing C(4:4,0:5)  */

      Areg0 = _mm_load_pd(&A[32]) ;

      Breg0 = _mm_load_pd(&B[0]) ;
      Breg2 = _mm_load_pd(&B[2]) ;
      Breg4 = _mm_load_pd(&B[4]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(1,1)) ;
      Creg00 = Breg0 ;
      Creg00 = _mm_mul_pd(Creg00,Atmp) ;
      Creg02 = Breg2 ;
      Creg02 = _mm_mul_pd(Creg02,Atmp) ;
      Creg04 = Breg4 ;
      Creg04 = _mm_mul_pd(Creg04,Atmp) ;


      Areg0 = _mm_load_pd(&A[34]) ;

      Breg0 = _mm_load_pd(&B[8]) ;
      Breg2 = _mm_load_pd(&B[10]) ;
      Breg4 = _mm_load_pd(&B[12]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg00 = _mm_add_pd(Creg00,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg02 = _mm_add_pd(Creg02,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg04 = _mm_add_pd(Creg04,Btmp) ;


      Breg0 = _mm_load_pd(&B[16]) ;
      Breg2 = _mm_load_pd(&B[18]) ;
      Breg4 = _mm_load_pd(&B[20]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg00 = _mm_add_pd(Creg00,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg02 = _mm_add_pd(Creg02,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg04 = _mm_add_pd(Creg04,Btmp) ;


      Areg0 = _mm_load_pd(&A[36]) ;

      Breg0 = _mm_load_pd(&B[24]) ;
      Breg2 = _mm_load_pd(&B[26]) ;
      Breg4 = _mm_load_pd(&B[28]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg00 = _mm_add_pd(Creg00,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg02 = _mm_add_pd(Creg02,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg04 = _mm_add_pd(Creg04,Btmp) ;


      Breg0 = _mm_load_pd(&B[32]) ;
      Breg2 = _mm_load_pd(&B[34]) ;
      Breg4 = _mm_load_pd(&B[36]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg00 = _mm_add_pd(Creg00,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg02 = _mm_add_pd(Creg02,Btmp) ;
      Btmp = Breg4 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg04 = _mm_add_pd(Creg04,Btmp) ;


      _mm_store_pd(&C[32],Creg00) ;
      _mm_store_pd(&C[34],Creg02) ;
      _mm_store_pd(&C[36],Creg04) ;


}
