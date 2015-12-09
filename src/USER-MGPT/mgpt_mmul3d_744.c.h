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

#include <xmmintrin.h>

void mmul3_7_8_4x4v2(const double * restrict A,
                     const double * restrict B,
                           double * restrict C) {
  __m128d
    Creg00,Creg02,
    Creg10,Creg12,
    Creg20,Creg22,
    Creg30,Creg32;
  __m128d Areg0,Areg1,Areg2,Areg3;
  __m128d Breg0,Breg2;
  __m128d Atmp,Btmp;


    /*  Computing C(0:3,0:3)  */

      Areg0 = _mm_load_pd(&A[0]) ;
      Areg1 = _mm_load_pd(&A[8]) ;
      Areg2 = _mm_load_pd(&A[16]) ;
      Areg3 = _mm_load_pd(&A[24]) ;

      Breg0 = _mm_load_pd(&B[0]) ;
      Breg2 = _mm_load_pd(&B[2]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(1,1)) ;
      Creg00 = Breg0 ;
      Creg00 = _mm_mul_pd(Creg00,Atmp) ;
      Creg02 = Breg2 ;
      Creg02 = _mm_mul_pd(Creg02,Atmp) ;

      Atmp = Areg1 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg1,_MM_SHUFFLE2(1,1)) ;
      Creg10 = Breg0 ;
      Creg10 = _mm_mul_pd(Creg10,Atmp) ;
      Creg12 = Breg2 ;
      Creg12 = _mm_mul_pd(Creg12,Atmp) ;

      Atmp = Areg2 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg2,_MM_SHUFFLE2(1,1)) ;
      Creg20 = Breg0 ;
      Creg20 = _mm_mul_pd(Creg20,Atmp) ;
      Creg22 = Breg2 ;
      Creg22 = _mm_mul_pd(Creg22,Atmp) ;

      Atmp = Areg3 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg3,_MM_SHUFFLE2(1,1)) ;
      Creg30 = Breg0 ;
      Creg30 = _mm_mul_pd(Creg30,Atmp) ;
      Creg32 = Breg2 ;
      Creg32 = _mm_mul_pd(Creg32,Atmp) ;


      Areg0 = _mm_load_pd(&A[2]) ;
      Areg1 = _mm_load_pd(&A[10]) ;
      Areg2 = _mm_load_pd(&A[18]) ;
      Areg3 = _mm_load_pd(&A[26]) ;

      Breg0 = _mm_load_pd(&B[8]) ;
      Breg2 = _mm_load_pd(&B[10]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg00 = _mm_add_pd(Creg00,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg02 = _mm_add_pd(Creg02,Btmp) ;

      Atmp = Areg1 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg1,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg10 = _mm_add_pd(Creg10,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg12 = _mm_add_pd(Creg12,Btmp) ;

      Atmp = Areg2 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg2,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg20 = _mm_add_pd(Creg20,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg22 = _mm_add_pd(Creg22,Btmp) ;

      Atmp = Areg3 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg3,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg30 = _mm_add_pd(Creg30,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg32 = _mm_add_pd(Creg32,Btmp) ;


      Breg0 = _mm_load_pd(&B[16]) ;
      Breg2 = _mm_load_pd(&B[18]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg00 = _mm_add_pd(Creg00,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg02 = _mm_add_pd(Creg02,Btmp) ;

      Atmp = Areg1 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg1,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg10 = _mm_add_pd(Creg10,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg12 = _mm_add_pd(Creg12,Btmp) ;

      Atmp = Areg2 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg2,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg20 = _mm_add_pd(Creg20,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg22 = _mm_add_pd(Creg22,Btmp) ;

      Atmp = Areg3 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg3,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg30 = _mm_add_pd(Creg30,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg32 = _mm_add_pd(Creg32,Btmp) ;


      Areg0 = _mm_load_pd(&A[4]) ;
      Areg1 = _mm_load_pd(&A[12]) ;
      Areg2 = _mm_load_pd(&A[20]) ;
      Areg3 = _mm_load_pd(&A[28]) ;

      Breg0 = _mm_load_pd(&B[24]) ;
      Breg2 = _mm_load_pd(&B[26]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg00 = _mm_add_pd(Creg00,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg02 = _mm_add_pd(Creg02,Btmp) ;

      Atmp = Areg1 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg1,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg10 = _mm_add_pd(Creg10,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg12 = _mm_add_pd(Creg12,Btmp) ;

      Atmp = Areg2 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg2,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg20 = _mm_add_pd(Creg20,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg22 = _mm_add_pd(Creg22,Btmp) ;

      Atmp = Areg3 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg3,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg30 = _mm_add_pd(Creg30,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg32 = _mm_add_pd(Creg32,Btmp) ;


      Breg0 = _mm_load_pd(&B[32]) ;
      Breg2 = _mm_load_pd(&B[34]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg00 = _mm_add_pd(Creg00,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg02 = _mm_add_pd(Creg02,Btmp) ;

      Atmp = Areg1 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg1,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg10 = _mm_add_pd(Creg10,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg12 = _mm_add_pd(Creg12,Btmp) ;

      Atmp = Areg2 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg2,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg20 = _mm_add_pd(Creg20,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg22 = _mm_add_pd(Creg22,Btmp) ;

      Atmp = Areg3 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg3,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg30 = _mm_add_pd(Creg30,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg32 = _mm_add_pd(Creg32,Btmp) ;


      Areg0 = _mm_load_pd(&A[6]) ;
      Areg1 = _mm_load_pd(&A[14]) ;
      Areg2 = _mm_load_pd(&A[22]) ;
      Areg3 = _mm_load_pd(&A[30]) ;

      Breg0 = _mm_load_pd(&B[40]) ;
      Breg2 = _mm_load_pd(&B[42]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg00 = _mm_add_pd(Creg00,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg02 = _mm_add_pd(Creg02,Btmp) ;

      Atmp = Areg1 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg1,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg10 = _mm_add_pd(Creg10,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg12 = _mm_add_pd(Creg12,Btmp) ;

      Atmp = Areg2 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg2,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg20 = _mm_add_pd(Creg20,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg22 = _mm_add_pd(Creg22,Btmp) ;

      Atmp = Areg3 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg3,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg30 = _mm_add_pd(Creg30,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg32 = _mm_add_pd(Creg32,Btmp) ;


      Breg0 = _mm_load_pd(&B[48]) ;
      Breg2 = _mm_load_pd(&B[50]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg00 = _mm_add_pd(Creg00,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg02 = _mm_add_pd(Creg02,Btmp) ;

      Atmp = Areg1 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg1,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg10 = _mm_add_pd(Creg10,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg12 = _mm_add_pd(Creg12,Btmp) ;

      Atmp = Areg2 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg2,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg20 = _mm_add_pd(Creg20,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg22 = _mm_add_pd(Creg22,Btmp) ;

      Atmp = Areg3 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg3,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg30 = _mm_add_pd(Creg30,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg32 = _mm_add_pd(Creg32,Btmp) ;


      _mm_store_pd(&C[0],Creg00) ;
      _mm_store_pd(&C[2],Creg02) ;
      _mm_store_pd(&C[8],Creg10) ;
      _mm_store_pd(&C[10],Creg12) ;
      _mm_store_pd(&C[16],Creg20) ;
      _mm_store_pd(&C[18],Creg22) ;
      _mm_store_pd(&C[24],Creg30) ;
      _mm_store_pd(&C[26],Creg32) ;


    /*  Computing C(0:3,4:7)  */

      Areg0 = _mm_load_pd(&A[0]) ;
      Areg1 = _mm_load_pd(&A[8]) ;
      Areg2 = _mm_load_pd(&A[16]) ;
      Areg3 = _mm_load_pd(&A[24]) ;

      Breg0 = _mm_load_pd(&B[4]) ;
      Breg2 = _mm_load_pd(&B[6]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(1,1)) ;
      Creg00 = Breg0 ;
      Creg00 = _mm_mul_pd(Creg00,Atmp) ;
      Creg02 = Breg2 ;
      Creg02 = _mm_mul_pd(Creg02,Atmp) ;

      Atmp = Areg1 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg1,_MM_SHUFFLE2(1,1)) ;
      Creg10 = Breg0 ;
      Creg10 = _mm_mul_pd(Creg10,Atmp) ;
      Creg12 = Breg2 ;
      Creg12 = _mm_mul_pd(Creg12,Atmp) ;

      Atmp = Areg2 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg2,_MM_SHUFFLE2(1,1)) ;
      Creg20 = Breg0 ;
      Creg20 = _mm_mul_pd(Creg20,Atmp) ;
      Creg22 = Breg2 ;
      Creg22 = _mm_mul_pd(Creg22,Atmp) ;

      Atmp = Areg3 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg3,_MM_SHUFFLE2(1,1)) ;
      Creg30 = Breg0 ;
      Creg30 = _mm_mul_pd(Creg30,Atmp) ;
      Creg32 = Breg2 ;
      Creg32 = _mm_mul_pd(Creg32,Atmp) ;


      Areg0 = _mm_load_pd(&A[2]) ;
      Areg1 = _mm_load_pd(&A[10]) ;
      Areg2 = _mm_load_pd(&A[18]) ;
      Areg3 = _mm_load_pd(&A[26]) ;

      Breg0 = _mm_load_pd(&B[12]) ;
      Breg2 = _mm_load_pd(&B[14]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg00 = _mm_add_pd(Creg00,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg02 = _mm_add_pd(Creg02,Btmp) ;

      Atmp = Areg1 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg1,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg10 = _mm_add_pd(Creg10,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg12 = _mm_add_pd(Creg12,Btmp) ;

      Atmp = Areg2 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg2,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg20 = _mm_add_pd(Creg20,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg22 = _mm_add_pd(Creg22,Btmp) ;

      Atmp = Areg3 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg3,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg30 = _mm_add_pd(Creg30,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg32 = _mm_add_pd(Creg32,Btmp) ;


      Breg0 = _mm_load_pd(&B[20]) ;
      Breg2 = _mm_load_pd(&B[22]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg00 = _mm_add_pd(Creg00,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg02 = _mm_add_pd(Creg02,Btmp) ;

      Atmp = Areg1 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg1,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg10 = _mm_add_pd(Creg10,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg12 = _mm_add_pd(Creg12,Btmp) ;

      Atmp = Areg2 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg2,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg20 = _mm_add_pd(Creg20,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg22 = _mm_add_pd(Creg22,Btmp) ;

      Atmp = Areg3 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg3,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg30 = _mm_add_pd(Creg30,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg32 = _mm_add_pd(Creg32,Btmp) ;


      Areg0 = _mm_load_pd(&A[4]) ;
      Areg1 = _mm_load_pd(&A[12]) ;
      Areg2 = _mm_load_pd(&A[20]) ;
      Areg3 = _mm_load_pd(&A[28]) ;

      Breg0 = _mm_load_pd(&B[28]) ;
      Breg2 = _mm_load_pd(&B[30]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg00 = _mm_add_pd(Creg00,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg02 = _mm_add_pd(Creg02,Btmp) ;

      Atmp = Areg1 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg1,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg10 = _mm_add_pd(Creg10,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg12 = _mm_add_pd(Creg12,Btmp) ;

      Atmp = Areg2 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg2,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg20 = _mm_add_pd(Creg20,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg22 = _mm_add_pd(Creg22,Btmp) ;

      Atmp = Areg3 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg3,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg30 = _mm_add_pd(Creg30,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg32 = _mm_add_pd(Creg32,Btmp) ;


      Breg0 = _mm_load_pd(&B[36]) ;
      Breg2 = _mm_load_pd(&B[38]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg00 = _mm_add_pd(Creg00,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg02 = _mm_add_pd(Creg02,Btmp) ;

      Atmp = Areg1 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg1,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg10 = _mm_add_pd(Creg10,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg12 = _mm_add_pd(Creg12,Btmp) ;

      Atmp = Areg2 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg2,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg20 = _mm_add_pd(Creg20,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg22 = _mm_add_pd(Creg22,Btmp) ;

      Atmp = Areg3 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg3,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg30 = _mm_add_pd(Creg30,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg32 = _mm_add_pd(Creg32,Btmp) ;


      Areg0 = _mm_load_pd(&A[6]) ;
      Areg1 = _mm_load_pd(&A[14]) ;
      Areg2 = _mm_load_pd(&A[22]) ;
      Areg3 = _mm_load_pd(&A[30]) ;

      Breg0 = _mm_load_pd(&B[44]) ;
      Breg2 = _mm_load_pd(&B[46]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg00 = _mm_add_pd(Creg00,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg02 = _mm_add_pd(Creg02,Btmp) ;

      Atmp = Areg1 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg1,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg10 = _mm_add_pd(Creg10,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg12 = _mm_add_pd(Creg12,Btmp) ;

      Atmp = Areg2 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg2,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg20 = _mm_add_pd(Creg20,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg22 = _mm_add_pd(Creg22,Btmp) ;

      Atmp = Areg3 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg3,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg30 = _mm_add_pd(Creg30,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg32 = _mm_add_pd(Creg32,Btmp) ;


      Breg0 = _mm_load_pd(&B[52]) ;
      Breg2 = _mm_load_pd(&B[54]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg00 = _mm_add_pd(Creg00,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg02 = _mm_add_pd(Creg02,Btmp) ;

      Atmp = Areg1 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg1,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg10 = _mm_add_pd(Creg10,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg12 = _mm_add_pd(Creg12,Btmp) ;

      Atmp = Areg2 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg2,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg20 = _mm_add_pd(Creg20,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg22 = _mm_add_pd(Creg22,Btmp) ;

      Atmp = Areg3 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg3,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg30 = _mm_add_pd(Creg30,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg32 = _mm_add_pd(Creg32,Btmp) ;


      _mm_store_pd(&C[4],Creg00) ;
      _mm_store_pd(&C[6],Creg02) ;
      _mm_store_pd(&C[12],Creg10) ;
      _mm_store_pd(&C[14],Creg12) ;
      _mm_store_pd(&C[20],Creg20) ;
      _mm_store_pd(&C[22],Creg22) ;
      _mm_store_pd(&C[28],Creg30) ;
      _mm_store_pd(&C[30],Creg32) ;


    /*  Computing C(4:6,0:3)  */

      Areg0 = _mm_load_pd(&A[32]) ;
      Areg1 = _mm_load_pd(&A[40]) ;
      Areg2 = _mm_load_pd(&A[48]) ;

      Breg0 = _mm_load_pd(&B[0]) ;
      Breg2 = _mm_load_pd(&B[2]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(1,1)) ;
      Creg00 = Breg0 ;
      Creg00 = _mm_mul_pd(Creg00,Atmp) ;
      Creg02 = Breg2 ;
      Creg02 = _mm_mul_pd(Creg02,Atmp) ;

      Atmp = Areg1 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg1,_MM_SHUFFLE2(1,1)) ;
      Creg10 = Breg0 ;
      Creg10 = _mm_mul_pd(Creg10,Atmp) ;
      Creg12 = Breg2 ;
      Creg12 = _mm_mul_pd(Creg12,Atmp) ;

      Atmp = Areg2 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg2,_MM_SHUFFLE2(1,1)) ;
      Creg20 = Breg0 ;
      Creg20 = _mm_mul_pd(Creg20,Atmp) ;
      Creg22 = Breg2 ;
      Creg22 = _mm_mul_pd(Creg22,Atmp) ;


      Areg0 = _mm_load_pd(&A[34]) ;
      Areg1 = _mm_load_pd(&A[42]) ;
      Areg2 = _mm_load_pd(&A[50]) ;

      Breg0 = _mm_load_pd(&B[8]) ;
      Breg2 = _mm_load_pd(&B[10]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg00 = _mm_add_pd(Creg00,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg02 = _mm_add_pd(Creg02,Btmp) ;

      Atmp = Areg1 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg1,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg10 = _mm_add_pd(Creg10,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg12 = _mm_add_pd(Creg12,Btmp) ;

      Atmp = Areg2 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg2,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg20 = _mm_add_pd(Creg20,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg22 = _mm_add_pd(Creg22,Btmp) ;


      Breg0 = _mm_load_pd(&B[16]) ;
      Breg2 = _mm_load_pd(&B[18]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg00 = _mm_add_pd(Creg00,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg02 = _mm_add_pd(Creg02,Btmp) ;

      Atmp = Areg1 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg1,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg10 = _mm_add_pd(Creg10,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg12 = _mm_add_pd(Creg12,Btmp) ;

      Atmp = Areg2 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg2,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg20 = _mm_add_pd(Creg20,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg22 = _mm_add_pd(Creg22,Btmp) ;


      Areg0 = _mm_load_pd(&A[36]) ;
      Areg1 = _mm_load_pd(&A[44]) ;
      Areg2 = _mm_load_pd(&A[52]) ;

      Breg0 = _mm_load_pd(&B[24]) ;
      Breg2 = _mm_load_pd(&B[26]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg00 = _mm_add_pd(Creg00,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg02 = _mm_add_pd(Creg02,Btmp) ;

      Atmp = Areg1 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg1,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg10 = _mm_add_pd(Creg10,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg12 = _mm_add_pd(Creg12,Btmp) ;

      Atmp = Areg2 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg2,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg20 = _mm_add_pd(Creg20,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg22 = _mm_add_pd(Creg22,Btmp) ;


      Breg0 = _mm_load_pd(&B[32]) ;
      Breg2 = _mm_load_pd(&B[34]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg00 = _mm_add_pd(Creg00,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg02 = _mm_add_pd(Creg02,Btmp) ;

      Atmp = Areg1 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg1,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg10 = _mm_add_pd(Creg10,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg12 = _mm_add_pd(Creg12,Btmp) ;

      Atmp = Areg2 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg2,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg20 = _mm_add_pd(Creg20,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg22 = _mm_add_pd(Creg22,Btmp) ;


      Areg0 = _mm_load_pd(&A[38]) ;
      Areg1 = _mm_load_pd(&A[46]) ;
      Areg2 = _mm_load_pd(&A[54]) ;

      Breg0 = _mm_load_pd(&B[40]) ;
      Breg2 = _mm_load_pd(&B[42]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg00 = _mm_add_pd(Creg00,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg02 = _mm_add_pd(Creg02,Btmp) ;

      Atmp = Areg1 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg1,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg10 = _mm_add_pd(Creg10,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg12 = _mm_add_pd(Creg12,Btmp) ;

      Atmp = Areg2 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg2,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg20 = _mm_add_pd(Creg20,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg22 = _mm_add_pd(Creg22,Btmp) ;


      Breg0 = _mm_load_pd(&B[48]) ;
      Breg2 = _mm_load_pd(&B[50]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg00 = _mm_add_pd(Creg00,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg02 = _mm_add_pd(Creg02,Btmp) ;

      Atmp = Areg1 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg1,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg10 = _mm_add_pd(Creg10,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg12 = _mm_add_pd(Creg12,Btmp) ;

      Atmp = Areg2 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg2,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg20 = _mm_add_pd(Creg20,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg22 = _mm_add_pd(Creg22,Btmp) ;


      _mm_store_pd(&C[32],Creg00) ;
      _mm_store_pd(&C[34],Creg02) ;
      _mm_store_pd(&C[40],Creg10) ;
      _mm_store_pd(&C[42],Creg12) ;
      _mm_store_pd(&C[48],Creg20) ;
      _mm_store_pd(&C[50],Creg22) ;


    /*  Computing C(4:6,4:7)  */

      Areg0 = _mm_load_pd(&A[32]) ;
      Areg1 = _mm_load_pd(&A[40]) ;
      Areg2 = _mm_load_pd(&A[48]) ;

      Breg0 = _mm_load_pd(&B[4]) ;
      Breg2 = _mm_load_pd(&B[6]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(1,1)) ;
      Creg00 = Breg0 ;
      Creg00 = _mm_mul_pd(Creg00,Atmp) ;
      Creg02 = Breg2 ;
      Creg02 = _mm_mul_pd(Creg02,Atmp) ;

      Atmp = Areg1 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg1,_MM_SHUFFLE2(1,1)) ;
      Creg10 = Breg0 ;
      Creg10 = _mm_mul_pd(Creg10,Atmp) ;
      Creg12 = Breg2 ;
      Creg12 = _mm_mul_pd(Creg12,Atmp) ;

      Atmp = Areg2 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg2,_MM_SHUFFLE2(1,1)) ;
      Creg20 = Breg0 ;
      Creg20 = _mm_mul_pd(Creg20,Atmp) ;
      Creg22 = Breg2 ;
      Creg22 = _mm_mul_pd(Creg22,Atmp) ;


      Areg0 = _mm_load_pd(&A[34]) ;
      Areg1 = _mm_load_pd(&A[42]) ;
      Areg2 = _mm_load_pd(&A[50]) ;

      Breg0 = _mm_load_pd(&B[12]) ;
      Breg2 = _mm_load_pd(&B[14]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg00 = _mm_add_pd(Creg00,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg02 = _mm_add_pd(Creg02,Btmp) ;

      Atmp = Areg1 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg1,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg10 = _mm_add_pd(Creg10,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg12 = _mm_add_pd(Creg12,Btmp) ;

      Atmp = Areg2 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg2,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg20 = _mm_add_pd(Creg20,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg22 = _mm_add_pd(Creg22,Btmp) ;


      Breg0 = _mm_load_pd(&B[20]) ;
      Breg2 = _mm_load_pd(&B[22]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg00 = _mm_add_pd(Creg00,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg02 = _mm_add_pd(Creg02,Btmp) ;

      Atmp = Areg1 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg1,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg10 = _mm_add_pd(Creg10,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg12 = _mm_add_pd(Creg12,Btmp) ;

      Atmp = Areg2 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg2,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg20 = _mm_add_pd(Creg20,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg22 = _mm_add_pd(Creg22,Btmp) ;


      Areg0 = _mm_load_pd(&A[36]) ;
      Areg1 = _mm_load_pd(&A[44]) ;
      Areg2 = _mm_load_pd(&A[52]) ;

      Breg0 = _mm_load_pd(&B[28]) ;
      Breg2 = _mm_load_pd(&B[30]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg00 = _mm_add_pd(Creg00,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg02 = _mm_add_pd(Creg02,Btmp) ;

      Atmp = Areg1 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg1,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg10 = _mm_add_pd(Creg10,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg12 = _mm_add_pd(Creg12,Btmp) ;

      Atmp = Areg2 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg2,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg20 = _mm_add_pd(Creg20,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg22 = _mm_add_pd(Creg22,Btmp) ;


      Breg0 = _mm_load_pd(&B[36]) ;
      Breg2 = _mm_load_pd(&B[38]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg00 = _mm_add_pd(Creg00,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg02 = _mm_add_pd(Creg02,Btmp) ;

      Atmp = Areg1 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg1,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg10 = _mm_add_pd(Creg10,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg12 = _mm_add_pd(Creg12,Btmp) ;

      Atmp = Areg2 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg2,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg20 = _mm_add_pd(Creg20,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg22 = _mm_add_pd(Creg22,Btmp) ;


      Areg0 = _mm_load_pd(&A[38]) ;
      Areg1 = _mm_load_pd(&A[46]) ;
      Areg2 = _mm_load_pd(&A[54]) ;

      Breg0 = _mm_load_pd(&B[44]) ;
      Breg2 = _mm_load_pd(&B[46]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg00 = _mm_add_pd(Creg00,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg02 = _mm_add_pd(Creg02,Btmp) ;

      Atmp = Areg1 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg1,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg10 = _mm_add_pd(Creg10,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg12 = _mm_add_pd(Creg12,Btmp) ;

      Atmp = Areg2 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg2,_MM_SHUFFLE2(0,0)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg20 = _mm_add_pd(Creg20,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg22 = _mm_add_pd(Creg22,Btmp) ;


      Breg0 = _mm_load_pd(&B[52]) ;
      Breg2 = _mm_load_pd(&B[54]) ;

      Atmp = Areg0 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg0,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg00 = _mm_add_pd(Creg00,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg02 = _mm_add_pd(Creg02,Btmp) ;

      Atmp = Areg1 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg1,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg10 = _mm_add_pd(Creg10,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg12 = _mm_add_pd(Creg12,Btmp) ;

      Atmp = Areg2 ;
      Atmp = _mm_shuffle_pd(Atmp,Areg2,_MM_SHUFFLE2(1,1)) ;
      Btmp = Breg0 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg20 = _mm_add_pd(Creg20,Btmp) ;
      Btmp = Breg2 ;
      Btmp = _mm_mul_pd(Btmp,Atmp) ;
      Creg22 = _mm_add_pd(Creg22,Btmp) ;


      _mm_store_pd(&C[36],Creg00) ;
      _mm_store_pd(&C[38],Creg02) ;
      _mm_store_pd(&C[44],Creg10) ;
      _mm_store_pd(&C[46],Creg12) ;
      _mm_store_pd(&C[52],Creg20) ;
      _mm_store_pd(&C[54],Creg22) ;


}
