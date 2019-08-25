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

#include <pmmintrin.h>

void ttr_7_8_3_v2r2(const double * restrict A,
    const double * restrict B0,double * restrict tout0,
    const double * restrict B1,double * restrict tout1,
    const double * restrict B2,double * restrict tout2) {
__m128d Areg1,Areg2;
__m128d B0reg1,B0reg2,B1reg1,B1reg2,B2reg1,B2reg2;
__m128d T0reg1,T0reg2,T1reg1,T1reg2,T2reg1,T2reg2;

Areg1 = _mm_load_pd(&A[0]) ;
T0reg1 = _mm_load_pd(&B0[0]) ;
T1reg1 = _mm_load_pd(&B1[0]) ;
T2reg1 = _mm_load_pd(&B2[0]) ;
T0reg1 = _mm_mul_pd(T0reg1,Areg1) ;
T1reg1 = _mm_mul_pd(T1reg1,Areg1) ;
T2reg1 = _mm_mul_pd(T2reg1,Areg1) ;

Areg2 = _mm_load_pd(&A[2]) ;
T0reg2 = _mm_load_pd(&B0[2]) ;
T1reg2 = _mm_load_pd(&B1[2]) ;
T2reg2 = _mm_load_pd(&B2[2]) ;
T0reg2 = _mm_mul_pd(T0reg2,Areg2) ;
T1reg2 = _mm_mul_pd(T1reg2,Areg2) ;
T2reg2 = _mm_mul_pd(T2reg2,Areg2) ;

Areg1 = _mm_load_pd(&A[4]) ;
B0reg1 = _mm_load_pd(&B0[4]) ;
B1reg1 = _mm_load_pd(&B1[4]) ;
B2reg1 = _mm_load_pd(&B2[4]) ;
B0reg1 = _mm_mul_pd(B0reg1,Areg1) ;
T0reg1 = _mm_add_pd(T0reg1,B0reg1) ;
B1reg1 = _mm_mul_pd(B1reg1,Areg1) ;
T1reg1 = _mm_add_pd(T1reg1,B1reg1) ;
B2reg1 = _mm_mul_pd(B2reg1,Areg1) ;
T2reg1 = _mm_add_pd(T2reg1,B2reg1) ;

Areg2 = _mm_load_pd(&A[6]) ;
B0reg2 = _mm_load_pd(&B0[6]) ;
B1reg2 = _mm_load_pd(&B1[6]) ;
B2reg2 = _mm_load_pd(&B2[6]) ;
B0reg2 = _mm_mul_pd(B0reg2,Areg2) ;
T0reg2 = _mm_add_pd(T0reg2,B0reg2) ;
B1reg2 = _mm_mul_pd(B1reg2,Areg2) ;
T1reg2 = _mm_add_pd(T1reg2,B1reg2) ;
B2reg2 = _mm_mul_pd(B2reg2,Areg2) ;
T2reg2 = _mm_add_pd(T2reg2,B2reg2) ;

Areg1 = _mm_load_pd(&A[8]) ;
B0reg1 = _mm_load_pd(&B0[8]) ;
B1reg1 = _mm_load_pd(&B1[8]) ;
B2reg1 = _mm_load_pd(&B2[8]) ;
B0reg1 = _mm_mul_pd(B0reg1,Areg1) ;
T0reg1 = _mm_add_pd(T0reg1,B0reg1) ;
B1reg1 = _mm_mul_pd(B1reg1,Areg1) ;
T1reg1 = _mm_add_pd(T1reg1,B1reg1) ;
B2reg1 = _mm_mul_pd(B2reg1,Areg1) ;
T2reg1 = _mm_add_pd(T2reg1,B2reg1) ;

Areg2 = _mm_load_pd(&A[10]) ;
B0reg2 = _mm_load_pd(&B0[10]) ;
B1reg2 = _mm_load_pd(&B1[10]) ;
B2reg2 = _mm_load_pd(&B2[10]) ;
B0reg2 = _mm_mul_pd(B0reg2,Areg2) ;
T0reg2 = _mm_add_pd(T0reg2,B0reg2) ;
B1reg2 = _mm_mul_pd(B1reg2,Areg2) ;
T1reg2 = _mm_add_pd(T1reg2,B1reg2) ;
B2reg2 = _mm_mul_pd(B2reg2,Areg2) ;
T2reg2 = _mm_add_pd(T2reg2,B2reg2) ;

Areg1 = _mm_load_pd(&A[12]) ;
B0reg1 = _mm_load_pd(&B0[12]) ;
B1reg1 = _mm_load_pd(&B1[12]) ;
B2reg1 = _mm_load_pd(&B2[12]) ;
B0reg1 = _mm_mul_pd(B0reg1,Areg1) ;
T0reg1 = _mm_add_pd(T0reg1,B0reg1) ;
B1reg1 = _mm_mul_pd(B1reg1,Areg1) ;
T1reg1 = _mm_add_pd(T1reg1,B1reg1) ;
B2reg1 = _mm_mul_pd(B2reg1,Areg1) ;
T2reg1 = _mm_add_pd(T2reg1,B2reg1) ;

Areg2 = _mm_load_pd(&A[14]) ;
B0reg2 = _mm_load_pd(&B0[14]) ;
B1reg2 = _mm_load_pd(&B1[14]) ;
B2reg2 = _mm_load_pd(&B2[14]) ;
B0reg2 = _mm_mul_pd(B0reg2,Areg2) ;
T0reg2 = _mm_add_pd(T0reg2,B0reg2) ;
B1reg2 = _mm_mul_pd(B1reg2,Areg2) ;
T1reg2 = _mm_add_pd(T1reg2,B1reg2) ;
B2reg2 = _mm_mul_pd(B2reg2,Areg2) ;
T2reg2 = _mm_add_pd(T2reg2,B2reg2) ;

Areg1 = _mm_load_pd(&A[16]) ;
B0reg1 = _mm_load_pd(&B0[16]) ;
B1reg1 = _mm_load_pd(&B1[16]) ;
B2reg1 = _mm_load_pd(&B2[16]) ;
B0reg1 = _mm_mul_pd(B0reg1,Areg1) ;
T0reg1 = _mm_add_pd(T0reg1,B0reg1) ;
B1reg1 = _mm_mul_pd(B1reg1,Areg1) ;
T1reg1 = _mm_add_pd(T1reg1,B1reg1) ;
B2reg1 = _mm_mul_pd(B2reg1,Areg1) ;
T2reg1 = _mm_add_pd(T2reg1,B2reg1) ;

Areg2 = _mm_load_pd(&A[18]) ;
B0reg2 = _mm_load_pd(&B0[18]) ;
B1reg2 = _mm_load_pd(&B1[18]) ;
B2reg2 = _mm_load_pd(&B2[18]) ;
B0reg2 = _mm_mul_pd(B0reg2,Areg2) ;
T0reg2 = _mm_add_pd(T0reg2,B0reg2) ;
B1reg2 = _mm_mul_pd(B1reg2,Areg2) ;
T1reg2 = _mm_add_pd(T1reg2,B1reg2) ;
B2reg2 = _mm_mul_pd(B2reg2,Areg2) ;
T2reg2 = _mm_add_pd(T2reg2,B2reg2) ;

Areg1 = _mm_load_pd(&A[20]) ;
B0reg1 = _mm_load_pd(&B0[20]) ;
B1reg1 = _mm_load_pd(&B1[20]) ;
B2reg1 = _mm_load_pd(&B2[20]) ;
B0reg1 = _mm_mul_pd(B0reg1,Areg1) ;
T0reg1 = _mm_add_pd(T0reg1,B0reg1) ;
B1reg1 = _mm_mul_pd(B1reg1,Areg1) ;
T1reg1 = _mm_add_pd(T1reg1,B1reg1) ;
B2reg1 = _mm_mul_pd(B2reg1,Areg1) ;
T2reg1 = _mm_add_pd(T2reg1,B2reg1) ;

Areg2 = _mm_load_pd(&A[22]) ;
B0reg2 = _mm_load_pd(&B0[22]) ;
B1reg2 = _mm_load_pd(&B1[22]) ;
B2reg2 = _mm_load_pd(&B2[22]) ;
B0reg2 = _mm_mul_pd(B0reg2,Areg2) ;
T0reg2 = _mm_add_pd(T0reg2,B0reg2) ;
B1reg2 = _mm_mul_pd(B1reg2,Areg2) ;
T1reg2 = _mm_add_pd(T1reg2,B1reg2) ;
B2reg2 = _mm_mul_pd(B2reg2,Areg2) ;
T2reg2 = _mm_add_pd(T2reg2,B2reg2) ;

Areg1 = _mm_load_pd(&A[24]) ;
B0reg1 = _mm_load_pd(&B0[24]) ;
B1reg1 = _mm_load_pd(&B1[24]) ;
B2reg1 = _mm_load_pd(&B2[24]) ;
B0reg1 = _mm_mul_pd(B0reg1,Areg1) ;
T0reg1 = _mm_add_pd(T0reg1,B0reg1) ;
B1reg1 = _mm_mul_pd(B1reg1,Areg1) ;
T1reg1 = _mm_add_pd(T1reg1,B1reg1) ;
B2reg1 = _mm_mul_pd(B2reg1,Areg1) ;
T2reg1 = _mm_add_pd(T2reg1,B2reg1) ;

Areg2 = _mm_load_pd(&A[26]) ;
B0reg2 = _mm_load_pd(&B0[26]) ;
B1reg2 = _mm_load_pd(&B1[26]) ;
B2reg2 = _mm_load_pd(&B2[26]) ;
B0reg2 = _mm_mul_pd(B0reg2,Areg2) ;
T0reg2 = _mm_add_pd(T0reg2,B0reg2) ;
B1reg2 = _mm_mul_pd(B1reg2,Areg2) ;
T1reg2 = _mm_add_pd(T1reg2,B1reg2) ;
B2reg2 = _mm_mul_pd(B2reg2,Areg2) ;
T2reg2 = _mm_add_pd(T2reg2,B2reg2) ;

Areg1 = _mm_load_pd(&A[28]) ;
B0reg1 = _mm_load_pd(&B0[28]) ;
B1reg1 = _mm_load_pd(&B1[28]) ;
B2reg1 = _mm_load_pd(&B2[28]) ;
B0reg1 = _mm_mul_pd(B0reg1,Areg1) ;
T0reg1 = _mm_add_pd(T0reg1,B0reg1) ;
B1reg1 = _mm_mul_pd(B1reg1,Areg1) ;
T1reg1 = _mm_add_pd(T1reg1,B1reg1) ;
B2reg1 = _mm_mul_pd(B2reg1,Areg1) ;
T2reg1 = _mm_add_pd(T2reg1,B2reg1) ;

Areg2 = _mm_load_pd(&A[30]) ;
B0reg2 = _mm_load_pd(&B0[30]) ;
B1reg2 = _mm_load_pd(&B1[30]) ;
B2reg2 = _mm_load_pd(&B2[30]) ;
B0reg2 = _mm_mul_pd(B0reg2,Areg2) ;
T0reg2 = _mm_add_pd(T0reg2,B0reg2) ;
B1reg2 = _mm_mul_pd(B1reg2,Areg2) ;
T1reg2 = _mm_add_pd(T1reg2,B1reg2) ;
B2reg2 = _mm_mul_pd(B2reg2,Areg2) ;
T2reg2 = _mm_add_pd(T2reg2,B2reg2) ;

Areg1 = _mm_load_pd(&A[32]) ;
B0reg1 = _mm_load_pd(&B0[32]) ;
B1reg1 = _mm_load_pd(&B1[32]) ;
B2reg1 = _mm_load_pd(&B2[32]) ;
B0reg1 = _mm_mul_pd(B0reg1,Areg1) ;
T0reg1 = _mm_add_pd(T0reg1,B0reg1) ;
B1reg1 = _mm_mul_pd(B1reg1,Areg1) ;
T1reg1 = _mm_add_pd(T1reg1,B1reg1) ;
B2reg1 = _mm_mul_pd(B2reg1,Areg1) ;
T2reg1 = _mm_add_pd(T2reg1,B2reg1) ;

Areg2 = _mm_load_pd(&A[34]) ;
B0reg2 = _mm_load_pd(&B0[34]) ;
B1reg2 = _mm_load_pd(&B1[34]) ;
B2reg2 = _mm_load_pd(&B2[34]) ;
B0reg2 = _mm_mul_pd(B0reg2,Areg2) ;
T0reg2 = _mm_add_pd(T0reg2,B0reg2) ;
B1reg2 = _mm_mul_pd(B1reg2,Areg2) ;
T1reg2 = _mm_add_pd(T1reg2,B1reg2) ;
B2reg2 = _mm_mul_pd(B2reg2,Areg2) ;
T2reg2 = _mm_add_pd(T2reg2,B2reg2) ;

Areg1 = _mm_load_pd(&A[36]) ;
B0reg1 = _mm_load_pd(&B0[36]) ;
B1reg1 = _mm_load_pd(&B1[36]) ;
B2reg1 = _mm_load_pd(&B2[36]) ;
B0reg1 = _mm_mul_pd(B0reg1,Areg1) ;
T0reg1 = _mm_add_pd(T0reg1,B0reg1) ;
B1reg1 = _mm_mul_pd(B1reg1,Areg1) ;
T1reg1 = _mm_add_pd(T1reg1,B1reg1) ;
B2reg1 = _mm_mul_pd(B2reg1,Areg1) ;
T2reg1 = _mm_add_pd(T2reg1,B2reg1) ;

Areg2 = _mm_load_pd(&A[38]) ;
B0reg2 = _mm_load_pd(&B0[38]) ;
B1reg2 = _mm_load_pd(&B1[38]) ;
B2reg2 = _mm_load_pd(&B2[38]) ;
B0reg2 = _mm_mul_pd(B0reg2,Areg2) ;
T0reg2 = _mm_add_pd(T0reg2,B0reg2) ;
B1reg2 = _mm_mul_pd(B1reg2,Areg2) ;
T1reg2 = _mm_add_pd(T1reg2,B1reg2) ;
B2reg2 = _mm_mul_pd(B2reg2,Areg2) ;
T2reg2 = _mm_add_pd(T2reg2,B2reg2) ;

Areg1 = _mm_load_pd(&A[40]) ;
B0reg1 = _mm_load_pd(&B0[40]) ;
B1reg1 = _mm_load_pd(&B1[40]) ;
B2reg1 = _mm_load_pd(&B2[40]) ;
B0reg1 = _mm_mul_pd(B0reg1,Areg1) ;
T0reg1 = _mm_add_pd(T0reg1,B0reg1) ;
B1reg1 = _mm_mul_pd(B1reg1,Areg1) ;
T1reg1 = _mm_add_pd(T1reg1,B1reg1) ;
B2reg1 = _mm_mul_pd(B2reg1,Areg1) ;
T2reg1 = _mm_add_pd(T2reg1,B2reg1) ;

Areg2 = _mm_load_pd(&A[42]) ;
B0reg2 = _mm_load_pd(&B0[42]) ;
B1reg2 = _mm_load_pd(&B1[42]) ;
B2reg2 = _mm_load_pd(&B2[42]) ;
B0reg2 = _mm_mul_pd(B0reg2,Areg2) ;
T0reg2 = _mm_add_pd(T0reg2,B0reg2) ;
B1reg2 = _mm_mul_pd(B1reg2,Areg2) ;
T1reg2 = _mm_add_pd(T1reg2,B1reg2) ;
B2reg2 = _mm_mul_pd(B2reg2,Areg2) ;
T2reg2 = _mm_add_pd(T2reg2,B2reg2) ;

Areg1 = _mm_load_pd(&A[44]) ;
B0reg1 = _mm_load_pd(&B0[44]) ;
B1reg1 = _mm_load_pd(&B1[44]) ;
B2reg1 = _mm_load_pd(&B2[44]) ;
B0reg1 = _mm_mul_pd(B0reg1,Areg1) ;
T0reg1 = _mm_add_pd(T0reg1,B0reg1) ;
B1reg1 = _mm_mul_pd(B1reg1,Areg1) ;
T1reg1 = _mm_add_pd(T1reg1,B1reg1) ;
B2reg1 = _mm_mul_pd(B2reg1,Areg1) ;
T2reg1 = _mm_add_pd(T2reg1,B2reg1) ;

Areg2 = _mm_load_pd(&A[46]) ;
B0reg2 = _mm_load_pd(&B0[46]) ;
B1reg2 = _mm_load_pd(&B1[46]) ;
B2reg2 = _mm_load_pd(&B2[46]) ;
B0reg2 = _mm_mul_pd(B0reg2,Areg2) ;
T0reg2 = _mm_add_pd(T0reg2,B0reg2) ;
B1reg2 = _mm_mul_pd(B1reg2,Areg2) ;
T1reg2 = _mm_add_pd(T1reg2,B1reg2) ;
B2reg2 = _mm_mul_pd(B2reg2,Areg2) ;
T2reg2 = _mm_add_pd(T2reg2,B2reg2) ;

Areg1 = _mm_load_pd(&A[48]) ;
B0reg1 = _mm_load_pd(&B0[48]) ;
B1reg1 = _mm_load_pd(&B1[48]) ;
B2reg1 = _mm_load_pd(&B2[48]) ;
B0reg1 = _mm_mul_pd(B0reg1,Areg1) ;
T0reg1 = _mm_add_pd(T0reg1,B0reg1) ;
B1reg1 = _mm_mul_pd(B1reg1,Areg1) ;
T1reg1 = _mm_add_pd(T1reg1,B1reg1) ;
B2reg1 = _mm_mul_pd(B2reg1,Areg1) ;
T2reg1 = _mm_add_pd(T2reg1,B2reg1) ;

Areg2 = _mm_load_pd(&A[50]) ;
B0reg2 = _mm_load_pd(&B0[50]) ;
B1reg2 = _mm_load_pd(&B1[50]) ;
B2reg2 = _mm_load_pd(&B2[50]) ;
B0reg2 = _mm_mul_pd(B0reg2,Areg2) ;
T0reg2 = _mm_add_pd(T0reg2,B0reg2) ;
B1reg2 = _mm_mul_pd(B1reg2,Areg2) ;
T1reg2 = _mm_add_pd(T1reg2,B1reg2) ;
B2reg2 = _mm_mul_pd(B2reg2,Areg2) ;
T2reg2 = _mm_add_pd(T2reg2,B2reg2) ;

Areg1 = _mm_load_pd(&A[52]) ;
B0reg1 = _mm_load_pd(&B0[52]) ;
B1reg1 = _mm_load_pd(&B1[52]) ;
B2reg1 = _mm_load_pd(&B2[52]) ;
B0reg1 = _mm_mul_pd(B0reg1,Areg1) ;
T0reg1 = _mm_add_pd(T0reg1,B0reg1) ;
B1reg1 = _mm_mul_pd(B1reg1,Areg1) ;
T1reg1 = _mm_add_pd(T1reg1,B1reg1) ;
B2reg1 = _mm_mul_pd(B2reg1,Areg1) ;
T2reg1 = _mm_add_pd(T2reg1,B2reg1) ;

Areg2 = _mm_load_pd(&A[54]) ;
B0reg2 = _mm_load_pd(&B0[54]) ;
B1reg2 = _mm_load_pd(&B1[54]) ;
B2reg2 = _mm_load_pd(&B2[54]) ;
B0reg2 = _mm_mul_pd(B0reg2,Areg2) ;
T0reg2 = _mm_add_pd(T0reg2,B0reg2) ;
B1reg2 = _mm_mul_pd(B1reg2,Areg2) ;
T1reg2 = _mm_add_pd(T1reg2,B1reg2) ;
B2reg2 = _mm_mul_pd(B2reg2,Areg2) ;
T2reg2 = _mm_add_pd(T2reg2,B2reg2) ;


T0reg1 = _mm_add_pd(T0reg1,T0reg2) ;
T1reg1 = _mm_add_pd(T1reg1,T1reg2) ;
T2reg1 = _mm_add_pd(T2reg1,T2reg2) ;

T0reg1 = _mm_hadd_pd(T0reg1,T0reg1) ;
_mm_store_sd(tout0,T0reg1) ;
T1reg1 = _mm_hadd_pd(T1reg1,T1reg1) ;
_mm_store_sd(tout1,T1reg1) ;
T2reg1 = _mm_hadd_pd(T2reg1,T2reg1) ;
_mm_store_sd(tout2,T2reg1) ;

}
