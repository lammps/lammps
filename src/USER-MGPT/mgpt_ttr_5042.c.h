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

void ttr_5_8_3_v4r2(const float * restrict A,
    const float * restrict B0,float * restrict tout0,
    const float * restrict B1,float * restrict tout1,
    const float * restrict B2,float * restrict tout2) {
__m128 Areg1,Areg2;
__m128 B0reg1,B0reg2,B1reg1,B1reg2,B2reg1,B2reg2;
__m128 T0reg1,T0reg2,T1reg1,T1reg2,T2reg1,T2reg2;

Areg1 = _mm_load_ps(&A[0]) ;
T0reg1 = _mm_load_ps(&B0[0]) ;
T1reg1 = _mm_load_ps(&B1[0]) ;
T2reg1 = _mm_load_ps(&B2[0]) ;
T0reg1 = _mm_mul_ps(T0reg1,Areg1) ;
T1reg1 = _mm_mul_ps(T1reg1,Areg1) ;
T2reg1 = _mm_mul_ps(T2reg1,Areg1) ;

Areg2 = _mm_load_ps(&A[4]) ;
T0reg2 = _mm_load_ps(&B0[4]) ;
T1reg2 = _mm_load_ps(&B1[4]) ;
T2reg2 = _mm_load_ps(&B2[4]) ;
T0reg2 = _mm_mul_ps(T0reg2,Areg2) ;
T1reg2 = _mm_mul_ps(T1reg2,Areg2) ;
T2reg2 = _mm_mul_ps(T2reg2,Areg2) ;

Areg1 = _mm_load_ps(&A[8]) ;
B0reg1 = _mm_load_ps(&B0[8]) ;
B1reg1 = _mm_load_ps(&B1[8]) ;
B2reg1 = _mm_load_ps(&B2[8]) ;
B0reg1 = _mm_mul_ps(B0reg1,Areg1) ;
T0reg1 = _mm_add_ps(T0reg1,B0reg1) ;
B1reg1 = _mm_mul_ps(B1reg1,Areg1) ;
T1reg1 = _mm_add_ps(T1reg1,B1reg1) ;
B2reg1 = _mm_mul_ps(B2reg1,Areg1) ;
T2reg1 = _mm_add_ps(T2reg1,B2reg1) ;

Areg2 = _mm_load_ps(&A[12]) ;
B0reg2 = _mm_load_ps(&B0[12]) ;
B1reg2 = _mm_load_ps(&B1[12]) ;
B2reg2 = _mm_load_ps(&B2[12]) ;
B0reg2 = _mm_mul_ps(B0reg2,Areg2) ;
T0reg2 = _mm_add_ps(T0reg2,B0reg2) ;
B1reg2 = _mm_mul_ps(B1reg2,Areg2) ;
T1reg2 = _mm_add_ps(T1reg2,B1reg2) ;
B2reg2 = _mm_mul_ps(B2reg2,Areg2) ;
T2reg2 = _mm_add_ps(T2reg2,B2reg2) ;

Areg1 = _mm_load_ps(&A[16]) ;
B0reg1 = _mm_load_ps(&B0[16]) ;
B1reg1 = _mm_load_ps(&B1[16]) ;
B2reg1 = _mm_load_ps(&B2[16]) ;
B0reg1 = _mm_mul_ps(B0reg1,Areg1) ;
T0reg1 = _mm_add_ps(T0reg1,B0reg1) ;
B1reg1 = _mm_mul_ps(B1reg1,Areg1) ;
T1reg1 = _mm_add_ps(T1reg1,B1reg1) ;
B2reg1 = _mm_mul_ps(B2reg1,Areg1) ;
T2reg1 = _mm_add_ps(T2reg1,B2reg1) ;

Areg2 = _mm_load_ps(&A[20]) ;
B0reg2 = _mm_load_ps(&B0[20]) ;
B1reg2 = _mm_load_ps(&B1[20]) ;
B2reg2 = _mm_load_ps(&B2[20]) ;
B0reg2 = _mm_mul_ps(B0reg2,Areg2) ;
T0reg2 = _mm_add_ps(T0reg2,B0reg2) ;
B1reg2 = _mm_mul_ps(B1reg2,Areg2) ;
T1reg2 = _mm_add_ps(T1reg2,B1reg2) ;
B2reg2 = _mm_mul_ps(B2reg2,Areg2) ;
T2reg2 = _mm_add_ps(T2reg2,B2reg2) ;

Areg1 = _mm_load_ps(&A[24]) ;
B0reg1 = _mm_load_ps(&B0[24]) ;
B1reg1 = _mm_load_ps(&B1[24]) ;
B2reg1 = _mm_load_ps(&B2[24]) ;
B0reg1 = _mm_mul_ps(B0reg1,Areg1) ;
T0reg1 = _mm_add_ps(T0reg1,B0reg1) ;
B1reg1 = _mm_mul_ps(B1reg1,Areg1) ;
T1reg1 = _mm_add_ps(T1reg1,B1reg1) ;
B2reg1 = _mm_mul_ps(B2reg1,Areg1) ;
T2reg1 = _mm_add_ps(T2reg1,B2reg1) ;

Areg2 = _mm_load_ps(&A[28]) ;
B0reg2 = _mm_load_ps(&B0[28]) ;
B1reg2 = _mm_load_ps(&B1[28]) ;
B2reg2 = _mm_load_ps(&B2[28]) ;
B0reg2 = _mm_mul_ps(B0reg2,Areg2) ;
T0reg2 = _mm_add_ps(T0reg2,B0reg2) ;
B1reg2 = _mm_mul_ps(B1reg2,Areg2) ;
T1reg2 = _mm_add_ps(T1reg2,B1reg2) ;
B2reg2 = _mm_mul_ps(B2reg2,Areg2) ;
T2reg2 = _mm_add_ps(T2reg2,B2reg2) ;

Areg1 = _mm_load_ps(&A[32]) ;
B0reg1 = _mm_load_ps(&B0[32]) ;
B1reg1 = _mm_load_ps(&B1[32]) ;
B2reg1 = _mm_load_ps(&B2[32]) ;
B0reg1 = _mm_mul_ps(B0reg1,Areg1) ;
T0reg1 = _mm_add_ps(T0reg1,B0reg1) ;
B1reg1 = _mm_mul_ps(B1reg1,Areg1) ;
T1reg1 = _mm_add_ps(T1reg1,B1reg1) ;
B2reg1 = _mm_mul_ps(B2reg1,Areg1) ;
T2reg1 = _mm_add_ps(T2reg1,B2reg1) ;

Areg2 = _mm_load_ps(&A[36]) ;
B0reg2 = _mm_load_ps(&B0[36]) ;
B1reg2 = _mm_load_ps(&B1[36]) ;
B2reg2 = _mm_load_ps(&B2[36]) ;
B0reg2 = _mm_mul_ps(B0reg2,Areg2) ;
T0reg2 = _mm_add_ps(T0reg2,B0reg2) ;
B1reg2 = _mm_mul_ps(B1reg2,Areg2) ;
T1reg2 = _mm_add_ps(T1reg2,B1reg2) ;
B2reg2 = _mm_mul_ps(B2reg2,Areg2) ;
T2reg2 = _mm_add_ps(T2reg2,B2reg2) ;


T0reg1 = _mm_add_ps(T0reg1,T0reg2) ;
T1reg1 = _mm_add_ps(T1reg1,T1reg2) ;
T2reg1 = _mm_add_ps(T2reg1,T2reg2) ;

T0reg1 = _mm_hadd_ps(T0reg1,T0reg1) ;
T0reg1 = _mm_hadd_ps(T0reg1,T0reg1) ;
_mm_store_ss(tout0,T0reg1) ;
T1reg1 = _mm_hadd_ps(T1reg1,T1reg1) ;
T1reg1 = _mm_hadd_ps(T1reg1,T1reg1) ;
_mm_store_ss(tout1,T1reg1) ;
T2reg1 = _mm_hadd_ps(T2reg1,T2reg1) ;
T2reg1 = _mm_hadd_ps(T2reg1,T2reg1) ;
_mm_store_ss(tout2,T2reg1) ;

}
