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


//#define TESTING

#ifdef TESTING

typedef struct {
  double x,y;
} pair;


#define CPLX pair

static double my_creal(CPLX x) {
  return ((double *) &x)[0];
}
static double my_cimag(CPLX x) {
  return ((double *) &x)[1];
}
static CPLX my_lfpd(const double *p) {
  return ((CPLX *) p)[0];
}

/*
// Not needed for trace calculations.
static void my_stfpd(double *p,CPLX x) {
  ((CPLX *) p)[0] = x;
}
static CPLX my_fxsmul(CPLX x,double a) {
  double y[2];
  y[0] = a * my_creal(x);
  y[1] = a * my_cimag(x);
  return ((CPLX *) y)[0];
}
static CPLX my_fxcpmadd(CPLX t,CPLX x,double a) {
  double y[2];
  y[0] = my_creal(t) + a * my_creal(x);
  y[1] = my_cimag(t) + a * my_cimag(x);
  return ((CPLX *) y)[0];
}
*/

static CPLX my_fpmul(CPLX x,CPLX y) {
  union {
    double z[2];
    CPLX c;
  } U;
  U.z[0] = my_creal(y) * my_creal(x);
  U.z[1] = my_cimag(y) * my_cimag(x);
  return U.c;
}
static CPLX my_fpadd(CPLX x,CPLX y) {
  union {
    double z[2];
    CPLX c;
  } U;
  U.z[0] = my_creal(y) + my_creal(x);
  U.z[1] = my_cimag(y) + my_cimag(x);
  return U.c;
}
static CPLX my_fpmadd(CPLX t,CPLX x,CPLX y) {
  union {
    double z[2];
    CPLX c;
  } U;
  U.z[0] = my_creal(t) + my_creal(y) * my_creal(x);
  U.z[1] = my_cimag(t) + my_cimag(y) * my_cimag(x);
  return U.c;
}

#define __creal    my_creal
#define __cimag    my_cimag
#define __lfpd     my_lfpd
#define __stfpd    my_stfpd
#define __fxsmul   my_fxsmul
#define __fxcpmadd my_fxcpmadd

#define __fpadd my_fpadd
#define __fpmul my_fpmul
#define __fpmadd my_fpmadd

#else

#define CPLX double _Complex

#endif

void ttr_bg_7_8_3_v2r3(const double * restrict A,
    const double * restrict B0,double * restrict tout0,
    const double * restrict B1,double * restrict tout1,
    const double * restrict B2,double * restrict tout2) {
CPLX Areg1,Areg2,Areg3;
CPLX B0reg1,B0reg2,B0reg3,B1reg1,B1reg2,B1reg3,B2reg1,B2reg2,B2reg3;
CPLX T0reg1,T0reg2,T0reg3,T1reg1,T1reg2,T1reg3,T2reg1,T2reg2,T2reg3;

Areg1 = __lfpd(&A[0]) ;
B0reg1 = __lfpd(&B0[0]) ;
B1reg1 = __lfpd(&B1[0]) ;
B2reg1 = __lfpd(&B2[0]) ;
T0reg1 = __fpmul(Areg1,B0reg1) ;
T1reg1 = __fpmul(Areg1,B1reg1) ;
T2reg1 = __fpmul(Areg1,B2reg1) ;

Areg2 = __lfpd(&A[2]) ;
B0reg2 = __lfpd(&B0[2]) ;
B1reg2 = __lfpd(&B1[2]) ;
B2reg2 = __lfpd(&B2[2]) ;
T0reg2 = __fpmul(Areg2,B0reg2) ;
T1reg2 = __fpmul(Areg2,B1reg2) ;
T2reg2 = __fpmul(Areg2,B2reg2) ;

Areg3 = __lfpd(&A[4]) ;
B0reg3 = __lfpd(&B0[4]) ;
B1reg3 = __lfpd(&B1[4]) ;
B2reg3 = __lfpd(&B2[4]) ;
T0reg3 = __fpmul(Areg3,B0reg3) ;
T1reg3 = __fpmul(Areg3,B1reg3) ;
T2reg3 = __fpmul(Areg3,B2reg3) ;

Areg1 = __lfpd(&A[6]) ;
B0reg1 = __lfpd(&B0[6]) ;
B1reg1 = __lfpd(&B1[6]) ;
B2reg1 = __lfpd(&B2[6]) ;
T0reg1 = __fpmadd(T0reg1,Areg1,B0reg1) ;
T1reg1 = __fpmadd(T1reg1,Areg1,B1reg1) ;
T2reg1 = __fpmadd(T2reg1,Areg1,B2reg1) ;

Areg2 = __lfpd(&A[8]) ;
B0reg2 = __lfpd(&B0[8]) ;
B1reg2 = __lfpd(&B1[8]) ;
B2reg2 = __lfpd(&B2[8]) ;
T0reg2 = __fpmadd(T0reg2,Areg2,B0reg2) ;
T1reg2 = __fpmadd(T1reg2,Areg2,B1reg2) ;
T2reg2 = __fpmadd(T2reg2,Areg2,B2reg2) ;

Areg3 = __lfpd(&A[10]) ;
B0reg3 = __lfpd(&B0[10]) ;
B1reg3 = __lfpd(&B1[10]) ;
B2reg3 = __lfpd(&B2[10]) ;
T0reg3 = __fpmadd(T0reg3,Areg3,B0reg3) ;
T1reg3 = __fpmadd(T1reg3,Areg3,B1reg3) ;
T2reg3 = __fpmadd(T2reg3,Areg3,B2reg3) ;

Areg1 = __lfpd(&A[12]) ;
B0reg1 = __lfpd(&B0[12]) ;
B1reg1 = __lfpd(&B1[12]) ;
B2reg1 = __lfpd(&B2[12]) ;
T0reg1 = __fpmadd(T0reg1,Areg1,B0reg1) ;
T1reg1 = __fpmadd(T1reg1,Areg1,B1reg1) ;
T2reg1 = __fpmadd(T2reg1,Areg1,B2reg1) ;

Areg2 = __lfpd(&A[14]) ;
B0reg2 = __lfpd(&B0[14]) ;
B1reg2 = __lfpd(&B1[14]) ;
B2reg2 = __lfpd(&B2[14]) ;
T0reg2 = __fpmadd(T0reg2,Areg2,B0reg2) ;
T1reg2 = __fpmadd(T1reg2,Areg2,B1reg2) ;
T2reg2 = __fpmadd(T2reg2,Areg2,B2reg2) ;

Areg3 = __lfpd(&A[16]) ;
B0reg3 = __lfpd(&B0[16]) ;
B1reg3 = __lfpd(&B1[16]) ;
B2reg3 = __lfpd(&B2[16]) ;
T0reg3 = __fpmadd(T0reg3,Areg3,B0reg3) ;
T1reg3 = __fpmadd(T1reg3,Areg3,B1reg3) ;
T2reg3 = __fpmadd(T2reg3,Areg3,B2reg3) ;

Areg1 = __lfpd(&A[18]) ;
B0reg1 = __lfpd(&B0[18]) ;
B1reg1 = __lfpd(&B1[18]) ;
B2reg1 = __lfpd(&B2[18]) ;
T0reg1 = __fpmadd(T0reg1,Areg1,B0reg1) ;
T1reg1 = __fpmadd(T1reg1,Areg1,B1reg1) ;
T2reg1 = __fpmadd(T2reg1,Areg1,B2reg1) ;

Areg2 = __lfpd(&A[20]) ;
B0reg2 = __lfpd(&B0[20]) ;
B1reg2 = __lfpd(&B1[20]) ;
B2reg2 = __lfpd(&B2[20]) ;
T0reg2 = __fpmadd(T0reg2,Areg2,B0reg2) ;
T1reg2 = __fpmadd(T1reg2,Areg2,B1reg2) ;
T2reg2 = __fpmadd(T2reg2,Areg2,B2reg2) ;

Areg3 = __lfpd(&A[22]) ;
B0reg3 = __lfpd(&B0[22]) ;
B1reg3 = __lfpd(&B1[22]) ;
B2reg3 = __lfpd(&B2[22]) ;
T0reg3 = __fpmadd(T0reg3,Areg3,B0reg3) ;
T1reg3 = __fpmadd(T1reg3,Areg3,B1reg3) ;
T2reg3 = __fpmadd(T2reg3,Areg3,B2reg3) ;

Areg1 = __lfpd(&A[24]) ;
B0reg1 = __lfpd(&B0[24]) ;
B1reg1 = __lfpd(&B1[24]) ;
B2reg1 = __lfpd(&B2[24]) ;
T0reg1 = __fpmadd(T0reg1,Areg1,B0reg1) ;
T1reg1 = __fpmadd(T1reg1,Areg1,B1reg1) ;
T2reg1 = __fpmadd(T2reg1,Areg1,B2reg1) ;

Areg2 = __lfpd(&A[26]) ;
B0reg2 = __lfpd(&B0[26]) ;
B1reg2 = __lfpd(&B1[26]) ;
B2reg2 = __lfpd(&B2[26]) ;
T0reg2 = __fpmadd(T0reg2,Areg2,B0reg2) ;
T1reg2 = __fpmadd(T1reg2,Areg2,B1reg2) ;
T2reg2 = __fpmadd(T2reg2,Areg2,B2reg2) ;

Areg3 = __lfpd(&A[28]) ;
B0reg3 = __lfpd(&B0[28]) ;
B1reg3 = __lfpd(&B1[28]) ;
B2reg3 = __lfpd(&B2[28]) ;
T0reg3 = __fpmadd(T0reg3,Areg3,B0reg3) ;
T1reg3 = __fpmadd(T1reg3,Areg3,B1reg3) ;
T2reg3 = __fpmadd(T2reg3,Areg3,B2reg3) ;

Areg1 = __lfpd(&A[30]) ;
B0reg1 = __lfpd(&B0[30]) ;
B1reg1 = __lfpd(&B1[30]) ;
B2reg1 = __lfpd(&B2[30]) ;
T0reg1 = __fpmadd(T0reg1,Areg1,B0reg1) ;
T1reg1 = __fpmadd(T1reg1,Areg1,B1reg1) ;
T2reg1 = __fpmadd(T2reg1,Areg1,B2reg1) ;

Areg2 = __lfpd(&A[32]) ;
B0reg2 = __lfpd(&B0[32]) ;
B1reg2 = __lfpd(&B1[32]) ;
B2reg2 = __lfpd(&B2[32]) ;
T0reg2 = __fpmadd(T0reg2,Areg2,B0reg2) ;
T1reg2 = __fpmadd(T1reg2,Areg2,B1reg2) ;
T2reg2 = __fpmadd(T2reg2,Areg2,B2reg2) ;

Areg3 = __lfpd(&A[34]) ;
B0reg3 = __lfpd(&B0[34]) ;
B1reg3 = __lfpd(&B1[34]) ;
B2reg3 = __lfpd(&B2[34]) ;
T0reg3 = __fpmadd(T0reg3,Areg3,B0reg3) ;
T1reg3 = __fpmadd(T1reg3,Areg3,B1reg3) ;
T2reg3 = __fpmadd(T2reg3,Areg3,B2reg3) ;

Areg1 = __lfpd(&A[36]) ;
B0reg1 = __lfpd(&B0[36]) ;
B1reg1 = __lfpd(&B1[36]) ;
B2reg1 = __lfpd(&B2[36]) ;
T0reg1 = __fpmadd(T0reg1,Areg1,B0reg1) ;
T1reg1 = __fpmadd(T1reg1,Areg1,B1reg1) ;
T2reg1 = __fpmadd(T2reg1,Areg1,B2reg1) ;

Areg2 = __lfpd(&A[38]) ;
B0reg2 = __lfpd(&B0[38]) ;
B1reg2 = __lfpd(&B1[38]) ;
B2reg2 = __lfpd(&B2[38]) ;
T0reg2 = __fpmadd(T0reg2,Areg2,B0reg2) ;
T1reg2 = __fpmadd(T1reg2,Areg2,B1reg2) ;
T2reg2 = __fpmadd(T2reg2,Areg2,B2reg2) ;

Areg3 = __lfpd(&A[40]) ;
B0reg3 = __lfpd(&B0[40]) ;
B1reg3 = __lfpd(&B1[40]) ;
B2reg3 = __lfpd(&B2[40]) ;
T0reg3 = __fpmadd(T0reg3,Areg3,B0reg3) ;
T1reg3 = __fpmadd(T1reg3,Areg3,B1reg3) ;
T2reg3 = __fpmadd(T2reg3,Areg3,B2reg3) ;

Areg1 = __lfpd(&A[42]) ;
B0reg1 = __lfpd(&B0[42]) ;
B1reg1 = __lfpd(&B1[42]) ;
B2reg1 = __lfpd(&B2[42]) ;
T0reg1 = __fpmadd(T0reg1,Areg1,B0reg1) ;
T1reg1 = __fpmadd(T1reg1,Areg1,B1reg1) ;
T2reg1 = __fpmadd(T2reg1,Areg1,B2reg1) ;

Areg2 = __lfpd(&A[44]) ;
B0reg2 = __lfpd(&B0[44]) ;
B1reg2 = __lfpd(&B1[44]) ;
B2reg2 = __lfpd(&B2[44]) ;
T0reg2 = __fpmadd(T0reg2,Areg2,B0reg2) ;
T1reg2 = __fpmadd(T1reg2,Areg2,B1reg2) ;
T2reg2 = __fpmadd(T2reg2,Areg2,B2reg2) ;

Areg3 = __lfpd(&A[46]) ;
B0reg3 = __lfpd(&B0[46]) ;
B1reg3 = __lfpd(&B1[46]) ;
B2reg3 = __lfpd(&B2[46]) ;
T0reg3 = __fpmadd(T0reg3,Areg3,B0reg3) ;
T1reg3 = __fpmadd(T1reg3,Areg3,B1reg3) ;
T2reg3 = __fpmadd(T2reg3,Areg3,B2reg3) ;

Areg1 = __lfpd(&A[48]) ;
B0reg1 = __lfpd(&B0[48]) ;
B1reg1 = __lfpd(&B1[48]) ;
B2reg1 = __lfpd(&B2[48]) ;
T0reg1 = __fpmadd(T0reg1,Areg1,B0reg1) ;
T1reg1 = __fpmadd(T1reg1,Areg1,B1reg1) ;
T2reg1 = __fpmadd(T2reg1,Areg1,B2reg1) ;

Areg2 = __lfpd(&A[50]) ;
B0reg2 = __lfpd(&B0[50]) ;
B1reg2 = __lfpd(&B1[50]) ;
B2reg2 = __lfpd(&B2[50]) ;
T0reg2 = __fpmadd(T0reg2,Areg2,B0reg2) ;
T1reg2 = __fpmadd(T1reg2,Areg2,B1reg2) ;
T2reg2 = __fpmadd(T2reg2,Areg2,B2reg2) ;

Areg3 = __lfpd(&A[52]) ;
B0reg3 = __lfpd(&B0[52]) ;
B1reg3 = __lfpd(&B1[52]) ;
B2reg3 = __lfpd(&B2[52]) ;
T0reg3 = __fpmadd(T0reg3,Areg3,B0reg3) ;
T1reg3 = __fpmadd(T1reg3,Areg3,B1reg3) ;
T2reg3 = __fpmadd(T2reg3,Areg3,B2reg3) ;

Areg1 = __lfpd(&A[54]) ;
B0reg1 = __lfpd(&B0[54]) ;
B1reg1 = __lfpd(&B1[54]) ;
B2reg1 = __lfpd(&B2[54]) ;
T0reg1 = __fpmadd(T0reg1,Areg1,B0reg1) ;
T1reg1 = __fpmadd(T1reg1,Areg1,B1reg1) ;
T2reg1 = __fpmadd(T2reg1,Areg1,B2reg1) ;


T0reg1 = __fpadd(T0reg1,T0reg2) ;
T1reg1 = __fpadd(T1reg1,T1reg2) ;
T2reg1 = __fpadd(T2reg1,T2reg2) ;
T0reg1 = __fpadd(T0reg1,T0reg3) ;
T1reg1 = __fpadd(T1reg1,T1reg3) ;
T2reg1 = __fpadd(T2reg1,T2reg3) ;

*tout0 = __creal(T0reg1) + __cimag(T0reg1) ;
*tout1 = __creal(T1reg1) + __cimag(T1reg1) ;
*tout2 = __creal(T2reg1) + __cimag(T2reg1) ;

}
