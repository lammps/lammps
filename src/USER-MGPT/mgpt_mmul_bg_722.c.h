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


/* #define TESTING */

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

#define __creal    my_creal
#define __cimag    my_cimag
#define __lfpd     my_lfpd
#define __stfpd    my_stfpd
#define __fxsmul   my_fxsmul
#define __fxcpmadd my_fxcpmadd

#else

#define CPLX double _Complex

#endif

void mmul_bg_7_8_2x2v2(const double * restrict A,
                       const double * restrict B,
                             double * restrict C) {
  CPLX
    Creg00,
    Creg10;
  CPLX Areg0,Areg1;
  CPLX Breg0;


    /*  Computing C(0:1,0:1)  */

      Areg0 = __lfpd(&A[0]) ;
      Areg1 = __lfpd(&A[8]) ;

      Breg0 = __lfpd(&B[0]) ;

      Creg00 = __fxsmul(Breg0,__cimag(Areg0)) ;

      Creg10 = __fxsmul(Breg0,__cimag(Areg1)) ;


      Areg0 = __lfpd(&A[2]) ;
      Areg1 = __lfpd(&A[10]) ;

      Breg0 = __lfpd(&B[8]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;


      Breg0 = __lfpd(&B[16]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;


      Areg0 = __lfpd(&A[4]) ;
      Areg1 = __lfpd(&A[12]) ;

      Breg0 = __lfpd(&B[24]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;


      Breg0 = __lfpd(&B[32]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;


      Areg0 = __lfpd(&A[6]) ;
      Areg1 = __lfpd(&A[14]) ;

      Breg0 = __lfpd(&B[40]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;


      Breg0 = __lfpd(&B[48]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;


      __stfpd(&C[0],Creg00) ;
      __stfpd(&C[8],Creg10) ;


    /*  Computing C(0:1,2:3)  */

      Areg0 = __lfpd(&A[0]) ;
      Areg1 = __lfpd(&A[8]) ;

      Breg0 = __lfpd(&B[2]) ;

      Creg00 = __fxsmul(Breg0,__cimag(Areg0)) ;

      Creg10 = __fxsmul(Breg0,__cimag(Areg1)) ;


      Areg0 = __lfpd(&A[2]) ;
      Areg1 = __lfpd(&A[10]) ;

      Breg0 = __lfpd(&B[10]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;


      Breg0 = __lfpd(&B[18]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;


      Areg0 = __lfpd(&A[4]) ;
      Areg1 = __lfpd(&A[12]) ;

      Breg0 = __lfpd(&B[26]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;


      Breg0 = __lfpd(&B[34]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;


      Areg0 = __lfpd(&A[6]) ;
      Areg1 = __lfpd(&A[14]) ;

      Breg0 = __lfpd(&B[42]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;


      Breg0 = __lfpd(&B[50]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;


      __stfpd(&C[2],Creg00) ;
      __stfpd(&C[10],Creg10) ;


    /*  Computing C(0:1,4:5)  */

      Areg0 = __lfpd(&A[0]) ;
      Areg1 = __lfpd(&A[8]) ;

      Breg0 = __lfpd(&B[4]) ;

      Creg00 = __fxsmul(Breg0,__cimag(Areg0)) ;

      Creg10 = __fxsmul(Breg0,__cimag(Areg1)) ;


      Areg0 = __lfpd(&A[2]) ;
      Areg1 = __lfpd(&A[10]) ;

      Breg0 = __lfpd(&B[12]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;


      Breg0 = __lfpd(&B[20]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;


      Areg0 = __lfpd(&A[4]) ;
      Areg1 = __lfpd(&A[12]) ;

      Breg0 = __lfpd(&B[28]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;


      Breg0 = __lfpd(&B[36]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;


      Areg0 = __lfpd(&A[6]) ;
      Areg1 = __lfpd(&A[14]) ;

      Breg0 = __lfpd(&B[44]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;


      Breg0 = __lfpd(&B[52]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;


      __stfpd(&C[4],Creg00) ;
      __stfpd(&C[12],Creg10) ;


    /*  Computing C(0:1,6:7)  */

      Areg0 = __lfpd(&A[0]) ;
      Areg1 = __lfpd(&A[8]) ;

      Breg0 = __lfpd(&B[6]) ;

      Creg00 = __fxsmul(Breg0,__cimag(Areg0)) ;

      Creg10 = __fxsmul(Breg0,__cimag(Areg1)) ;


      Areg0 = __lfpd(&A[2]) ;
      Areg1 = __lfpd(&A[10]) ;

      Breg0 = __lfpd(&B[14]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;


      Breg0 = __lfpd(&B[22]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;


      Areg0 = __lfpd(&A[4]) ;
      Areg1 = __lfpd(&A[12]) ;

      Breg0 = __lfpd(&B[30]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;


      Breg0 = __lfpd(&B[38]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;


      Areg0 = __lfpd(&A[6]) ;
      Areg1 = __lfpd(&A[14]) ;

      Breg0 = __lfpd(&B[46]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;


      Breg0 = __lfpd(&B[54]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;


      __stfpd(&C[6],Creg00) ;
      __stfpd(&C[14],Creg10) ;


    /*  Computing C(2:3,0:1)  */

      Areg0 = __lfpd(&A[16]) ;
      Areg1 = __lfpd(&A[24]) ;

      Breg0 = __lfpd(&B[0]) ;

      Creg00 = __fxsmul(Breg0,__cimag(Areg0)) ;

      Creg10 = __fxsmul(Breg0,__cimag(Areg1)) ;


      Areg0 = __lfpd(&A[18]) ;
      Areg1 = __lfpd(&A[26]) ;

      Breg0 = __lfpd(&B[8]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;


      Breg0 = __lfpd(&B[16]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;


      Areg0 = __lfpd(&A[20]) ;
      Areg1 = __lfpd(&A[28]) ;

      Breg0 = __lfpd(&B[24]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;


      Breg0 = __lfpd(&B[32]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;


      Areg0 = __lfpd(&A[22]) ;
      Areg1 = __lfpd(&A[30]) ;

      Breg0 = __lfpd(&B[40]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;


      Breg0 = __lfpd(&B[48]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;


      __stfpd(&C[16],Creg00) ;
      __stfpd(&C[24],Creg10) ;


    /*  Computing C(2:3,2:3)  */

      Areg0 = __lfpd(&A[16]) ;
      Areg1 = __lfpd(&A[24]) ;

      Breg0 = __lfpd(&B[2]) ;

      Creg00 = __fxsmul(Breg0,__cimag(Areg0)) ;

      Creg10 = __fxsmul(Breg0,__cimag(Areg1)) ;


      Areg0 = __lfpd(&A[18]) ;
      Areg1 = __lfpd(&A[26]) ;

      Breg0 = __lfpd(&B[10]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;


      Breg0 = __lfpd(&B[18]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;


      Areg0 = __lfpd(&A[20]) ;
      Areg1 = __lfpd(&A[28]) ;

      Breg0 = __lfpd(&B[26]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;


      Breg0 = __lfpd(&B[34]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;


      Areg0 = __lfpd(&A[22]) ;
      Areg1 = __lfpd(&A[30]) ;

      Breg0 = __lfpd(&B[42]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;


      Breg0 = __lfpd(&B[50]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;


      __stfpd(&C[18],Creg00) ;
      __stfpd(&C[26],Creg10) ;


    /*  Computing C(2:3,4:5)  */

      Areg0 = __lfpd(&A[16]) ;
      Areg1 = __lfpd(&A[24]) ;

      Breg0 = __lfpd(&B[4]) ;

      Creg00 = __fxsmul(Breg0,__cimag(Areg0)) ;

      Creg10 = __fxsmul(Breg0,__cimag(Areg1)) ;


      Areg0 = __lfpd(&A[18]) ;
      Areg1 = __lfpd(&A[26]) ;

      Breg0 = __lfpd(&B[12]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;


      Breg0 = __lfpd(&B[20]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;


      Areg0 = __lfpd(&A[20]) ;
      Areg1 = __lfpd(&A[28]) ;

      Breg0 = __lfpd(&B[28]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;


      Breg0 = __lfpd(&B[36]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;


      Areg0 = __lfpd(&A[22]) ;
      Areg1 = __lfpd(&A[30]) ;

      Breg0 = __lfpd(&B[44]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;


      Breg0 = __lfpd(&B[52]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;


      __stfpd(&C[20],Creg00) ;
      __stfpd(&C[28],Creg10) ;


    /*  Computing C(2:3,6:7)  */

      Areg0 = __lfpd(&A[16]) ;
      Areg1 = __lfpd(&A[24]) ;

      Breg0 = __lfpd(&B[6]) ;

      Creg00 = __fxsmul(Breg0,__cimag(Areg0)) ;

      Creg10 = __fxsmul(Breg0,__cimag(Areg1)) ;


      Areg0 = __lfpd(&A[18]) ;
      Areg1 = __lfpd(&A[26]) ;

      Breg0 = __lfpd(&B[14]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;


      Breg0 = __lfpd(&B[22]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;


      Areg0 = __lfpd(&A[20]) ;
      Areg1 = __lfpd(&A[28]) ;

      Breg0 = __lfpd(&B[30]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;


      Breg0 = __lfpd(&B[38]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;


      Areg0 = __lfpd(&A[22]) ;
      Areg1 = __lfpd(&A[30]) ;

      Breg0 = __lfpd(&B[46]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;


      Breg0 = __lfpd(&B[54]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;


      __stfpd(&C[22],Creg00) ;
      __stfpd(&C[30],Creg10) ;


    /*  Computing C(4:5,0:1)  */

      Areg0 = __lfpd(&A[32]) ;
      Areg1 = __lfpd(&A[40]) ;

      Breg0 = __lfpd(&B[0]) ;

      Creg00 = __fxsmul(Breg0,__cimag(Areg0)) ;

      Creg10 = __fxsmul(Breg0,__cimag(Areg1)) ;


      Areg0 = __lfpd(&A[34]) ;
      Areg1 = __lfpd(&A[42]) ;

      Breg0 = __lfpd(&B[8]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;


      Breg0 = __lfpd(&B[16]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;


      Areg0 = __lfpd(&A[36]) ;
      Areg1 = __lfpd(&A[44]) ;

      Breg0 = __lfpd(&B[24]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;


      Breg0 = __lfpd(&B[32]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;


      Areg0 = __lfpd(&A[38]) ;
      Areg1 = __lfpd(&A[46]) ;

      Breg0 = __lfpd(&B[40]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;


      Breg0 = __lfpd(&B[48]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;


      __stfpd(&C[32],Creg00) ;
      __stfpd(&C[40],Creg10) ;


    /*  Computing C(4:5,2:3)  */

      Areg0 = __lfpd(&A[32]) ;
      Areg1 = __lfpd(&A[40]) ;

      Breg0 = __lfpd(&B[2]) ;

      Creg00 = __fxsmul(Breg0,__cimag(Areg0)) ;

      Creg10 = __fxsmul(Breg0,__cimag(Areg1)) ;


      Areg0 = __lfpd(&A[34]) ;
      Areg1 = __lfpd(&A[42]) ;

      Breg0 = __lfpd(&B[10]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;


      Breg0 = __lfpd(&B[18]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;


      Areg0 = __lfpd(&A[36]) ;
      Areg1 = __lfpd(&A[44]) ;

      Breg0 = __lfpd(&B[26]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;


      Breg0 = __lfpd(&B[34]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;


      Areg0 = __lfpd(&A[38]) ;
      Areg1 = __lfpd(&A[46]) ;

      Breg0 = __lfpd(&B[42]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;


      Breg0 = __lfpd(&B[50]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;


      __stfpd(&C[34],Creg00) ;
      __stfpd(&C[42],Creg10) ;


    /*  Computing C(4:5,4:5)  */

      Areg0 = __lfpd(&A[32]) ;
      Areg1 = __lfpd(&A[40]) ;

      Breg0 = __lfpd(&B[4]) ;

      Creg00 = __fxsmul(Breg0,__cimag(Areg0)) ;

      Creg10 = __fxsmul(Breg0,__cimag(Areg1)) ;


      Areg0 = __lfpd(&A[34]) ;
      Areg1 = __lfpd(&A[42]) ;

      Breg0 = __lfpd(&B[12]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;


      Breg0 = __lfpd(&B[20]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;


      Areg0 = __lfpd(&A[36]) ;
      Areg1 = __lfpd(&A[44]) ;

      Breg0 = __lfpd(&B[28]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;


      Breg0 = __lfpd(&B[36]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;


      Areg0 = __lfpd(&A[38]) ;
      Areg1 = __lfpd(&A[46]) ;

      Breg0 = __lfpd(&B[44]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;


      Breg0 = __lfpd(&B[52]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;


      __stfpd(&C[36],Creg00) ;
      __stfpd(&C[44],Creg10) ;


    /*  Computing C(4:5,6:7)  */

      Areg0 = __lfpd(&A[32]) ;
      Areg1 = __lfpd(&A[40]) ;

      Breg0 = __lfpd(&B[6]) ;

      Creg00 = __fxsmul(Breg0,__cimag(Areg0)) ;

      Creg10 = __fxsmul(Breg0,__cimag(Areg1)) ;


      Areg0 = __lfpd(&A[34]) ;
      Areg1 = __lfpd(&A[42]) ;

      Breg0 = __lfpd(&B[14]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;


      Breg0 = __lfpd(&B[22]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;


      Areg0 = __lfpd(&A[36]) ;
      Areg1 = __lfpd(&A[44]) ;

      Breg0 = __lfpd(&B[30]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;


      Breg0 = __lfpd(&B[38]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;


      Areg0 = __lfpd(&A[38]) ;
      Areg1 = __lfpd(&A[46]) ;

      Breg0 = __lfpd(&B[46]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;


      Breg0 = __lfpd(&B[54]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;


      __stfpd(&C[38],Creg00) ;
      __stfpd(&C[46],Creg10) ;


    /*  Computing C(6:6,0:1)  */

      Areg0 = __lfpd(&A[48]) ;

      Breg0 = __lfpd(&B[0]) ;

      Creg00 = __fxsmul(Breg0,__cimag(Areg0)) ;


      Areg0 = __lfpd(&A[50]) ;

      Breg0 = __lfpd(&B[8]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;


      Breg0 = __lfpd(&B[16]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;


      Areg0 = __lfpd(&A[52]) ;

      Breg0 = __lfpd(&B[24]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;


      Breg0 = __lfpd(&B[32]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;


      Areg0 = __lfpd(&A[54]) ;

      Breg0 = __lfpd(&B[40]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;


      Breg0 = __lfpd(&B[48]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;


      __stfpd(&C[48],Creg00) ;


    /*  Computing C(6:6,2:3)  */

      Areg0 = __lfpd(&A[48]) ;

      Breg0 = __lfpd(&B[2]) ;

      Creg00 = __fxsmul(Breg0,__cimag(Areg0)) ;


      Areg0 = __lfpd(&A[50]) ;

      Breg0 = __lfpd(&B[10]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;


      Breg0 = __lfpd(&B[18]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;


      Areg0 = __lfpd(&A[52]) ;

      Breg0 = __lfpd(&B[26]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;


      Breg0 = __lfpd(&B[34]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;


      Areg0 = __lfpd(&A[54]) ;

      Breg0 = __lfpd(&B[42]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;


      Breg0 = __lfpd(&B[50]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;


      __stfpd(&C[50],Creg00) ;


    /*  Computing C(6:6,4:5)  */

      Areg0 = __lfpd(&A[48]) ;

      Breg0 = __lfpd(&B[4]) ;

      Creg00 = __fxsmul(Breg0,__cimag(Areg0)) ;


      Areg0 = __lfpd(&A[50]) ;

      Breg0 = __lfpd(&B[12]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;


      Breg0 = __lfpd(&B[20]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;


      Areg0 = __lfpd(&A[52]) ;

      Breg0 = __lfpd(&B[28]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;


      Breg0 = __lfpd(&B[36]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;


      Areg0 = __lfpd(&A[54]) ;

      Breg0 = __lfpd(&B[44]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;


      Breg0 = __lfpd(&B[52]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;


      __stfpd(&C[52],Creg00) ;


    /*  Computing C(6:6,6:7)  */

      Areg0 = __lfpd(&A[48]) ;

      Breg0 = __lfpd(&B[6]) ;

      Creg00 = __fxsmul(Breg0,__cimag(Areg0)) ;


      Areg0 = __lfpd(&A[50]) ;

      Breg0 = __lfpd(&B[14]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;


      Breg0 = __lfpd(&B[22]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;


      Areg0 = __lfpd(&A[52]) ;

      Breg0 = __lfpd(&B[30]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;


      Breg0 = __lfpd(&B[38]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;


      Areg0 = __lfpd(&A[54]) ;

      Breg0 = __lfpd(&B[46]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;


      Breg0 = __lfpd(&B[54]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;


      __stfpd(&C[54],Creg00) ;


}
