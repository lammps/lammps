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

void mmul_bg_5_8_5x2v2(const double * restrict A,
                       const double * restrict B,
                             double * restrict C) {
  CPLX
    Creg00,
    Creg10,
    Creg20,
    Creg30,
    Creg40;
  CPLX Areg0,Areg1,Areg2,Areg3,Areg4;
  CPLX Breg0;


    /*  Computing C(0:4,0:1)  */

      Areg0 = __lfpd(&A[0]) ;
      Areg1 = __lfpd(&A[8]) ;
      Areg2 = __lfpd(&A[16]) ;
      Areg3 = __lfpd(&A[24]) ;
      Areg4 = __lfpd(&A[32]) ;

      Breg0 = __lfpd(&B[0]) ;

      Creg00 = __fxsmul(Breg0,__cimag(Areg0)) ;

      Creg10 = __fxsmul(Breg0,__cimag(Areg1)) ;

      Creg20 = __fxsmul(Breg0,__cimag(Areg2)) ;

      Creg30 = __fxsmul(Breg0,__cimag(Areg3)) ;

      Creg40 = __fxsmul(Breg0,__cimag(Areg4)) ;


      Areg0 = __lfpd(&A[2]) ;
      Areg1 = __lfpd(&A[10]) ;
      Areg2 = __lfpd(&A[18]) ;
      Areg3 = __lfpd(&A[26]) ;
      Areg4 = __lfpd(&A[34]) ;

      Breg0 = __lfpd(&B[8]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;

      Creg20 = __fxcpmadd(Creg20,Breg0,__creal(Areg2)) ;

      Creg30 = __fxcpmadd(Creg30,Breg0,__creal(Areg3)) ;

      Creg40 = __fxcpmadd(Creg40,Breg0,__creal(Areg4)) ;


      Breg0 = __lfpd(&B[16]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;

      Creg20 = __fxcpmadd(Creg20,Breg0,__cimag(Areg2)) ;

      Creg30 = __fxcpmadd(Creg30,Breg0,__cimag(Areg3)) ;

      Creg40 = __fxcpmadd(Creg40,Breg0,__cimag(Areg4)) ;


      Areg0 = __lfpd(&A[4]) ;
      Areg1 = __lfpd(&A[12]) ;
      Areg2 = __lfpd(&A[20]) ;
      Areg3 = __lfpd(&A[28]) ;
      Areg4 = __lfpd(&A[36]) ;

      Breg0 = __lfpd(&B[24]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;

      Creg20 = __fxcpmadd(Creg20,Breg0,__creal(Areg2)) ;

      Creg30 = __fxcpmadd(Creg30,Breg0,__creal(Areg3)) ;

      Creg40 = __fxcpmadd(Creg40,Breg0,__creal(Areg4)) ;


      Breg0 = __lfpd(&B[32]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;

      Creg20 = __fxcpmadd(Creg20,Breg0,__cimag(Areg2)) ;

      Creg30 = __fxcpmadd(Creg30,Breg0,__cimag(Areg3)) ;

      Creg40 = __fxcpmadd(Creg40,Breg0,__cimag(Areg4)) ;


      __stfpd(&C[0],Creg00) ;
      __stfpd(&C[8],Creg10) ;
      __stfpd(&C[16],Creg20) ;
      __stfpd(&C[24],Creg30) ;
      __stfpd(&C[32],Creg40) ;


    /*  Computing C(0:4,2:3)  */

      Areg0 = __lfpd(&A[0]) ;
      Areg1 = __lfpd(&A[8]) ;
      Areg2 = __lfpd(&A[16]) ;
      Areg3 = __lfpd(&A[24]) ;
      Areg4 = __lfpd(&A[32]) ;

      Breg0 = __lfpd(&B[2]) ;

      Creg00 = __fxsmul(Breg0,__cimag(Areg0)) ;

      Creg10 = __fxsmul(Breg0,__cimag(Areg1)) ;

      Creg20 = __fxsmul(Breg0,__cimag(Areg2)) ;

      Creg30 = __fxsmul(Breg0,__cimag(Areg3)) ;

      Creg40 = __fxsmul(Breg0,__cimag(Areg4)) ;


      Areg0 = __lfpd(&A[2]) ;
      Areg1 = __lfpd(&A[10]) ;
      Areg2 = __lfpd(&A[18]) ;
      Areg3 = __lfpd(&A[26]) ;
      Areg4 = __lfpd(&A[34]) ;

      Breg0 = __lfpd(&B[10]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;

      Creg20 = __fxcpmadd(Creg20,Breg0,__creal(Areg2)) ;

      Creg30 = __fxcpmadd(Creg30,Breg0,__creal(Areg3)) ;

      Creg40 = __fxcpmadd(Creg40,Breg0,__creal(Areg4)) ;


      Breg0 = __lfpd(&B[18]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;

      Creg20 = __fxcpmadd(Creg20,Breg0,__cimag(Areg2)) ;

      Creg30 = __fxcpmadd(Creg30,Breg0,__cimag(Areg3)) ;

      Creg40 = __fxcpmadd(Creg40,Breg0,__cimag(Areg4)) ;


      Areg0 = __lfpd(&A[4]) ;
      Areg1 = __lfpd(&A[12]) ;
      Areg2 = __lfpd(&A[20]) ;
      Areg3 = __lfpd(&A[28]) ;
      Areg4 = __lfpd(&A[36]) ;

      Breg0 = __lfpd(&B[26]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;

      Creg20 = __fxcpmadd(Creg20,Breg0,__creal(Areg2)) ;

      Creg30 = __fxcpmadd(Creg30,Breg0,__creal(Areg3)) ;

      Creg40 = __fxcpmadd(Creg40,Breg0,__creal(Areg4)) ;


      Breg0 = __lfpd(&B[34]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;

      Creg20 = __fxcpmadd(Creg20,Breg0,__cimag(Areg2)) ;

      Creg30 = __fxcpmadd(Creg30,Breg0,__cimag(Areg3)) ;

      Creg40 = __fxcpmadd(Creg40,Breg0,__cimag(Areg4)) ;


      __stfpd(&C[2],Creg00) ;
      __stfpd(&C[10],Creg10) ;
      __stfpd(&C[18],Creg20) ;
      __stfpd(&C[26],Creg30) ;
      __stfpd(&C[34],Creg40) ;


    /*  Computing C(0:4,4:5)  */

      Areg0 = __lfpd(&A[0]) ;
      Areg1 = __lfpd(&A[8]) ;
      Areg2 = __lfpd(&A[16]) ;
      Areg3 = __lfpd(&A[24]) ;
      Areg4 = __lfpd(&A[32]) ;

      Breg0 = __lfpd(&B[4]) ;

      Creg00 = __fxsmul(Breg0,__cimag(Areg0)) ;

      Creg10 = __fxsmul(Breg0,__cimag(Areg1)) ;

      Creg20 = __fxsmul(Breg0,__cimag(Areg2)) ;

      Creg30 = __fxsmul(Breg0,__cimag(Areg3)) ;

      Creg40 = __fxsmul(Breg0,__cimag(Areg4)) ;


      Areg0 = __lfpd(&A[2]) ;
      Areg1 = __lfpd(&A[10]) ;
      Areg2 = __lfpd(&A[18]) ;
      Areg3 = __lfpd(&A[26]) ;
      Areg4 = __lfpd(&A[34]) ;

      Breg0 = __lfpd(&B[12]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;

      Creg20 = __fxcpmadd(Creg20,Breg0,__creal(Areg2)) ;

      Creg30 = __fxcpmadd(Creg30,Breg0,__creal(Areg3)) ;

      Creg40 = __fxcpmadd(Creg40,Breg0,__creal(Areg4)) ;


      Breg0 = __lfpd(&B[20]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;

      Creg20 = __fxcpmadd(Creg20,Breg0,__cimag(Areg2)) ;

      Creg30 = __fxcpmadd(Creg30,Breg0,__cimag(Areg3)) ;

      Creg40 = __fxcpmadd(Creg40,Breg0,__cimag(Areg4)) ;


      Areg0 = __lfpd(&A[4]) ;
      Areg1 = __lfpd(&A[12]) ;
      Areg2 = __lfpd(&A[20]) ;
      Areg3 = __lfpd(&A[28]) ;
      Areg4 = __lfpd(&A[36]) ;

      Breg0 = __lfpd(&B[28]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__creal(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__creal(Areg1)) ;

      Creg20 = __fxcpmadd(Creg20,Breg0,__creal(Areg2)) ;

      Creg30 = __fxcpmadd(Creg30,Breg0,__creal(Areg3)) ;

      Creg40 = __fxcpmadd(Creg40,Breg0,__creal(Areg4)) ;


      Breg0 = __lfpd(&B[36]) ;

      Creg00 = __fxcpmadd(Creg00,Breg0,__cimag(Areg0)) ;

      Creg10 = __fxcpmadd(Creg10,Breg0,__cimag(Areg1)) ;

      Creg20 = __fxcpmadd(Creg20,Breg0,__cimag(Areg2)) ;

      Creg30 = __fxcpmadd(Creg30,Breg0,__cimag(Areg3)) ;

      Creg40 = __fxcpmadd(Creg40,Breg0,__cimag(Areg4)) ;


      __stfpd(&C[4],Creg00) ;
      __stfpd(&C[12],Creg10) ;
      __stfpd(&C[20],Creg20) ;
      __stfpd(&C[28],Creg30) ;
      __stfpd(&C[36],Creg40) ;


}
