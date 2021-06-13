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

#define __creal    my_creal
#define __cimag    my_cimag
#define __lfpd     my_lfpd
#define __stfpd    my_stfpd
#define __fxsmul   my_fxsmul
#define __fxcpmadd my_fxcpmadd

#else

#define CPLX double _Complex

#endif

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

void bgmul_7(double (* restrict Ain)[8],double (* restrict Bin)[8],double (* restrict Cin)[8]) {

CPLX C1,C2,C3,C4,C5,C6,C7;
CPLX A1_01,A1_23,A1_45,A1_67;
CPLX       A2_23,A2_45,A2_67;
CPLX       A3_23,A3_45,A3_67;
CPLX             A4_45,A4_67;
CPLX             A5_45,A5_67;
CPLX                   A6_67;
CPLX                   A7_67;
CPLX Bkj;

 double (* restrict A)[8] = &Ain[-1];
 double (* restrict B)[8] = &Bin[-1];
 double (* restrict C)[8] = &Cin[-1];


int j;
A1_01 = __lfpd(&A[1][0]); A1_23 = __lfpd(&A[1][2]); A1_45 = __lfpd(&A[1][4]); A1_67 = __lfpd(&A[1][6]);
                          A2_23 = __lfpd(&A[2][2]); A2_45 = __lfpd(&A[2][4]); A2_67 = __lfpd(&A[2][6]);
                          A3_23 = __lfpd(&A[3][2]); A3_45 = __lfpd(&A[3][4]); A3_67 = __lfpd(&A[3][6]);
                                                    A4_45 = __lfpd(&A[4][4]); A4_67 = __lfpd(&A[4][6]);
                                                    A5_45 = __lfpd(&A[5][4]); A5_67 = __lfpd(&A[5][6]);
                                                                              A6_67 = __lfpd(&A[6][6]);
                                                                              A7_67 = __lfpd(&A[7][6]);

for(j = 0; j<7; j+=2) {
  /* k = 1 */
  Bkj = __lfpd(&B[1][j]);

    C1 = __fxsmul(Bkj,__cimag(A1_01));
    C2 = __fxsmul(Bkj,__creal(A1_23));
    C3 = __fxsmul(Bkj,__cimag(A1_23));
    C4 = __fxsmul(Bkj,__creal(A1_45));
    C5 = __fxsmul(Bkj,__cimag(A1_45));
    C6 = __fxsmul(Bkj,__creal(A1_67));
    C7 = __fxsmul(Bkj,__cimag(A1_67));

  /* k = 2 */
  Bkj = __lfpd(&B[2][j]);

    C1 = __fxcpmadd(C1,Bkj,__creal(A1_23));
    C2 = __fxcpmadd(C2,Bkj,__creal(A2_23));
    C3 = __fxcpmadd(C3,Bkj,__cimag(A2_23));
    C4 = __fxcpmadd(C4,Bkj,__creal(A2_45));
    C5 = __fxcpmadd(C5,Bkj,__cimag(A2_45));
    C6 = __fxcpmadd(C6,Bkj,__creal(A2_67));
    C7 = __fxcpmadd(C7,Bkj,__cimag(A2_67));

  /* k = 3 */
  Bkj = __lfpd(&B[3][j]);

    C1 = __fxcpmadd(C1,Bkj,__cimag(A1_23));
    C2 = __fxcpmadd(C2,Bkj,__cimag(A2_23));
    C3 = __fxcpmadd(C3,Bkj,__cimag(A3_23));
    C4 = __fxcpmadd(C4,Bkj,__creal(A3_45));
    C5 = __fxcpmadd(C5,Bkj,__cimag(A3_45));
    C6 = __fxcpmadd(C6,Bkj,__creal(A3_67));
    C7 = __fxcpmadd(C7,Bkj,__cimag(A3_67));

  /* k = 4 */
  Bkj = __lfpd(&B[4][j]);

    C1 = __fxcpmadd(C1,Bkj,__creal(A1_45));
    C2 = __fxcpmadd(C2,Bkj,__creal(A2_45));
    C3 = __fxcpmadd(C3,Bkj,__creal(A3_45));
    C4 = __fxcpmadd(C4,Bkj,__creal(A4_45));
    C5 = __fxcpmadd(C5,Bkj,__cimag(A4_45));
    C6 = __fxcpmadd(C6,Bkj,__creal(A4_67));
    C7 = __fxcpmadd(C7,Bkj,__cimag(A4_67));

  /* k = 5 */
  Bkj = __lfpd(&B[5][j]);

    C1 = __fxcpmadd(C1,Bkj,__cimag(A1_45));
    C2 = __fxcpmadd(C2,Bkj,__cimag(A2_45));
    C3 = __fxcpmadd(C3,Bkj,__cimag(A3_45));
    C4 = __fxcpmadd(C4,Bkj,__cimag(A4_45));
    C5 = __fxcpmadd(C5,Bkj,__cimag(A5_45));
    C6 = __fxcpmadd(C6,Bkj,__creal(A5_67));
    C7 = __fxcpmadd(C7,Bkj,__cimag(A5_67));

  /* k = 6 */
  Bkj = __lfpd(&B[6][j]);

    C1 = __fxcpmadd(C1,Bkj,__creal(A1_67));
    C2 = __fxcpmadd(C2,Bkj,__creal(A2_67));
    C3 = __fxcpmadd(C3,Bkj,__creal(A3_67));
    C4 = __fxcpmadd(C4,Bkj,__creal(A4_67));
    C5 = __fxcpmadd(C5,Bkj,__creal(A5_67));
    C6 = __fxcpmadd(C6,Bkj,__creal(A6_67));
    C7 = __fxcpmadd(C7,Bkj,__cimag(A6_67));

  /* k = 7 */
  Bkj = __lfpd(&B[7][j]);

    C1 = __fxcpmadd(C1,Bkj,__cimag(A1_67));
    C2 = __fxcpmadd(C2,Bkj,__cimag(A2_67));
    C3 = __fxcpmadd(C3,Bkj,__cimag(A3_67));
    C4 = __fxcpmadd(C4,Bkj,__cimag(A4_67));
    C5 = __fxcpmadd(C5,Bkj,__cimag(A5_67));
    C6 = __fxcpmadd(C6,Bkj,__cimag(A6_67));
    C7 = __fxcpmadd(C7,Bkj,__cimag(A7_67));

  __stfpd(&C[1][j],C1);
  __stfpd(&C[2][j],C2);
  __stfpd(&C[3][j],C3);
  __stfpd(&C[4][j],C4);
  __stfpd(&C[5][j],C5);
  __stfpd(&C[6][j],C6);
  __stfpd(&C[7][j],C7);
}
}
