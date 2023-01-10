// clang-format off
/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

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
#define my_vec_add(a,b)     vec_add(a,b)
#define my_vec_mul(a,b)     vec_mul(a,b)
#define my_vec_fma(a,b,c)   vec_madd(b,c,a)
#define my_vec_ld(ptr)      vec_lda(0,ptr)
#define my_vec_sldw(x,y,n)  vec_sldw(x,y,n)
#define my_vec_sts(x,ptr)   vec_sts(x,0,ptr)

#define const
#define real double

void ttr_bg_5_8_3_v4r1(const real * restrict A,
    const real * restrict B0,real * restrict tout0,
    const real * restrict B1,real * restrict tout1,
    const real * restrict B2,real * restrict tout2) {
vector Areg1;
vector B0reg1,B1reg1,B2reg1;
vector T0reg1,T1reg1,T2reg1;

Areg1 = my_vec_ld(&A[0]) ;
B0reg1 = my_vec_ld(&B0[0]) ;
B1reg1 = my_vec_ld(&B1[0]) ;
B2reg1 = my_vec_ld(&B2[0]) ;
T0reg1 = my_vec_mul(Areg1,B0reg1) ;
T1reg1 = my_vec_mul(Areg1,B1reg1) ;
T2reg1 = my_vec_mul(Areg1,B2reg1) ;

Areg1 = my_vec_ld(&A[4]) ;
B0reg1 = my_vec_ld(&B0[4]) ;
B1reg1 = my_vec_ld(&B1[4]) ;
B2reg1 = my_vec_ld(&B2[4]) ;
T0reg1 = my_vec_fma(T0reg1,Areg1,B0reg1) ;
T1reg1 = my_vec_fma(T1reg1,Areg1,B1reg1) ;
T2reg1 = my_vec_fma(T2reg1,Areg1,B2reg1) ;

Areg1 = my_vec_ld(&A[8]) ;
B0reg1 = my_vec_ld(&B0[8]) ;
B1reg1 = my_vec_ld(&B1[8]) ;
B2reg1 = my_vec_ld(&B2[8]) ;
T0reg1 = my_vec_fma(T0reg1,Areg1,B0reg1) ;
T1reg1 = my_vec_fma(T1reg1,Areg1,B1reg1) ;
T2reg1 = my_vec_fma(T2reg1,Areg1,B2reg1) ;

Areg1 = my_vec_ld(&A[12]) ;
B0reg1 = my_vec_ld(&B0[12]) ;
B1reg1 = my_vec_ld(&B1[12]) ;
B2reg1 = my_vec_ld(&B2[12]) ;
T0reg1 = my_vec_fma(T0reg1,Areg1,B0reg1) ;
T1reg1 = my_vec_fma(T1reg1,Areg1,B1reg1) ;
T2reg1 = my_vec_fma(T2reg1,Areg1,B2reg1) ;

Areg1 = my_vec_ld(&A[16]) ;
B0reg1 = my_vec_ld(&B0[16]) ;
B1reg1 = my_vec_ld(&B1[16]) ;
B2reg1 = my_vec_ld(&B2[16]) ;
T0reg1 = my_vec_fma(T0reg1,Areg1,B0reg1) ;
T1reg1 = my_vec_fma(T1reg1,Areg1,B1reg1) ;
T2reg1 = my_vec_fma(T2reg1,Areg1,B2reg1) ;

Areg1 = my_vec_ld(&A[20]) ;
B0reg1 = my_vec_ld(&B0[20]) ;
B1reg1 = my_vec_ld(&B1[20]) ;
B2reg1 = my_vec_ld(&B2[20]) ;
T0reg1 = my_vec_fma(T0reg1,Areg1,B0reg1) ;
T1reg1 = my_vec_fma(T1reg1,Areg1,B1reg1) ;
T2reg1 = my_vec_fma(T2reg1,Areg1,B2reg1) ;

Areg1 = my_vec_ld(&A[24]) ;
B0reg1 = my_vec_ld(&B0[24]) ;
B1reg1 = my_vec_ld(&B1[24]) ;
B2reg1 = my_vec_ld(&B2[24]) ;
T0reg1 = my_vec_fma(T0reg1,Areg1,B0reg1) ;
T1reg1 = my_vec_fma(T1reg1,Areg1,B1reg1) ;
T2reg1 = my_vec_fma(T2reg1,Areg1,B2reg1) ;

Areg1 = my_vec_ld(&A[28]) ;
B0reg1 = my_vec_ld(&B0[28]) ;
B1reg1 = my_vec_ld(&B1[28]) ;
B2reg1 = my_vec_ld(&B2[28]) ;
T0reg1 = my_vec_fma(T0reg1,Areg1,B0reg1) ;
T1reg1 = my_vec_fma(T1reg1,Areg1,B1reg1) ;
T2reg1 = my_vec_fma(T2reg1,Areg1,B2reg1) ;

Areg1 = my_vec_ld(&A[32]) ;
B0reg1 = my_vec_ld(&B0[32]) ;
B1reg1 = my_vec_ld(&B1[32]) ;
B2reg1 = my_vec_ld(&B2[32]) ;
T0reg1 = my_vec_fma(T0reg1,Areg1,B0reg1) ;
T1reg1 = my_vec_fma(T1reg1,Areg1,B1reg1) ;
T2reg1 = my_vec_fma(T2reg1,Areg1,B2reg1) ;

Areg1 = my_vec_ld(&A[36]) ;
B0reg1 = my_vec_ld(&B0[36]) ;
B1reg1 = my_vec_ld(&B1[36]) ;
B2reg1 = my_vec_ld(&B2[36]) ;
T0reg1 = my_vec_fma(T0reg1,Areg1,B0reg1) ;
T1reg1 = my_vec_fma(T1reg1,Areg1,B1reg1) ;
T2reg1 = my_vec_fma(T2reg1,Areg1,B2reg1) ;



T0reg1 = my_vec_add(T0reg1,my_vec_sldw(T0reg1,T0reg1,2));
T0reg1 = my_vec_add(T0reg1,my_vec_sldw(T0reg1,T0reg1,1));
my_vec_sts(T0reg1,tout0);
T1reg1 = my_vec_add(T1reg1,my_vec_sldw(T1reg1,T1reg1,2));
T1reg1 = my_vec_add(T1reg1,my_vec_sldw(T1reg1,T1reg1,1));
my_vec_sts(T1reg1,tout1);
T2reg1 = my_vec_add(T2reg1,my_vec_sldw(T2reg1,T2reg1,2));
T2reg1 = my_vec_add(T2reg1,my_vec_sldw(T2reg1,T2reg1,1));
my_vec_sts(T2reg1,tout2);

}


#undef vector
#undef my_vec_add
#undef my_vec_mul
#undef my_vec_fma
#undef my_vec_ld
#undef my_vec_sldw
#undef my_vec_sts

#undef const
#undef real
