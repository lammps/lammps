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
#define const
#define vec_fma vec_madd

void mmul_bgq_n5_lda8_2x8(const double * restrict A,
                          const double * restrict B,
                                double * restrict C) {
  vector
    Creg00,Creg04,
    Creg10,Creg14;
  vector Areg0,Areg1;
  vector Breg0,Breg4;
  vector Atmp;


    /*  Computing C(0:1,0:5)  */

      Areg0 = vec_ld(8*0,A) ;
      Areg1 = vec_ld(8*8,A) ;

      Breg0 = vec_ld(8*0,B) ;
      Breg4 = vec_ld(8*4,B) ;

      Atmp = vec_splat(Areg0,1) ;
      Creg00 = vec_mul(Atmp,Breg0) ;
      Creg04 = vec_mul(Atmp,Breg4) ;

      Atmp = vec_splat(Areg1,1) ;
      Creg10 = vec_mul(Atmp,Breg0) ;
      Creg14 = vec_mul(Atmp,Breg4) ;


      Breg0 = vec_ld(8*8,B) ;
      Breg4 = vec_ld(8*12,B) ;

      Atmp = vec_splat(Areg0,2) ;
      Creg00 = vec_fma(Atmp,Breg0,Creg00) ;
      Creg04 = vec_fma(Atmp,Breg4,Creg04) ;

      Atmp = vec_splat(Areg1,2) ;
      Creg10 = vec_fma(Atmp,Breg0,Creg10) ;
      Creg14 = vec_fma(Atmp,Breg4,Creg14) ;


      Breg0 = vec_ld(8*16,B) ;
      Breg4 = vec_ld(8*20,B) ;

      Atmp = vec_splat(Areg0,3) ;
      Creg00 = vec_fma(Atmp,Breg0,Creg00) ;
      Creg04 = vec_fma(Atmp,Breg4,Creg04) ;

      Atmp = vec_splat(Areg1,3) ;
      Creg10 = vec_fma(Atmp,Breg0,Creg10) ;
      Creg14 = vec_fma(Atmp,Breg4,Creg14) ;


      Areg0 = vec_ld(8*4,A) ;
      Areg1 = vec_ld(8*12,A) ;

      Breg0 = vec_ld(8*24,B) ;
      Breg4 = vec_ld(8*28,B) ;

      Atmp = vec_splat(Areg0,0) ;
      Creg00 = vec_fma(Atmp,Breg0,Creg00) ;
      Creg04 = vec_fma(Atmp,Breg4,Creg04) ;

      Atmp = vec_splat(Areg1,0) ;
      Creg10 = vec_fma(Atmp,Breg0,Creg10) ;
      Creg14 = vec_fma(Atmp,Breg4,Creg14) ;


      Breg0 = vec_ld(8*32,B) ;
      Breg4 = vec_ld(8*36,B) ;

      Atmp = vec_splat(Areg0,1) ;
      Creg00 = vec_fma(Atmp,Breg0,Creg00) ;
      Creg04 = vec_fma(Atmp,Breg4,Creg04) ;

      Atmp = vec_splat(Areg1,1) ;
      Creg10 = vec_fma(Atmp,Breg0,Creg10) ;
      Creg14 = vec_fma(Atmp,Breg4,Creg14) ;


      vec_st(Creg00,8*0,C) ;
      vec_st(Creg04,8*4,C) ;
      vec_st(Creg10,8*8,C) ;
      vec_st(Creg14,8*12,C) ;


    /*  Computing C(2:3,0:5)  */

      Areg0 = vec_ld(8*16,A) ;
      Areg1 = vec_ld(8*24,A) ;

      Breg0 = vec_ld(8*0,B) ;
      Breg4 = vec_ld(8*4,B) ;

      Atmp = vec_splat(Areg0,1) ;
      Creg00 = vec_mul(Atmp,Breg0) ;
      Creg04 = vec_mul(Atmp,Breg4) ;

      Atmp = vec_splat(Areg1,1) ;
      Creg10 = vec_mul(Atmp,Breg0) ;
      Creg14 = vec_mul(Atmp,Breg4) ;


      Breg0 = vec_ld(8*8,B) ;
      Breg4 = vec_ld(8*12,B) ;

      Atmp = vec_splat(Areg0,2) ;
      Creg00 = vec_fma(Atmp,Breg0,Creg00) ;
      Creg04 = vec_fma(Atmp,Breg4,Creg04) ;

      Atmp = vec_splat(Areg1,2) ;
      Creg10 = vec_fma(Atmp,Breg0,Creg10) ;
      Creg14 = vec_fma(Atmp,Breg4,Creg14) ;


      Breg0 = vec_ld(8*16,B) ;
      Breg4 = vec_ld(8*20,B) ;

      Atmp = vec_splat(Areg0,3) ;
      Creg00 = vec_fma(Atmp,Breg0,Creg00) ;
      Creg04 = vec_fma(Atmp,Breg4,Creg04) ;

      Atmp = vec_splat(Areg1,3) ;
      Creg10 = vec_fma(Atmp,Breg0,Creg10) ;
      Creg14 = vec_fma(Atmp,Breg4,Creg14) ;


      Areg0 = vec_ld(8*20,A) ;
      Areg1 = vec_ld(8*28,A) ;

      Breg0 = vec_ld(8*24,B) ;
      Breg4 = vec_ld(8*28,B) ;

      Atmp = vec_splat(Areg0,0) ;
      Creg00 = vec_fma(Atmp,Breg0,Creg00) ;
      Creg04 = vec_fma(Atmp,Breg4,Creg04) ;

      Atmp = vec_splat(Areg1,0) ;
      Creg10 = vec_fma(Atmp,Breg0,Creg10) ;
      Creg14 = vec_fma(Atmp,Breg4,Creg14) ;


      Breg0 = vec_ld(8*32,B) ;
      Breg4 = vec_ld(8*36,B) ;

      Atmp = vec_splat(Areg0,1) ;
      Creg00 = vec_fma(Atmp,Breg0,Creg00) ;
      Creg04 = vec_fma(Atmp,Breg4,Creg04) ;

      Atmp = vec_splat(Areg1,1) ;
      Creg10 = vec_fma(Atmp,Breg0,Creg10) ;
      Creg14 = vec_fma(Atmp,Breg4,Creg14) ;


      vec_st(Creg00,8*16,C) ;
      vec_st(Creg04,8*20,C) ;
      vec_st(Creg10,8*24,C) ;
      vec_st(Creg14,8*28,C) ;


    /*  Computing C(4:4,0:5)  */

      Areg0 = vec_ld(8*32,A) ;

      Breg0 = vec_ld(8*0,B) ;
      Breg4 = vec_ld(8*4,B) ;

      Atmp = vec_splat(Areg0,1) ;
      Creg00 = vec_mul(Atmp,Breg0) ;
      Creg04 = vec_mul(Atmp,Breg4) ;


      Breg0 = vec_ld(8*8,B) ;
      Breg4 = vec_ld(8*12,B) ;

      Atmp = vec_splat(Areg0,2) ;
      Creg00 = vec_fma(Atmp,Breg0,Creg00) ;
      Creg04 = vec_fma(Atmp,Breg4,Creg04) ;


      Breg0 = vec_ld(8*16,B) ;
      Breg4 = vec_ld(8*20,B) ;

      Atmp = vec_splat(Areg0,3) ;
      Creg00 = vec_fma(Atmp,Breg0,Creg00) ;
      Creg04 = vec_fma(Atmp,Breg4,Creg04) ;


      Areg0 = vec_ld(8*36,A) ;

      Breg0 = vec_ld(8*24,B) ;
      Breg4 = vec_ld(8*28,B) ;

      Atmp = vec_splat(Areg0,0) ;
      Creg00 = vec_fma(Atmp,Breg0,Creg00) ;
      Creg04 = vec_fma(Atmp,Breg4,Creg04) ;


      Breg0 = vec_ld(8*32,B) ;
      Breg4 = vec_ld(8*36,B) ;

      Atmp = vec_splat(Areg0,1) ;
      Creg00 = vec_fma(Atmp,Breg0,Creg00) ;
      Creg04 = vec_fma(Atmp,Breg4,Creg04) ;


      vec_st(Creg00,8*32,C) ;
      vec_st(Creg04,8*36,C) ;


}

#undef vector
#undef const
#undef vec_fma
