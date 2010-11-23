/***************************************************************************
                              ucl_arg_kludge.h
                             -------------------
                               W. Michael Brown

  Allow multiple arguments to be added for a kernel call at a single time

 __________________________________________________________________________
    This file is part of the Geryon Unified Coprocessor Library (UCL)
 __________________________________________________________________________

    begin                : Sun Feb 7 2010
    copyright            : (C) 2010 by W. Michael Brown
    email                : brownw@ornl.gov
 ***************************************************************************/
 
/* -----------------------------------------------------------------------
   Copyright (2010) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the Simplified BSD License.
   ----------------------------------------------------------------------- */

  template <class t1, class t2>
  inline void add_args(t1 *a1, t2 *a2) {
    add_arg(a1); add_arg(a2);
  }

  template <class t1, class t2, class t3>
  inline void add_args(t1 *a1, t2 *a2, t3 *a3) {
    add_arg(a1); add_arg(a2); add_arg(a3);
  }

  template <class t1, class t2, class t3, class t4>
  inline void add_args(t1 *a1, t2 *a2, t3 *a3, t4 *a4) {
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4);
  }

  template <class t1, class t2, class t3, class t4, class t5>
  inline void add_args(t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5) {
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6>
  inline void add_args(t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6) {
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6);  
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7>
  inline void add_args(t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7) {
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); 
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8>
  inline void add_args(t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8) {
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8);  
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8, class t9>
  inline void add_args(t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8, t9 *a9) {
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8); add_arg(a9); 
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8, class t9, class t10>
  inline void add_args(t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8, t9 *a9, t10 *a10) {
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8); add_arg(a9); add_arg(a10); 
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8, class t9, class t10,
            class t11>
  inline void add_args(t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8, t9 *a9, t10 *a10,
                       t11 *a11) {
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8); add_arg(a9); add_arg(a10); 
    add_arg(a11); 
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8, class t9, class t10,
            class t11, class t12>
  inline void add_args(t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8, t9 *a9, t10 *a10,
                       t11 *a11, t12 *a12) {
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8); add_arg(a9); add_arg(a10); 
    add_arg(a11); add_arg(a12);
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8, class t9, class t10,
            class t11, class t12, class t13>
  inline void add_args(t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8, t9 *a9, t10 *a10,
                       t11 *a11, t12 *a12, t13 *a13) {
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8); add_arg(a9); add_arg(a10); 
    add_arg(a11); add_arg(a12); add_arg(a13);  
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8, class t9, class t10,
            class t11, class t12, class t13, class t14>
  inline void add_args(t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8, t9 *a9, t10 *a10,
                       t11 *a11, t12 *a12, t13 *a13, t14 *a14) {
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8); add_arg(a9); add_arg(a10); 
    add_arg(a11); add_arg(a12); add_arg(a13); add_arg(a14);  
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8, class t9, class t10,
            class t11, class t12, class t13, class t14, class t15>
  inline void add_args(t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8, t9 *a9, t10 *a10,
                       t11 *a11, t12 *a12, t13 *a13, t14 *a14, t15 *a15) {
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8); add_arg(a9); add_arg(a10); 
    add_arg(a11); add_arg(a12); add_arg(a13); add_arg(a14); add_arg(a15); 
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8, class t9, class t10,
            class t11, class t12, class t13, class t14, class t15,
            class t16>
  inline void add_args(t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8, t9 *a9, t10 *a10,
                       t11 *a11, t12 *a12, t13 *a13, t14 *a14, t15 *a15,
                       t16 *a16) {
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8); add_arg(a9); add_arg(a10); 
    add_arg(a11); add_arg(a12); add_arg(a13); add_arg(a14); add_arg(a15); 
    add_arg(a16); 
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8, class t9, class t10,
            class t11, class t12, class t13, class t14, class t15,
            class t16, class t17>
  inline void add_args(t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8, t9 *a9, t10 *a10,
                       t11 *a11, t12 *a12, t13 *a13, t14 *a14, t15 *a15,
                       t16 *a16, t17 *a17) {
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8); add_arg(a9); add_arg(a10); 
    add_arg(a11); add_arg(a12); add_arg(a13); add_arg(a14); add_arg(a15); 
    add_arg(a16); add_arg(a17); 
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8, class t9, class t10,
            class t11, class t12, class t13, class t14, class t15,
            class t16, class t17, class t18>
  inline void add_args(t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8, t9 *a9, t10 *a10,
                       t11 *a11, t12 *a12, t13 *a13, t14 *a14, t15 *a15,
                       t16 *a16, t17 *a17, t18 *a18) {
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8); add_arg(a9); add_arg(a10); 
    add_arg(a11); add_arg(a12); add_arg(a13); add_arg(a14); add_arg(a15); 
    add_arg(a16); add_arg(a17); add_arg(a18); 
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8, class t9, class t10,
            class t11, class t12, class t13, class t14, class t15,
            class t16, class t17, class t18, class t19>
  inline void add_args(t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8, t9 *a9, t10 *a10,
                       t11 *a11, t12 *a12, t13 *a13, t14 *a14, t15 *a15,
                       t16 *a16, t17 *a17, t18 *a18, t19 *a19) {
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8); add_arg(a9); add_arg(a10); 
    add_arg(a11); add_arg(a12); add_arg(a13); add_arg(a14); add_arg(a15); 
    add_arg(a16); add_arg(a17); add_arg(a18); add_arg(a19); 
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8, class t9, class t10,
            class t11, class t12, class t13, class t14, class t15,
            class t16, class t17, class t18, class t19, class t20>
  inline void add_args(t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8, t9 *a9, t10 *a10,
                       t11 *a11, t12 *a12, t13 *a13, t14 *a14, t15 *a15,
                       t16 *a16, t17 *a17, t18 *a18, t19 *a19, t20 *a20) {
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8); add_arg(a9); add_arg(a10); 
    add_arg(a11); add_arg(a12); add_arg(a13); add_arg(a14); add_arg(a15); 
    add_arg(a16); add_arg(a17); add_arg(a18); add_arg(a19); add_arg(a20); 
  }


// ---------------------------------------------------------------------------

  template <class t1>
  inline void run(t1 *a1) {
    clear_args();
    add_arg(a1);
    run();
  }

  template <class t1, class t2>
  inline void run(t1 *a1, t2 *a2) {
    clear_args();
    add_arg(a1); add_arg(a2);
    run();
  }

  template <class t1, class t2, class t3>
  inline void run(t1 *a1, t2 *a2, t3 *a3) {
    clear_args();
    add_arg(a1); add_arg(a2); add_arg(a3);
    run();
  }

  template <class t1, class t2, class t3, class t4>
  inline void run(t1 *a1, t2 *a2, t3 *a3, t4 *a4) {
    clear_args();
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4);
    run();
  }

  template <class t1, class t2, class t3, class t4, class t5>
  inline void run(t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5) {
    clear_args();
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    run();
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6>
  inline void run(t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6) {
    clear_args();
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6);  
    run();
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7>
  inline void run(t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7) {
    clear_args();
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); 
    run();
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8>
  inline void run(t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8) {
    clear_args();
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8);  
    run();
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8, class t9>
  inline void run(t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8, t9 *a9) {
    clear_args();
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8); add_arg(a9); 
    run();
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8, class t9, class t10>
  inline void run(t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8, t9 *a9, t10 *a10) {
    clear_args();
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8); add_arg(a9); add_arg(a10); 
    run();
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8, class t9, class t10,
            class t11>
  inline void run(t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8, t9 *a9, t10 *a10,
                       t11 *a11) {
    clear_args();
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8); add_arg(a9); add_arg(a10); 
    add_arg(a11); 
    run();
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8, class t9, class t10,
            class t11, class t12>
  inline void run(t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8, t9 *a9, t10 *a10,
                       t11 *a11, t12 *a12) {
    clear_args();
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8); add_arg(a9); add_arg(a10); 
    add_arg(a11); add_arg(a12);
    run();
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8, class t9, class t10,
            class t11, class t12, class t13>
  inline void run(t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8, t9 *a9, t10 *a10,
                       t11 *a11, t12 *a12, t13 *a13) {
    clear_args();
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8); add_arg(a9); add_arg(a10); 
    add_arg(a11); add_arg(a12); add_arg(a13);  
    run();
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8, class t9, class t10,
            class t11, class t12, class t13, class t14>
  inline void run(t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8, t9 *a9, t10 *a10,
                       t11 *a11, t12 *a12, t13 *a13, t14 *a14) {
    clear_args();
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8); add_arg(a9); add_arg(a10); 
    add_arg(a11); add_arg(a12); add_arg(a13); add_arg(a14);  
    run();
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8, class t9, class t10,
            class t11, class t12, class t13, class t14, class t15>
  inline void run(t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8, t9 *a9, t10 *a10,
                       t11 *a11, t12 *a12, t13 *a13, t14 *a14, t15 *a15) {
    clear_args();
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8); add_arg(a9); add_arg(a10); 
    add_arg(a11); add_arg(a12); add_arg(a13); add_arg(a14); add_arg(a15); 
    run();
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8, class t9, class t10,
            class t11, class t12, class t13, class t14, class t15,
            class t16>
  inline void run(t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8, t9 *a9, t10 *a10,
                       t11 *a11, t12 *a12, t13 *a13, t14 *a14, t15 *a15,
                       t16 *a16) {
    clear_args();
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8); add_arg(a9); add_arg(a10); 
    add_arg(a11); add_arg(a12); add_arg(a13); add_arg(a14); add_arg(a15); 
    add_arg(a16); 
    run();
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8, class t9, class t10,
            class t11, class t12, class t13, class t14, class t15,
            class t16, class t17>
  inline void run(t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8, t9 *a9, t10 *a10,
                       t11 *a11, t12 *a12, t13 *a13, t14 *a14, t15 *a15,
                       t16 *a16, t17 *a17) {
    clear_args();
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8); add_arg(a9); add_arg(a10); 
    add_arg(a11); add_arg(a12); add_arg(a13); add_arg(a14); add_arg(a15); 
    add_arg(a16); add_arg(a17); 
    run();
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8, class t9, class t10,
            class t11, class t12, class t13, class t14, class t15,
            class t16, class t17, class t18>
  inline void run(t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8, t9 *a9, t10 *a10,
                       t11 *a11, t12 *a12, t13 *a13, t14 *a14, t15 *a15,
                       t16 *a16, t17 *a17, t18 *a18) {
    clear_args();
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8); add_arg(a9); add_arg(a10); 
    add_arg(a11); add_arg(a12); add_arg(a13); add_arg(a14); add_arg(a15); 
    add_arg(a16); add_arg(a17); add_arg(a18); 
    run();
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8, class t9, class t10,
            class t11, class t12, class t13, class t14, class t15,
            class t16, class t17, class t18, class t19>
  inline void run(t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8, t9 *a9, t10 *a10,
                       t11 *a11, t12 *a12, t13 *a13, t14 *a14, t15 *a15,
                       t16 *a16, t17 *a17, t18 *a18, t19 *a19) {
    clear_args();
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8); add_arg(a9); add_arg(a10); 
    add_arg(a11); add_arg(a12); add_arg(a13); add_arg(a14); add_arg(a15); 
    add_arg(a16); add_arg(a17); add_arg(a18); add_arg(a19); 
    run();
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8, class t9, class t10,
            class t11, class t12, class t13, class t14, class t15,
            class t16, class t17, class t18, class t19, class t20>
  inline void run(t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8, t9 *a9, t10 *a10,
                       t11 *a11, t12 *a12, t13 *a13, t14 *a14, t15 *a15,
                       t16 *a16, t17 *a17, t18 *a18, t19 *a19, t20 *a20) {
    clear_args();
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8); add_arg(a9); add_arg(a10); 
    add_arg(a11); add_arg(a12); add_arg(a13); add_arg(a14); add_arg(a15); 
    add_arg(a16); add_arg(a17); add_arg(a18); add_arg(a19); add_arg(a20); 
    run();
  }

// ---------------------------------------------------------------------------

  template <class t1>
  inline void run_cq(command_queue &cq, t1 *a1) {
    clear_args();
    add_arg(a1);
    run(cq);
  }

  template <class t1, class t2>
  inline void run_cq(command_queue &cq, t1 *a1, t2 *a2) {
    clear_args();
    add_arg(a1); add_arg(a2);
    run(cq);
  }

  template <class t1, class t2, class t3>
  inline void run_cq(command_queue &cq, t1 *a1, t2 *a2, t3 *a3) {
    clear_args();
    add_arg(a1); add_arg(a2); add_arg(a3);
    run(cq);
  }

  template <class t1, class t2, class t3, class t4>
  inline void run_cq(command_queue &cq, t1 *a1, t2 *a2, t3 *a3, t4 *a4) {
    clear_args();
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4);
    run(cq);
  }

  template <class t1, class t2, class t3, class t4, class t5>
  inline void run_cq(command_queue &cq, t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5) {
    clear_args();
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    run(cq);
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6>
  inline void run_cq(command_queue &cq, t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6) {
    clear_args();
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6);  
    run(cq);
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7>
  inline void run_cq(command_queue &cq, t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7) {
    clear_args();
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); 
    run(cq);
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8>
  inline void run_cq(command_queue &cq, t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8) {
    clear_args();
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8);  
    run(cq);
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8, class t9>
  inline void run_cq(command_queue &cq, t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8, t9 *a9) {
    clear_args();
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8); add_arg(a9); 
    run(cq);
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8, class t9, class t10>
  inline void run_cq(command_queue &cq, t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8, t9 *a9, t10 *a10) {
    clear_args();
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8); add_arg(a9); add_arg(a10); 
    run(cq);
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8, class t9, class t10,
            class t11>
  inline void run_cq(command_queue &cq, t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8, t9 *a9, t10 *a10,
                       t11 *a11) {
    clear_args();
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8); add_arg(a9); add_arg(a10); 
    add_arg(a11); 
    run(cq);
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8, class t9, class t10,
            class t11, class t12>
  inline void run_cq(command_queue &cq, t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8, t9 *a9, t10 *a10,
                       t11 *a11, t12 *a12) {
    clear_args();
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8); add_arg(a9); add_arg(a10); 
    add_arg(a11); add_arg(a12);
    run(cq);
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8, class t9, class t10,
            class t11, class t12, class t13>
  inline void run_cq(command_queue &cq, t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8, t9 *a9, t10 *a10,
                       t11 *a11, t12 *a12, t13 *a13) {
    clear_args();
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8); add_arg(a9); add_arg(a10); 
    add_arg(a11); add_arg(a12); add_arg(a13);  
    run(cq);
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8, class t9, class t10,
            class t11, class t12, class t13, class t14>
  inline void run_cq(command_queue &cq, t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8, t9 *a9, t10 *a10,
                       t11 *a11, t12 *a12, t13 *a13, t14 *a14) {
    clear_args();
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8); add_arg(a9); add_arg(a10); 
    add_arg(a11); add_arg(a12); add_arg(a13); add_arg(a14);  
    run(cq);
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8, class t9, class t10,
            class t11, class t12, class t13, class t14, class t15>
  inline void run_cq(command_queue &cq, t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8, t9 *a9, t10 *a10,
                       t11 *a11, t12 *a12, t13 *a13, t14 *a14, t15 *a15) {
    clear_args();
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8); add_arg(a9); add_arg(a10); 
    add_arg(a11); add_arg(a12); add_arg(a13); add_arg(a14); add_arg(a15); 
    run(cq);
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8, class t9, class t10,
            class t11, class t12, class t13, class t14, class t15,
            class t16>
  inline void run_cq(command_queue &cq, t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8, t9 *a9, t10 *a10,
                       t11 *a11, t12 *a12, t13 *a13, t14 *a14, t15 *a15,
                       t16 *a16) {
    clear_args();
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8); add_arg(a9); add_arg(a10); 
    add_arg(a11); add_arg(a12); add_arg(a13); add_arg(a14); add_arg(a15); 
    add_arg(a16); 
    run(cq);
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8, class t9, class t10,
            class t11, class t12, class t13, class t14, class t15,
            class t16, class t17>
  inline void run_cq(command_queue &cq, t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8, t9 *a9, t10 *a10,
                       t11 *a11, t12 *a12, t13 *a13, t14 *a14, t15 *a15,
                       t16 *a16, t17 *a17) {
    clear_args();
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8); add_arg(a9); add_arg(a10); 
    add_arg(a11); add_arg(a12); add_arg(a13); add_arg(a14); add_arg(a15); 
    add_arg(a16); add_arg(a17); 
    run(cq);
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8, class t9, class t10,
            class t11, class t12, class t13, class t14, class t15,
            class t16, class t17, class t18>
  inline void run_cq(command_queue &cq, t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8, t9 *a9, t10 *a10,
                       t11 *a11, t12 *a12, t13 *a13, t14 *a14, t15 *a15,
                       t16 *a16, t17 *a17, t18 *a18) {
    clear_args();
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8); add_arg(a9); add_arg(a10); 
    add_arg(a11); add_arg(a12); add_arg(a13); add_arg(a14); add_arg(a15); 
    add_arg(a16); add_arg(a17); add_arg(a18); 
    run(cq);
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8, class t9, class t10,
            class t11, class t12, class t13, class t14, class t15,
            class t16, class t17, class t18, class t19>
  inline void run_cq(command_queue &cq, t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8, t9 *a9, t10 *a10,
                       t11 *a11, t12 *a12, t13 *a13, t14 *a14, t15 *a15,
                       t16 *a16, t17 *a17, t18 *a18, t19 *a19) {
    clear_args();
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8); add_arg(a9); add_arg(a10); 
    add_arg(a11); add_arg(a12); add_arg(a13); add_arg(a14); add_arg(a15); 
    add_arg(a16); add_arg(a17); add_arg(a18); add_arg(a19); 
    run(cq);
  }

  template <class t1, class t2, class t3, class t4, class t5,
            class t6, class t7, class t8, class t9, class t10,
            class t11, class t12, class t13, class t14, class t15,
            class t16, class t17, class t18, class t19, class t20>
  inline void run_cq(command_queue &cq, t1 *a1, t2 *a2, t3 *a3, t4 *a4, t5 *a5,
                       t6 *a6, t7 *a7, t8 *a8, t9 *a9, t10 *a10,
                       t11 *a11, t12 *a12, t13 *a13, t14 *a14, t15 *a15,
                       t16 *a16, t17 *a17, t18 *a18, t19 *a19, t20 *a20) {
    clear_args();
    add_arg(a1); add_arg(a2); add_arg(a3); add_arg(a4); add_arg(a5); 
    add_arg(a6); add_arg(a7); add_arg(a8); add_arg(a9); add_arg(a10); 
    add_arg(a11); add_arg(a12); add_arg(a13); add_arg(a14); add_arg(a15); 
    add_arg(a16); add_arg(a17); add_arg(a18); add_arg(a19); add_arg(a20); 
    run(cq);
  }

