/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef TYPES_H_
#define TYPES_H_

/* Determine default device type and necessary includes */

#include <Kokkos_Core.hpp>

typedef Kokkos::DefaultExecutionSpace execution_space ;

#if ! defined( KOKKOS_ENABLE_CUDA )
  struct double2 {
    double x, y;
    KOKKOS_INLINE_FUNCTION
    double2(double xinit, double yinit) {
      x = xinit;
      y = yinit;
    }
    KOKKOS_INLINE_FUNCTION
    double2() {
      x = 0.0;
      y = 0.0;
    }
    KOKKOS_INLINE_FUNCTION
    double2& operator += (const double2& src) {
      x+=src.x;
      y+=src.y;
      return *this;
    }

    KOKKOS_INLINE_FUNCTION
    volatile double2& operator += (const volatile double2& src) volatile {
      x+=src.x;
      y+=src.y;
      return *this;
    }

  };
#endif

#include <impl/Kokkos_Timer.hpp>

/* Define types used throughout the code */

//Position arrays
typedef Kokkos::View<double*[3], Kokkos::LayoutRight, execution_space>                                   t_x_array ;
typedef t_x_array::HostMirror                                                                        t_x_array_host ;
typedef Kokkos::View<const double*[3], Kokkos::LayoutRight, execution_space>                             t_x_array_const ;
typedef Kokkos::View<const double*[3], Kokkos::LayoutRight, execution_space, Kokkos::MemoryRandomAccess >  t_x_array_randomread ;

//Force array
typedef Kokkos::View<double*[3],  execution_space>                                                       t_f_array ;


//Neighborlist
typedef Kokkos::View<int**, execution_space >                                                            t_neighbors ;
typedef Kokkos::View<const int**, execution_space >                                                      t_neighbors_const ;
typedef Kokkos::View<int*, execution_space, Kokkos::MemoryUnmanaged >                                    t_neighbors_sub ;
typedef Kokkos::View<const int*, execution_space, Kokkos::MemoryUnmanaged >                              t_neighbors_const_sub ;

//1d int array
typedef Kokkos::View<int*, execution_space >                                                             t_int_1d ;
typedef t_int_1d::HostMirror                                                                         t_int_1d_host ;
typedef Kokkos::View<const int*, execution_space >                                                       t_int_1d_const ;
typedef Kokkos::View<int*, execution_space , Kokkos::MemoryUnmanaged>                                    t_int_1d_um ;
typedef Kokkos::View<const int* , execution_space , Kokkos::MemoryUnmanaged>                             t_int_1d_const_um ;

//2d int array
typedef Kokkos::View<int**, Kokkos::LayoutRight, execution_space >                                       t_int_2d ;
typedef t_int_2d::HostMirror                                                                         t_int_2d_host ;

//Scalar ints
typedef Kokkos::View<int[1], Kokkos::LayoutLeft, execution_space>                                        t_int_scalar ;
typedef t_int_scalar::HostMirror                                                                     t_int_scalar_host ;

#endif /* TYPES_H_ */
