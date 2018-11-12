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
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <cstdio>

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

namespace Test {

namespace {

template < typename ExecSpace >
struct TestViewCopy {

  using InExecSpace = ExecSpace;

  static void test_view_copy()
  {
#if defined( KOKKOS_ENABLE_CUDA ) || defined( KOKKOS_ENABLE_ROCM )
   // ExecSpace = CudaUVM, CudaHostPinned
   // This test will fail at runtime with an illegal memory access if something goes wrong
   // Test 1: deep_copy from host_mirror_space to ExecSpace and ExecSpace back to host_mirror_space
   {
    const int dim0 = 4;
    const int dim1 = 2;
    const int dim2 = 3;

    typedef Kokkos::View<double****,InExecSpace> Rank4ViewType;
    Rank4ViewType view_4;
    view_4 = Rank4ViewType("view_4", dim0, dim1, dim2, dim2);

    typedef typename Kokkos::Impl::is_space<InExecSpace>::host_mirror_space::execution_space host_space_type;
    Kokkos::View<double**,Kokkos::LayoutLeft,host_space_type> srcView("srcView", dim2, dim2);

    // Strided dst view
    auto dstView = Kokkos::subview(view_4, 0, 0, Kokkos::ALL(), Kokkos::ALL());

    // host_mirror_space to ExecSpace
    Kokkos::deep_copy( dstView, srcView );
    Kokkos::fence();

    // ExecSpace to host_mirror_space 
    Kokkos::deep_copy( srcView, dstView );
    Kokkos::fence();
   }

   // Test 2: deep_copy from Cuda to ExecSpace and ExecSpace back to Cuda
   {
    const int dim0 = 4;
    const int dim1 = 2;
    const int dim2 = 3;

    typedef Kokkos::View<double****,InExecSpace> Rank4ViewType;
    Rank4ViewType view_4;
    view_4 = Rank4ViewType("view_4", dim0, dim1, dim2, dim2);

#if defined( KOKKOS_ENABLE_CUDA )
    typedef Kokkos::Cuda space_type;
#endif
#if defined( KOKKOS_ENABLE_ROCM )
    typedef Kokkos::Experimental::ROCm space_type;
#endif
    Kokkos::View<double**,Kokkos::LayoutLeft,space_type> srcView("srcView", dim2, dim2);

    // Strided dst view
    auto dstView = Kokkos::subview(view_4, 0, 0, Kokkos::ALL(), Kokkos::ALL());

    // Cuda to ExecSpace
    Kokkos::deep_copy( dstView, srcView );
    Kokkos::fence();

    // ExecSpace to Cuda
    Kokkos::deep_copy( srcView, dstView );
    Kokkos::fence();
   }

   // Test 3: deep_copy from host_space to ExecSpace and ExecSpace back to host_space
   {
    const int dim0 = 4;
    const int dim1 = 2;
    const int dim2 = 3;

    typedef Kokkos::View<double****,InExecSpace> Rank4ViewType;
    Rank4ViewType view_4;
    view_4 = Rank4ViewType("view_4", dim0, dim1, dim2, dim2);

    typedef Kokkos::HostSpace host_space_type;
    Kokkos::View<double**,Kokkos::LayoutLeft,host_space_type> srcView("srcView", dim2, dim2);

    // Strided dst view
    auto dstView = Kokkos::subview(view_4, 0, 0, Kokkos::ALL(), Kokkos::ALL());

    // host_space to ExecSpace
    Kokkos::deep_copy( dstView, srcView );
    Kokkos::fence();

    // ExecSpace to host_space 
    Kokkos::deep_copy( srcView, dstView );
    Kokkos::fence();
   }
#endif
  } // end test_view_copy

}; // end struct

} // namespace

TEST_F( TEST_CATEGORY , view_copy_tests ) {
  //Only include this file to be compiled with CudaUVM and CudaHostPinned
  TestViewCopy< TEST_EXECSPACE >::test_view_copy();
}

} // namespace Test
