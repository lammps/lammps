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

#include <gtest/gtest.h>

#include <iostream>

#include <Kokkos_Core.hpp>

//----------------------------------------------------------------------------

#include <Cuda/Kokkos_Cuda_TaskPolicy.hpp>
#include <impl/Kokkos_ViewTileLeft.hpp>
#include <TestTile.hpp>

//----------------------------------------------------------------------------

#include <TestSharedAlloc.hpp>
#include <TestViewMapping.hpp>

#include <TestViewImpl.hpp>
#include <TestAtomic.hpp>

#include <TestViewAPI.hpp>
#include <TestViewSubview.hpp>
#include <TestViewOfClass.hpp>

#include <TestReduce.hpp>
#include <TestScan.hpp>
#include <TestRange.hpp>
#include <TestTeam.hpp>
#include <TestAggregate.hpp>
#include <TestAggregateReduction.hpp>
#include <TestCompilerMacros.hpp>
#include <TestMemorySpaceTracking.hpp>
#include <TestMemoryPool.hpp>
#include <TestTeamVector.hpp>
#include <TestTemplateMetaFunctions.hpp>
#include <TestCXX11Deduction.hpp>

#include <TestTaskPolicy.hpp>
#include <TestPolicyConstruction.hpp>

//----------------------------------------------------------------------------

class cuda : public ::testing::Test {
protected:
  static void SetUpTestCase();
  static void TearDownTestCase();
};

//----------------------------------------------------------------------------

namespace Test {

TEST_F( cuda, view_impl )
{
  // test_abort<<<32,32>>>(); // Aborts the kernel with CUDA version 4.1 or greater

  test_view_impl< Kokkos::Cuda >();
}

TEST_F( cuda, view_api )
{
  typedef Kokkos::View< const int * , Kokkos::Cuda , Kokkos::MemoryTraits< Kokkos::RandomAccess > > view_texture_managed ;
  typedef Kokkos::View< const int * , Kokkos::Cuda , Kokkos::MemoryTraits< Kokkos::RandomAccess | Kokkos::Unmanaged > > view_texture_unmanaged ;

  TestViewAPI< double , Kokkos::Cuda >();
  TestViewAPI< double , Kokkos::CudaUVMSpace >();

#if 0
  Kokkos::View<double, Kokkos::Cuda > x("x");
  Kokkos::View<double[1], Kokkos::Cuda > y("y");
  // *x = 10 ;
  // x() = 10 ;
  // y[0] = 10 ;
  // y(0) = 10 ;
#endif
}

TEST_F( cuda , view_nested_view )
{
  ::Test::view_nested_view< Kokkos::Cuda >();
}

TEST_F( cuda, view_subview_auto_1d_left ) {
  TestViewSubview::test_auto_1d< Kokkos::LayoutLeft,Kokkos::Cuda >();
}

TEST_F( cuda, view_subview_auto_1d_right ) {
  TestViewSubview::test_auto_1d< Kokkos::LayoutRight,Kokkos::Cuda >();
}

TEST_F( cuda, view_subview_auto_1d_stride ) {
  TestViewSubview::test_auto_1d< Kokkos::LayoutStride,Kokkos::Cuda >();
}

TEST_F( cuda, view_subview_assign_strided ) {
  TestViewSubview::test_1d_strided_assignment< Kokkos::Cuda >();
}

TEST_F( cuda, view_subview_left_0 ) {
  TestViewSubview::test_left_0< Kokkos::CudaUVMSpace >();
}

TEST_F( cuda, view_subview_left_1 ) {
  TestViewSubview::test_left_1< Kokkos::CudaUVMSpace >();
}

TEST_F( cuda, view_subview_left_2 ) {
  TestViewSubview::test_left_2< Kokkos::CudaUVMSpace >();
}

TEST_F( cuda, view_subview_left_3 ) {
  TestViewSubview::test_left_3< Kokkos::CudaUVMSpace >();
}

TEST_F( cuda, view_subview_right_0 ) {
  TestViewSubview::test_right_0< Kokkos::CudaUVMSpace >();
}

TEST_F( cuda, view_subview_right_1 ) {
  TestViewSubview::test_right_1< Kokkos::CudaUVMSpace >();
}

TEST_F( cuda, view_subview_right_3 ) {
  TestViewSubview::test_right_3< Kokkos::CudaUVMSpace >();
}

TEST_F( cuda, view_subview_1d_assign ) {
  TestViewSubview::test_1d_assign< Kokkos::CudaUVMSpace >();
}

TEST_F( cuda, view_subview_2d_from_3d ) {
  TestViewSubview::test_2d_subview_3d< Kokkos::CudaUVMSpace >();
}

TEST_F( cuda, view_subview_2d_from_5d ) {
  TestViewSubview::test_2d_subview_5d< Kokkos::CudaUVMSpace >();
}

}
