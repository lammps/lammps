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

#include <TestMDRange.hpp>

//----------------------------------------------------------------------------

class cuda : public ::testing::Test {
protected:
  static void SetUpTestCase();
  static void TearDownTestCase();
};

void cuda::SetUpTestCase()
  {
    Kokkos::Cuda::print_configuration( std::cout );
    Kokkos::HostSpace::execution_space::initialize();
    Kokkos::Cuda::initialize( Kokkos::Cuda::SelectDevice(0) );
  }

void cuda::TearDownTestCase()
  {
    Kokkos::Cuda::finalize();
    Kokkos::HostSpace::execution_space::finalize();
  }

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Test {

__global__
void test_abort()
{
  Kokkos::Impl::VerifyExecutionCanAccessMemorySpace<
    Kokkos::CudaSpace ,
    Kokkos::HostSpace >::verify();
}

__global__
void test_cuda_spaces_int_value( int * ptr )
{
  if ( *ptr == 42 ) { *ptr = 2 * 42 ; }
}

TEST_F( cuda , md_range ) {
  TestMDRange_2D< Kokkos::Cuda >::test_for2(100,100);

  TestMDRange_3D< Kokkos::Cuda >::test_for3(100,100,100);
}

TEST_F( cuda , compiler_macros )
{
  ASSERT_TRUE( ( TestCompilerMacros::Test< Kokkos::Cuda >() ) );
}

TEST_F( cuda , memory_space )
{
  TestMemorySpace< Kokkos::Cuda >();
}

TEST_F( cuda, uvm )
{
  if ( Kokkos::CudaUVMSpace::available() ) {

    int * uvm_ptr = (int*) Kokkos::kokkos_malloc< Kokkos::CudaUVMSpace >("uvm_ptr",sizeof(int));

    *uvm_ptr = 42 ;

    Kokkos::Cuda::fence();
    test_cuda_spaces_int_value<<<1,1>>>(uvm_ptr);
    Kokkos::Cuda::fence();

    EXPECT_EQ( *uvm_ptr, int(2*42) );

    Kokkos::kokkos_free< Kokkos::CudaUVMSpace >(uvm_ptr );
  }
}

//----------------------------------------------------------------------------

TEST_F( cuda , impl_shared_alloc )
{
  test_shared_alloc< Kokkos::CudaSpace , Kokkos::HostSpace::execution_space >();
  test_shared_alloc< Kokkos::CudaUVMSpace , Kokkos::HostSpace::execution_space >();
  test_shared_alloc< Kokkos::CudaHostPinnedSpace , Kokkos::HostSpace::execution_space >();
}

TEST_F( cuda, policy_construction) {
  TestRangePolicyConstruction< Kokkos::Cuda >();
  TestTeamPolicyConstruction< Kokkos::Cuda >();
}

TEST_F( cuda , impl_view_mapping )
{
  test_view_mapping< Kokkos::Cuda >();
  test_view_mapping< Kokkos::CudaUVMSpace >();
  test_view_mapping_subview< Kokkos::Cuda >();
  test_view_mapping_subview< Kokkos::CudaUVMSpace >();
  test_view_mapping_operator< Kokkos::Cuda >();
  test_view_mapping_operator< Kokkos::CudaUVMSpace >();
  TestViewMappingAtomic< Kokkos::Cuda >::run();
}

TEST_F( cuda , view_of_class )
{
  TestViewMappingClassValue< Kokkos::CudaSpace >::run();
  TestViewMappingClassValue< Kokkos::CudaUVMSpace >::run();
}

template< class MemSpace >
struct TestViewCudaTexture {

  enum { N = 1000 };

  using V = Kokkos::Experimental::View<double*,MemSpace> ;
  using T = Kokkos::Experimental::View<const double*, MemSpace, Kokkos::MemoryRandomAccess > ;

  V m_base ;
  T m_tex ;

  struct TagInit {};
  struct TagTest {};

  KOKKOS_INLINE_FUNCTION
  void operator()( const TagInit & , const int i ) const { m_base[i] = i + 1 ; }

  KOKKOS_INLINE_FUNCTION
  void operator()( const TagTest & , const int i , long & error_count ) const
    { if ( m_tex[i] != i + 1 ) ++error_count ; }

  TestViewCudaTexture()
    : m_base("base",N)
    , m_tex( m_base )
    {}

  static void run()
    {
      EXPECT_TRUE( ( std::is_same< typename V::reference_type
                                 , double &
                                 >::value ) );

      EXPECT_TRUE( ( std::is_same< typename T::reference_type
                                 , const double
                                 >::value ) );

      EXPECT_TRUE(  V::reference_type_is_lvalue_reference ); // An ordinary view
      EXPECT_FALSE( T::reference_type_is_lvalue_reference ); // Texture fetch returns by value

      TestViewCudaTexture self ;
      Kokkos::parallel_for( Kokkos::RangePolicy< Kokkos::Cuda , TagInit >(0,N) , self );
      long error_count = -1 ;
      Kokkos::parallel_reduce( Kokkos::RangePolicy< Kokkos::Cuda , TagTest >(0,N) , self , error_count );
      EXPECT_EQ( error_count , 0 );
    }
};

TEST_F( cuda , impl_view_texture )
{
  TestViewCudaTexture< Kokkos::CudaSpace >::run();
  TestViewCudaTexture< Kokkos::CudaUVMSpace >::run();
}

template< class MemSpace , class ExecSpace >
struct TestViewCudaAccessible {

  enum { N = 1000 };

  using V = Kokkos::Experimental::View<double*,MemSpace> ;

  V m_base ;

  struct TagInit {};
  struct TagTest {};

  KOKKOS_INLINE_FUNCTION
  void operator()( const TagInit & , const int i ) const { m_base[i] = i + 1 ; }

  KOKKOS_INLINE_FUNCTION
  void operator()( const TagTest & , const int i , long & error_count ) const
    { if ( m_base[i] != i + 1 ) ++error_count ; }

  TestViewCudaAccessible()
    : m_base("base",N)
    {}

  static void run()
    {
      TestViewCudaAccessible self ;
      Kokkos::parallel_for( Kokkos::RangePolicy< typename MemSpace::execution_space , TagInit >(0,N) , self );
      MemSpace::execution_space::fence();
      // Next access is a different execution space, must complete prior kernel.
      long error_count = -1 ;
      Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace , TagTest >(0,N) , self , error_count );
      EXPECT_EQ( error_count , 0 );
    }
};

TEST_F( cuda , impl_view_accessible )
{
  TestViewCudaAccessible< Kokkos::CudaSpace , Kokkos::Cuda >::run();

  TestViewCudaAccessible< Kokkos::CudaUVMSpace , Kokkos::Cuda >::run();
  TestViewCudaAccessible< Kokkos::CudaUVMSpace , Kokkos::HostSpace::execution_space >::run();

  TestViewCudaAccessible< Kokkos::CudaHostPinnedSpace , Kokkos::Cuda >::run();
  TestViewCudaAccessible< Kokkos::CudaHostPinnedSpace , Kokkos::HostSpace::execution_space >::run();
}

}
