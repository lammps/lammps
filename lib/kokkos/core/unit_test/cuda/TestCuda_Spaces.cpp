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

#include <Kokkos_Core.hpp>
#include <cuda/TestCuda_Category.hpp>

namespace Test {

__global__
void test_abort()
{
  Kokkos::abort( "test_abort" );
}

__global__
void test_cuda_spaces_int_value( int * ptr )
{
  if ( *ptr == 42 ) { *ptr = 2 * 42; }
}

TEST_F( cuda, space_access )
{
  static_assert(
    Kokkos::Impl::MemorySpaceAccess< Kokkos::HostSpace, Kokkos::HostSpace >::assignable, "" );

  static_assert(
    Kokkos::Impl::MemorySpaceAccess< Kokkos::HostSpace, Kokkos::CudaHostPinnedSpace >::assignable, "" );

  static_assert(
    ! Kokkos::Impl::MemorySpaceAccess< Kokkos::HostSpace, Kokkos::CudaSpace >::assignable, "" );

  static_assert(
    ! Kokkos::Impl::MemorySpaceAccess< Kokkos::HostSpace, Kokkos::CudaSpace >::accessible, "" );

  static_assert(
    ! Kokkos::Impl::MemorySpaceAccess< Kokkos::HostSpace, Kokkos::CudaUVMSpace >::assignable, "" );

  static_assert(
    Kokkos::Impl::MemorySpaceAccess< Kokkos::HostSpace, Kokkos::CudaUVMSpace >::accessible, "" );

  //--------------------------------------

  static_assert(
    Kokkos::Impl::MemorySpaceAccess< Kokkos::CudaSpace, Kokkos::CudaSpace >::assignable, "" );

  static_assert(
    Kokkos::Impl::MemorySpaceAccess< Kokkos::CudaSpace, Kokkos::CudaUVMSpace >::assignable, "" );

  static_assert(
    ! Kokkos::Impl::MemorySpaceAccess< Kokkos::CudaSpace, Kokkos::CudaHostPinnedSpace >::assignable, "" );

  static_assert(
    Kokkos::Impl::MemorySpaceAccess< Kokkos::CudaSpace, Kokkos::CudaHostPinnedSpace >::accessible, "" );

  static_assert(
    ! Kokkos::Impl::MemorySpaceAccess< Kokkos::CudaSpace, Kokkos::HostSpace >::assignable, "" );

  static_assert(
    ! Kokkos::Impl::MemorySpaceAccess< Kokkos::CudaSpace, Kokkos::HostSpace >::accessible, "" );

  //--------------------------------------

  static_assert(
    Kokkos::Impl::MemorySpaceAccess< Kokkos::CudaUVMSpace, Kokkos::CudaUVMSpace >::assignable, "" );

  static_assert(
    ! Kokkos::Impl::MemorySpaceAccess< Kokkos::CudaUVMSpace, Kokkos::CudaSpace >::assignable, "" );

  static_assert(
    Kokkos::Impl::MemorySpaceAccess< Kokkos::CudaUVMSpace, Kokkos::CudaSpace >::accessible, "" );

  static_assert(
    ! Kokkos::Impl::MemorySpaceAccess< Kokkos::CudaUVMSpace, Kokkos::HostSpace >::assignable, "" );

  static_assert(
    ! Kokkos::Impl::MemorySpaceAccess< Kokkos::CudaUVMSpace, Kokkos::HostSpace >::accessible, "" );

  static_assert(
    ! Kokkos::Impl::MemorySpaceAccess< Kokkos::CudaUVMSpace, Kokkos::CudaHostPinnedSpace >::assignable, "" );

  static_assert(
    Kokkos::Impl::MemorySpaceAccess< Kokkos::CudaUVMSpace, Kokkos::CudaHostPinnedSpace >::accessible, "" );

  //--------------------------------------

  static_assert(
    Kokkos::Impl::MemorySpaceAccess< Kokkos::CudaHostPinnedSpace, Kokkos::CudaHostPinnedSpace >::assignable, "" );

  static_assert(
    ! Kokkos::Impl::MemorySpaceAccess< Kokkos::CudaHostPinnedSpace, Kokkos::HostSpace >::assignable, "" );

  static_assert(
    Kokkos::Impl::MemorySpaceAccess< Kokkos::CudaHostPinnedSpace, Kokkos::HostSpace >::accessible, "" );

  static_assert(
    ! Kokkos::Impl::MemorySpaceAccess< Kokkos::CudaHostPinnedSpace, Kokkos::CudaSpace >::assignable, "" );

  static_assert(
    ! Kokkos::Impl::MemorySpaceAccess< Kokkos::CudaHostPinnedSpace, Kokkos::CudaSpace >::accessible, "" );

  static_assert(
    ! Kokkos::Impl::MemorySpaceAccess< Kokkos::CudaHostPinnedSpace, Kokkos::CudaUVMSpace >::assignable, "" );

  static_assert(
    Kokkos::Impl::MemorySpaceAccess< Kokkos::CudaHostPinnedSpace, Kokkos::CudaUVMSpace >::accessible, "" );

  //--------------------------------------

  static_assert(
    ! Kokkos::Impl::SpaceAccessibility< Kokkos::Cuda, Kokkos::HostSpace >::accessible, "" );

  static_assert(
    Kokkos::Impl::SpaceAccessibility< Kokkos::Cuda, Kokkos::CudaSpace >::accessible, "" );

  static_assert(
    Kokkos::Impl::SpaceAccessibility< Kokkos::Cuda, Kokkos::CudaUVMSpace >::accessible, "" );

  static_assert(
    Kokkos::Impl::SpaceAccessibility< Kokkos::Cuda, Kokkos::CudaHostPinnedSpace >::accessible, "" );

  static_assert(
    ! Kokkos::Impl::SpaceAccessibility< Kokkos::HostSpace, Kokkos::CudaSpace >::accessible, "" );

  static_assert(
    Kokkos::Impl::SpaceAccessibility< Kokkos::HostSpace, Kokkos::CudaUVMSpace >::accessible, "" );

  static_assert(
    Kokkos::Impl::SpaceAccessibility< Kokkos::HostSpace, Kokkos::CudaHostPinnedSpace >::accessible, "" );

  static_assert(
    std::is_same< Kokkos::Impl::HostMirror< Kokkos::CudaSpace >::Space
                , Kokkos::HostSpace >::value, "" );

  static_assert(
    std::is_same< Kokkos::Impl::HostMirror< Kokkos::CudaUVMSpace >::Space
                , Kokkos::Device< Kokkos::HostSpace::execution_space
                                , Kokkos::CudaUVMSpace > >::value, "" );

  static_assert(
    std::is_same< Kokkos::Impl::HostMirror< Kokkos::CudaHostPinnedSpace >::Space
                , Kokkos::CudaHostPinnedSpace >::value, "" );

  static_assert(
    std::is_same< Kokkos::Device< Kokkos::HostSpace::execution_space
                                , Kokkos::CudaUVMSpace >
                , Kokkos::Device< Kokkos::HostSpace::execution_space
                                , Kokkos::CudaUVMSpace > >::value, "" );

  static_assert(
    Kokkos::Impl::SpaceAccessibility
      < Kokkos::Impl::HostMirror< Kokkos::Cuda >::Space
      , Kokkos::HostSpace
      >::accessible, "" );

  static_assert(
    Kokkos::Impl::SpaceAccessibility
      < Kokkos::Impl::HostMirror< Kokkos::CudaSpace >::Space
      , Kokkos::HostSpace
      >::accessible, "" );

  static_assert(
    Kokkos::Impl::SpaceAccessibility
      < Kokkos::Impl::HostMirror< Kokkos::CudaUVMSpace >::Space
      , Kokkos::HostSpace
      >::accessible, "" );

  static_assert(
    Kokkos::Impl::SpaceAccessibility
      < Kokkos::Impl::HostMirror< Kokkos::CudaHostPinnedSpace >::Space
      , Kokkos::HostSpace
      >::accessible, "" );
}

TEST_F( cuda, uvm )
{
  if ( Kokkos::CudaUVMSpace::available() ) {
    int * uvm_ptr = (int*) Kokkos::kokkos_malloc< Kokkos::CudaUVMSpace >( "uvm_ptr", sizeof( int ) );

    *uvm_ptr = 42;

    Kokkos::Cuda().fence();
    test_cuda_spaces_int_value<<< 1, 1 >>>( uvm_ptr );
    Kokkos::Cuda().fence();

    EXPECT_EQ( *uvm_ptr, int( 2 * 42 ) );

    Kokkos::kokkos_free< Kokkos::CudaUVMSpace >( uvm_ptr );
  }
}

/* Removing UVM Allocs Test due to added time to complete overall unit test
 * The issue verified with this unit test appears to no longer be an 
 * problem.  Refer to github issue 1880 for more details
 *
TEST_F( cuda, uvm_num_allocs )
{
  // The max number of UVM allocations allowed is 65536.
  #define MAX_NUM_ALLOCS 65536

  if ( Kokkos::CudaUVMSpace::available() ) {
    struct TestMaxUVMAllocs {

      using view_type         = Kokkos::View< double*, Kokkos::CudaUVMSpace >;
      using view_of_view_type = Kokkos::View< view_type[ MAX_NUM_ALLOCS ]
                                            , Kokkos::CudaUVMSpace >;

      TestMaxUVMAllocs() : view_allocs_test( "view_allocs_test" )
      {
        for ( auto i = 0; i < MAX_NUM_ALLOCS; ++i ) {

          // Kokkos will throw a runtime exception if an attempt is made to
          // allocate more than the maximum number of uvm allocations.

          // In this test, the max num of allocs occurs when i = MAX_NUM_ALLOCS - 1
          // since the 'outer' view counts as one UVM allocation, leaving
          // 65535 possible UVM allocations, that is 'i in [0, 65535)'.

          // The test will catch the exception thrown in this case and continue.

          if ( i == ( MAX_NUM_ALLOCS - 1 ) ) {
            EXPECT_ANY_THROW( { view_allocs_test( i ) = view_type( "inner_view", 1 ); } );
          }
          else {
            if ( i < MAX_NUM_ALLOCS - 1000 ) {
              EXPECT_NO_THROW( { view_allocs_test( i ) = view_type( "inner_view", 1 ); } );
            } else { // This might or might not throw depending on compilation options.
              try {
                view_allocs_test( i ) = view_type( "inner_view", 1 );
              }
              catch ( ... ) {}
            }
          }

        } // End allocation for loop.

        for ( auto i = 0; i < MAX_NUM_ALLOCS - 1; ++i ) {

          view_allocs_test( i ) = view_type();

        } // End deallocation for loop.

        view_allocs_test = view_of_view_type(); // Deallocate the view of views.
      }

      // Member.
      view_of_view_type view_allocs_test;
    };

    // Trigger the test via the TestMaxUVMAllocs constructor.
    TestMaxUVMAllocs();
  }

  #undef MAX_NUM_ALLOCS
}
*/

template< class MemSpace, class ExecSpace >
struct TestViewCudaAccessible {
  enum { N = 1000 };

  using V = Kokkos::View< double*, MemSpace >;

  V m_base;

  struct TagInit {};
  struct TagTest {};

  KOKKOS_INLINE_FUNCTION
  void operator()( const TagInit &, const int i ) const { m_base[i] = i + 1; }

  KOKKOS_INLINE_FUNCTION
  void operator()( const TagTest &, const int i, long & error_count ) const
  { if ( m_base[i] != i + 1 ) ++error_count; }

  TestViewCudaAccessible()
    : m_base( "base", N )
    {}

  static void run()
  {
    TestViewCudaAccessible self;
    Kokkos::parallel_for( Kokkos::RangePolicy< typename MemSpace::execution_space, TagInit >( 0, N ), self );
    typename MemSpace::execution_space().fence();

    // Next access is a different execution space, must complete prior kernel.
    long error_count = -1;
    Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace, TagTest >( 0, N ), self, error_count );
    EXPECT_EQ( error_count, 0 );
  }
};

TEST_F( cuda, impl_view_accessible )
{
  TestViewCudaAccessible< Kokkos::CudaSpace, Kokkos::Cuda >::run();

  TestViewCudaAccessible< Kokkos::CudaUVMSpace, Kokkos::Cuda >::run();
  TestViewCudaAccessible< Kokkos::CudaUVMSpace, Kokkos::HostSpace::execution_space >::run();

  TestViewCudaAccessible< Kokkos::CudaHostPinnedSpace, Kokkos::Cuda >::run();
  TestViewCudaAccessible< Kokkos::CudaHostPinnedSpace, Kokkos::HostSpace::execution_space >::run();
}

template< class MemSpace >
struct TestViewCudaTexture {
  enum { N = 1000 };

  using V = Kokkos::View< double*, MemSpace >;
  using T = Kokkos::View< const double*, MemSpace, Kokkos::MemoryRandomAccess >;

  V m_base;
  T m_tex;

  struct TagInit {};
  struct TagTest {};

  KOKKOS_INLINE_FUNCTION
  void operator()( const TagInit &, const int i ) const { m_base[i] = i + 1; }

  KOKKOS_INLINE_FUNCTION
  void operator()( const TagTest &, const int i, long & error_count ) const
  { if ( m_tex[i] != i + 1 ) ++error_count; }

  TestViewCudaTexture()
    : m_base( "base", N )
    , m_tex( m_base )
    {}

  static void run()
  {
    EXPECT_TRUE( ( std::is_same< typename V::reference_type, double & >::value ) );
    EXPECT_TRUE( ( std::is_same< typename T::reference_type, const double >::value ) );

    EXPECT_TRUE(  V::reference_type_is_lvalue_reference ); // An ordinary view.
    EXPECT_FALSE( T::reference_type_is_lvalue_reference ); // Texture fetch returns by value.

    TestViewCudaTexture self;
    Kokkos::parallel_for( Kokkos::RangePolicy< Kokkos::Cuda, TagInit >( 0, N ), self );

    long error_count = -1;
    Kokkos::parallel_reduce( Kokkos::RangePolicy< Kokkos::Cuda, TagTest >( 0, N ), self, error_count );
    EXPECT_EQ( error_count, 0 );
  }
};

TEST_F( cuda, impl_view_texture )
{
  TestViewCudaTexture< Kokkos::CudaSpace >::run();
  TestViewCudaTexture< Kokkos::CudaUVMSpace >::run();
}

} // namespace Test
