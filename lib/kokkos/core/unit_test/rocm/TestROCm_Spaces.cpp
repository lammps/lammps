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
#include <rocm/TestROCm_Category.hpp>

namespace Test {

KOKKOS_INLINE_FUNCTION
void test_abort()
{
  Kokkos::abort( "test_abort" );
}

KOKKOS_INLINE_FUNCTION
void test_rocm_spaces_int_value( int * ptr )
{
  if ( *ptr == 42 ) { *ptr = 2 * 42; }
}

TEST_F( rocm, space_access )
{
  static_assert(
    Kokkos::Impl::MemorySpaceAccess< Kokkos::HostSpace, Kokkos::HostSpace >::assignable, "" );

  static_assert(
    Kokkos::Impl::MemorySpaceAccess< Kokkos::HostSpace, Kokkos::Experimental::ROCmHostPinnedSpace >::assignable, "" );

  static_assert(
    ! Kokkos::Impl::MemorySpaceAccess< Kokkos::HostSpace, Kokkos::Experimental::ROCmSpace >::assignable, "" );

  static_assert(
    ! Kokkos::Impl::MemorySpaceAccess< Kokkos::HostSpace, Kokkos::Experimental::ROCmSpace >::accessible, "" );

  //--------------------------------------

  static_assert(
    Kokkos::Impl::MemorySpaceAccess< Kokkos::Experimental::ROCmSpace, Kokkos::Experimental::ROCmSpace >::assignable, "" );

  static_assert(
    ! Kokkos::Impl::MemorySpaceAccess< Kokkos::Experimental::ROCmSpace, Kokkos::Experimental::ROCmHostPinnedSpace >::assignable, "" );

  static_assert(
    Kokkos::Impl::MemorySpaceAccess< Kokkos::Experimental::ROCmSpace, Kokkos::Experimental::ROCmHostPinnedSpace >::accessible, "" );

  static_assert(
    ! Kokkos::Impl::MemorySpaceAccess< Kokkos::Experimental::ROCmSpace, Kokkos::HostSpace >::assignable, "" );

  static_assert(
    ! Kokkos::Impl::MemorySpaceAccess< Kokkos::Experimental::ROCmSpace, Kokkos::HostSpace >::accessible, "" );

  //--------------------------------------

  static_assert(
    Kokkos::Impl::MemorySpaceAccess< Kokkos::Experimental::ROCmHostPinnedSpace, Kokkos::Experimental::ROCmHostPinnedSpace >::assignable, "" );

  static_assert(
    ! Kokkos::Impl::MemorySpaceAccess< Kokkos::Experimental::ROCmHostPinnedSpace, Kokkos::HostSpace >::assignable, "" );

  static_assert(
    Kokkos::Impl::MemorySpaceAccess< Kokkos::Experimental::ROCmHostPinnedSpace, Kokkos::HostSpace >::accessible, "" );

  static_assert(
    ! Kokkos::Impl::MemorySpaceAccess< Kokkos::Experimental::ROCmHostPinnedSpace, Kokkos::Experimental::ROCmSpace >::assignable, "" );

  static_assert(
    ! Kokkos::Impl::MemorySpaceAccess< Kokkos::Experimental::ROCmHostPinnedSpace, Kokkos::Experimental::ROCmSpace >::accessible, "" );

  //--------------------------------------

  static_assert(
    ! Kokkos::Impl::SpaceAccessibility< Kokkos::Experimental::ROCm, Kokkos::HostSpace >::accessible, "" );

  static_assert(
    Kokkos::Impl::SpaceAccessibility< Kokkos::Experimental::ROCm, Kokkos::Experimental::ROCmSpace >::accessible, "" );

  static_assert(
    Kokkos::Impl::SpaceAccessibility< Kokkos::Experimental::ROCm, Kokkos::Experimental::ROCmHostPinnedSpace >::accessible, "" );

  static_assert(
    ! Kokkos::Impl::SpaceAccessibility< Kokkos::HostSpace, Kokkos::Experimental::ROCmSpace >::accessible, "" );

  static_assert(
    Kokkos::Impl::SpaceAccessibility< Kokkos::HostSpace, Kokkos::Experimental::ROCmHostPinnedSpace >::accessible, "" );

  static_assert(
    std::is_same< Kokkos::Impl::HostMirror< Kokkos::Experimental::ROCmSpace >::Space
                , Kokkos::HostSpace >::value, "" );

  static_assert(
    std::is_same< Kokkos::Impl::HostMirror< Kokkos::Experimental::ROCmHostPinnedSpace >::Space
                , Kokkos::Experimental::ROCmHostPinnedSpace >::value, "" );

  static_assert(
    Kokkos::Impl::SpaceAccessibility
      < Kokkos::Impl::HostMirror< Kokkos::Experimental::ROCm >::Space
      , Kokkos::HostSpace
      >::accessible, "" );

  static_assert(
    Kokkos::Impl::SpaceAccessibility
      < Kokkos::Impl::HostMirror< Kokkos::Experimental::ROCmSpace >::Space
      , Kokkos::HostSpace
      >::accessible, "" );

  static_assert(
    Kokkos::Impl::SpaceAccessibility
      < Kokkos::Impl::HostMirror< Kokkos::Experimental::ROCmHostPinnedSpace >::Space
      , Kokkos::HostSpace
      >::accessible, "" );
}

template< class MemSpace, class ExecSpace >
struct TestViewROCmAccessible {
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

  TestViewROCmAccessible()
    : m_base( "base", N )
    {}

  static void run()
  {
    TestViewROCmAccessible self;
    Kokkos::parallel_for( Kokkos::RangePolicy< typename MemSpace::execution_space, TagInit >( 0, N ), self );
    typename MemSpace::execution_space().fence();

    // Next access is a different execution space, must complete prior kernel.
    long error_count = -1;
    Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace, TagTest >( 0, N ), self, error_count );
    EXPECT_EQ( error_count, 0 );
  }
};

TEST_F( rocm, impl_view_accessible )
{
  TestViewROCmAccessible< Kokkos::Experimental::ROCmSpace, Kokkos::Experimental::ROCm >::run();

  TestViewROCmAccessible< Kokkos::Experimental::ROCmHostPinnedSpace, Kokkos::Experimental::ROCm >::run();
  TestViewROCmAccessible< Kokkos::Experimental::ROCmHostPinnedSpace, Kokkos::HostSpace::execution_space >::run();
}

} // namespace Test
