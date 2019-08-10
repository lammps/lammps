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

#include <gtest/gtest.h>

#include <stdexcept>
#include <sstream>
#include <iostream>

#include <Kokkos_Core.hpp>

namespace Test {

/*--------------------------------------------------------------------------*/

template< class Space >
struct TestViewMappingAtomic {
  typedef typename Space::execution_space ExecSpace;
  typedef typename Space::memory_space    MemSpace;

  typedef Kokkos::MemoryTraits< Kokkos::Atomic >  mem_trait;

  typedef Kokkos::View< int *, ExecSpace > T;
  typedef Kokkos::View< int *, ExecSpace, mem_trait >  T_atom;

  T      x;
  T_atom x_atom;

  enum { N = 100000};

  struct TagInit {};
  struct TagUpdate {};
  struct TagVerify {};

  KOKKOS_INLINE_FUNCTION
  void operator()( const TagInit &, const int i ) const
  { x( i ) = i; }

  KOKKOS_INLINE_FUNCTION
  void operator()( const TagUpdate &, const int i ) const
  { x_atom( i % 2 ) += 1; }

  KOKKOS_INLINE_FUNCTION
  void operator()( const TagVerify &, const int i, long & error_count ) const
  {
     if ( i < 2 ) { if ( x( i ) != int( i + N / 2 ) ) ++error_count; }
     else         { if ( x( i ) != int( i ) ) ++error_count; }
  }

  TestViewMappingAtomic()
    : x( "x", N )
    , x_atom( x )
    {}

  void run() {

    ASSERT_TRUE( T::reference_type_is_lvalue_reference );
    ASSERT_FALSE( T_atom::reference_type_is_lvalue_reference );

    Kokkos::parallel_for( Kokkos::RangePolicy< ExecSpace, TagInit >  ( 0, N ), *this );
    Kokkos::parallel_for( Kokkos::RangePolicy< ExecSpace, TagUpdate >( 0, N ), *this );

    long error_count = -1;

    Kokkos::parallel_reduce( Kokkos::RangePolicy< ExecSpace, TagVerify >( 0, N ), *this, error_count );

    ASSERT_EQ( 0, error_count );

    typename T_atom::HostMirror x_host = Kokkos::create_mirror_view( x );
    Kokkos::deep_copy( x_host, x );

    error_count = -1;

    Kokkos::parallel_reduce( Kokkos::RangePolicy< Kokkos::DefaultHostExecutionSpace, TagVerify >( 0, N ), 
      [=] ( const TagVerify &, const int i, long & tmp_error_count )
    {
      if ( i < 2 ) {
        if ( x_host( i ) != int( i + N / 2 ) ) ++tmp_error_count ;
      }
      else {
        if ( x_host( i ) != int( i ) ) ++tmp_error_count ;
      }
    }, error_count);

    ASSERT_EQ( 0 , error_count );
    Kokkos::deep_copy( x, x_host );
  }
};

TEST_F( TEST_CATEGORY , view_mapping_atomic )
{
  TestViewMappingAtomic< TEST_EXECSPACE > f;
  f.run();
}

}

/*--------------------------------------------------------------------------*/

namespace Test {

struct MappingClassValueType {
    KOKKOS_INLINE_FUNCTION
    MappingClassValueType() 
    {
#if 0
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_CUDA )
      printf( "TestViewMappingClassValue construct on Cuda\n" );
#elif defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
      printf( "TestViewMappingClassValue construct on Host\n" );
#else
      printf( "TestViewMappingClassValue construct unknown\n" );
#endif
#endif
    }
    KOKKOS_INLINE_FUNCTION
    ~MappingClassValueType()
    {
#if 0
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_CUDA )
      printf( "TestViewMappingClassValue destruct on Cuda\n" );
#elif defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
      printf( "TestViewMappingClassValue destruct on Host\n" );
#else
      printf( "TestViewMappingClassValue destruct unknown\n" );
#endif
#endif
    }
  };

template< class Space >
void test_view_mapping_class_value()
{
  typedef typename Space::execution_space ExecSpace;

  ExecSpace().fence();
  {
    Kokkos::View< MappingClassValueType, ExecSpace > a( "a" );
    ExecSpace().fence();
  }
  ExecSpace().fence();
}

TEST_F( TEST_CATEGORY , view_mapping_class_value )
{
  test_view_mapping_class_value< TEST_EXECSPACE >();
}

}

/*--------------------------------------------------------------------------*/

namespace Test {

TEST_F( TEST_CATEGORY , view_mapping_assignable )
{
  typedef TEST_EXECSPACE exec_space ;

  { // Assignment of rank-0 Left = Right
    typedef Kokkos::ViewTraits<int,Kokkos::LayoutLeft, exec_space> dst_traits ;
    typedef Kokkos::ViewTraits<int,Kokkos::LayoutRight,exec_space> src_traits ;
    typedef Kokkos::Impl::ViewMapping<dst_traits,src_traits,void> mapping ;
    static_assert( mapping::is_assignable , "" );

    Kokkos::View<int,Kokkos::LayoutRight,exec_space> src ;
    Kokkos::View<int,Kokkos::LayoutLeft,exec_space> dst( src );
    dst = src ;
  }

  { // Assignment of rank-0 Right = Left
    typedef Kokkos::ViewTraits<int,Kokkos::LayoutRight,exec_space> dst_traits ;
    typedef Kokkos::ViewTraits<int,Kokkos::LayoutLeft, exec_space> src_traits ;
    typedef Kokkos::Impl::ViewMapping<dst_traits,src_traits,void> mapping ;
    static_assert( mapping::is_assignable , "" );

    Kokkos::View<int,Kokkos::LayoutLeft,exec_space> src ;
    Kokkos::View<int,Kokkos::LayoutRight,exec_space> dst( src );
    dst = src ;
  }

  { // Assignment of rank-1 Left = Right
    typedef Kokkos::ViewTraits<int*,Kokkos::LayoutLeft, exec_space> dst_traits ;
    typedef Kokkos::ViewTraits<int*,Kokkos::LayoutRight,exec_space> src_traits ;
    typedef Kokkos::Impl::ViewMapping<dst_traits,src_traits,void> mapping ;
    static_assert( mapping::is_assignable , "" );

    Kokkos::View<int*,Kokkos::LayoutRight,exec_space> src ;
    Kokkos::View<int*,Kokkos::LayoutLeft,exec_space> dst( src );
    dst = src ;
  }

  { // Assignment of rank-1 Right = Left
    typedef Kokkos::ViewTraits<int*,Kokkos::LayoutRight,exec_space> dst_traits ;
    typedef Kokkos::ViewTraits<int*,Kokkos::LayoutLeft, exec_space> src_traits ;
    typedef Kokkos::Impl::ViewMapping<dst_traits,src_traits,void> mapping ;
    static_assert( mapping::is_assignable , "" );

    Kokkos::View<int*,Kokkos::LayoutLeft,exec_space> src ;
    Kokkos::View<int*,Kokkos::LayoutRight,exec_space> dst( src );
    dst = src ;
  }

  { // Assignment of rank-2 Left = Right
    typedef Kokkos::ViewTraits<int**,Kokkos::LayoutLeft, exec_space> dst_traits ;
    typedef Kokkos::ViewTraits<int**,Kokkos::LayoutRight,exec_space> src_traits ;
    typedef Kokkos::Impl::ViewMapping<dst_traits,src_traits,void> mapping ;
    static_assert( ! mapping::is_assignable , "" );
  }

  { // Assignment of rank-2 Right = Left
    typedef Kokkos::ViewTraits<int**,Kokkos::LayoutRight,exec_space> dst_traits ;
    typedef Kokkos::ViewTraits<int**,Kokkos::LayoutLeft, exec_space> src_traits ;
    typedef Kokkos::Impl::ViewMapping<dst_traits,src_traits,void> mapping ;
    static_assert( ! mapping::is_assignable , "" );
  }

}

}

