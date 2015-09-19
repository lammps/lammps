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

#include <Kokkos_Core.hpp>
#include <stdexcept>
#include <sstream>
#include <iostream>

/*--------------------------------------------------------------------------*/

namespace Test {

namespace {
volatile int nested_view_count ;
}

template< class Space >
class NestedView {
private:
  Kokkos::View<int*,Space> member ;

public:

  KOKKOS_INLINE_FUNCTION
  NestedView()
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
    : member("member",2)
    { Kokkos::atomic_increment( & nested_view_count ); }
#else
    : member(){}
#endif

  ~NestedView()
#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
    { Kokkos::atomic_decrement( & nested_view_count ); }
#else
    {}
#endif

};


template< class Space >
void view_nested_view()
{
  ASSERT_EQ( 0 , nested_view_count );
  {
    Kokkos::View< NestedView<Space> * , Space > a("a_nested_view",2);
    ASSERT_EQ( 2 , nested_view_count );
    Kokkos::View< NestedView<Space> * , Space > b("b_nested_view",2);
    ASSERT_EQ( 4 , nested_view_count );
  }
  // ASSERT_EQ( 0 , nested_view_count );
}

}

namespace Kokkos {
namespace Impl {

template< class ExecSpace , class S >
struct ViewDefaultConstruct< ExecSpace , Test::NestedView<S> , true >
{
  typedef Test::NestedView<S> type ;
  type * const m_ptr ;

  KOKKOS_FORCEINLINE_FUNCTION
  void operator()( const typename ExecSpace::size_type& i ) const
    { new(m_ptr+i) type(); }

  ViewDefaultConstruct( type * pointer , size_t capacity )
    : m_ptr( pointer )
    {
      Kokkos::RangePolicy< ExecSpace > range( 0 , capacity );
      parallel_for( range , *this );
      ExecSpace::fence();
    }
};

} // namespace Impl
} // namespace Kokkos

/*--------------------------------------------------------------------------*/

