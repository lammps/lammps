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

#ifndef KOKKOS_TEST_DYNAMICVIEW_HPP
#define KOKKOS_TEST_DYNAMICVIEW_HPP

#include <gtest/gtest.h>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <Kokkos_Core.hpp>

#include <Kokkos_DynamicView.hpp>
#include <impl/Kokkos_Timer.hpp>

namespace Test {

template< typename Scalar , class Space >
struct TestDynamicView
{
  typedef typename Space::execution_space  execution_space ;
  typedef typename Space::memory_space     memory_space ;

  typedef Kokkos::MemoryPool<typename Space::device_type> memory_pool_type;

  typedef Kokkos::Experimental::DynamicView<Scalar*,Space> view_type;
  typedef typename view_type::const_type const_view_type ;

  typedef typename Kokkos::TeamPolicy<execution_space>::member_type member_type ;
  typedef double value_type;

  struct TEST {};
  struct VERIFY {};

  view_type a;
  const unsigned total_size ;

  TestDynamicView( const view_type & arg_a , const unsigned arg_total )
    : a(arg_a), total_size( arg_total ) {}

  KOKKOS_INLINE_FUNCTION
  void operator() ( const TEST , member_type team_member, double& value) const
  {
    const unsigned int team_idx = team_member.league_rank() * team_member.team_size();

    if ( team_member.team_rank() == 0 ) {
      unsigned n = team_idx + team_member.team_size();

      if ( total_size < n ) n = total_size ;

      a.resize_parallel( n );

      if ( a.extent(0) < n ) {
        Kokkos::abort("GrowTest TEST failed resize_parallel");
      }
    }

    // Make sure resize is done for all team members:
    team_member.team_barrier();

    const unsigned int val = team_idx + team_member.team_rank();

    if ( val < total_size ) {
      value += val ;

      a( val ) = val ;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator() ( const VERIFY , member_type team_member, double& value) const
  {
    const unsigned int val =
      team_member.team_rank() + 
      team_member.league_rank() * team_member.team_size();

    if ( val < total_size ) {
    
      if ( val != a(val) ) {
        Kokkos::abort("GrowTest VERIFY failed resize_parallel");
      }

      value += a(val);
    }
  }

  static void run( unsigned arg_total_size )
  {
    typedef Kokkos::TeamPolicy<execution_space,TEST> TestPolicy ;
    typedef Kokkos::TeamPolicy<execution_space,VERIFY> VerifyPolicy ;

// printf("TestDynamicView::run(%d) construct memory pool\n",arg_total_size);

    const size_t total_alloc_size = arg_total_size * sizeof(Scalar) * 1.2 ;
    const size_t superblock = std::min( total_alloc_size , size_t(1000000) );

    memory_pool_type pool( memory_space()
                         , total_alloc_size
                         ,     500 /* min block size in bytes */
                         ,   30000 /* max block size in bytes */
                         , superblock
                         );

// printf("TestDynamicView::run(%d) construct dynamic view\n",arg_total_size);

    view_type da("A",pool,arg_total_size);

    const_view_type ca(da);

// printf("TestDynamicView::run(%d) construct test functor\n",arg_total_size);

    TestDynamicView functor(da,arg_total_size);

    const unsigned team_size = TestPolicy::team_size_recommended(functor);
    const unsigned league_size = ( arg_total_size + team_size - 1 ) / team_size ;

    double reference = 0;
    double result = 0;

// printf("TestDynamicView::run(%d) run functor test\n",arg_total_size);

    Kokkos::parallel_reduce( TestPolicy(league_size,team_size) , functor , reference);
    execution_space::fence();


// printf("TestDynamicView::run(%d) run functor verify\n",arg_total_size);

    Kokkos::parallel_reduce( VerifyPolicy(league_size,team_size) , functor , result );
    execution_space::fence();

// printf("TestDynamicView::run(%d) done\n",arg_total_size);

  }
};

} // namespace Test

#endif /* #ifndef KOKKOS_TEST_DYNAMICVIEW_HPP */

