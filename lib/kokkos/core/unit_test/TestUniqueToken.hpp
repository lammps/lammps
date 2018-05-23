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

#include <iostream>

#include <Kokkos_Core.hpp>

namespace Test {

template< class Space >
class TestUniqueToken
{
public:
  typedef typename Space::execution_space  execution_space;
  typedef Kokkos::View< int * , execution_space > view_type ;

  Kokkos::Experimental::UniqueToken< execution_space , Kokkos::Experimental::UniqueTokenScope::Global > tokens ;

  view_type verify ;
  view_type counts ;
  view_type errors ;

  KOKKOS_INLINE_FUNCTION
  void operator()( long ) const
  {
    const int32_t t = tokens.acquire();

    bool ok = true ;

    ok = ok && 0 <= t ;
    ok = ok && t < tokens.size();
    ok = ok && 0 == Kokkos::atomic_fetch_add( & verify(t) , 1 );

    Kokkos::atomic_fetch_add( & counts(t) , 1 );

    ok = ok && 1 == Kokkos::atomic_fetch_add( & verify(t) , -1 );

    if ( ! ok ) { Kokkos::atomic_fetch_add( & errors(0) , 1 ) ; }

    tokens.release(t);
  }

  TestUniqueToken()
    : tokens( execution_space() )
    , verify( "TestUniqueTokenVerify" , tokens.size() )
    , counts( "TestUniqueTokenCounts" , tokens.size() )
    , errors( "TestUniqueTokenErrors" , 1 )
    {}

  static void run()
    {
      using policy = Kokkos::RangePolicy<execution_space> ;

      TestUniqueToken self ;

      {
        const int duplicate = 100 ;
        const long n = duplicate * self.tokens.size();

        Kokkos::parallel_for( policy(0,n) , self );
        Kokkos::parallel_for( policy(0,n) , self );
        Kokkos::parallel_for( policy(0,n) , self );
        Kokkos::fence();
      }

      typename view_type::HostMirror host_counts =
        Kokkos::create_mirror_view( self.counts );

      Kokkos::deep_copy( host_counts , self.counts );

      int32_t max = 0 ;

      {
        const long n = host_counts.extent(0);
        for ( long i = 0 ; i < n ; ++i ) {
          if ( max < host_counts[i] ) max = host_counts[i] ;
        }
      }

      std::cout << "TestUniqueToken max reuse = " << max << std::endl ;

      typename view_type::HostMirror host_errors =
        Kokkos::create_mirror_view( self.errors );

      Kokkos::deep_copy( host_errors , self.errors );

      ASSERT_EQ( host_errors(0) , 0 );
    }
};


TEST_F( TEST_CATEGORY, unique_token )
{
  TestUniqueToken< TEST_EXECSPACE >::run();
}

} // namespace Test

