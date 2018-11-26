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
#include <cstdio>

namespace Test {

template< class Device, class WorkSpec = size_t >
struct TestScan {
  typedef  Device    execution_space;
  typedef  long int  value_type;

  Kokkos::View< int, Device, Kokkos::MemoryTraits<Kokkos::Atomic> > errors;

  KOKKOS_INLINE_FUNCTION
  void operator()( const int iwork, value_type & update, const bool final_pass ) const
  {
    const value_type n = iwork + 1;
    const value_type imbalance = ( ( 1000 <= n ) && ( 0 == n % 1000 ) ) ? 1000 : 0;

    // Insert an artificial load imbalance

    for ( value_type i = 0; i < imbalance; ++i ) { ++update; }

    update += n - imbalance;

    if ( final_pass ) {
      const value_type answer = n & 1 ? ( n * ( ( n + 1 ) / 2 ) ) : ( ( n / 2 ) * ( n + 1 ) );

      if ( answer != update ) {
        int fail = errors()++;

        if ( fail < 20 ) {
          printf( "TestScan(%d,%ld) != %ld\n", iwork, update, answer );
        }
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void init( value_type & update ) const { update = 0; }

  KOKKOS_INLINE_FUNCTION
  void join( volatile       value_type & update,
             volatile const value_type & input ) const
  { update += input; }

  TestScan( const WorkSpec & N )
  {
    Kokkos::View< int, Device > errors_a( "Errors" );
    Kokkos::deep_copy( errors_a, 0 );
    errors = errors_a;

    Kokkos::parallel_scan( N , *this );

    long long int total = 0;
    Kokkos::parallel_scan( N, *this, total );
    run_check( size_t( ( N+1 )*N/2 ), size_t( total ) );
    check_error();
  }

  TestScan( const WorkSpec & Start , const WorkSpec & N )
  {
    typedef Kokkos::RangePolicy< execution_space > exec_policy ;

    Kokkos::View< int, Device > errors_a( "Errors" );
    Kokkos::deep_copy( errors_a, 0 );
    errors = errors_a;
    
    Kokkos::parallel_scan( exec_policy( Start , N ) , *this );
    check_error();
  }

  void check_error() {
    int total_errors;
    Kokkos::deep_copy(total_errors, errors);
    ASSERT_EQ(total_errors,0);
  }

  static void test_range( const WorkSpec & begin, const WorkSpec & end )
  {
    for ( WorkSpec i = begin; i < end; ++i ) {
      (void) TestScan( i );
    }
  }

  void run_check( const size_t & expected, const size_t & actual )
  { 
    ASSERT_EQ( expected, actual ); 
  }

};

TEST_F( TEST_CATEGORY, scan )
{
  TestScan< TEST_EXECSPACE >::test_range( 1, 1000 );
  TestScan< TEST_EXECSPACE >( 0 );
  TestScan< TEST_EXECSPACE >( 100000 );
  TestScan< TEST_EXECSPACE >( 10000000 );
  TEST_EXECSPACE::fence();
}


/*TEST_F( TEST_CATEGORY, scan_small )
{
  typedef TestScan< TEST_EXECSPACE, Kokkos::Impl::ThreadsExecUseScanSmall > TestScanFunctor;

  for ( int i = 0; i < 1000; ++i ) {
    TestScanFunctor( 10 );
    TestScanFunctor( 10000 );
  }
  TestScanFunctor( 1000000 );
  TestScanFunctor( 10000000 );

  TEST_EXECSPACE::fence();
}*/


} // namespace Test
