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

#ifndef EXAMPLE_GROW_ARRAY
#define EXAMPLE_GROW_ARRAY

#include <cstdlib>

#include <Kokkos_Core.hpp>

#include <algorithm>

#if defined(KOKKOS_ENABLE_CUDA)
#include <thrust/device_ptr.h>
#include <thrust/sort.h>
#endif

namespace Example {

//----------------------------------------------------------------------------

template< class ExecSpace >
struct SortView {

  template< typename ValueType >
  SortView( const Kokkos::View<ValueType*,ExecSpace> v , int begin , int end )
    {
      std::sort( v.ptr_on_device() + begin , v.ptr_on_device() + end );
    }
};

#if defined(KOKKOS_ENABLE_CUDA)
template<>
struct SortView< Kokkos::Cuda > {
  template< typename ValueType >
  SortView( const Kokkos::View<ValueType*,Kokkos::Cuda> v , int begin , int end )
    {
      thrust::sort( thrust::device_ptr<ValueType>( v.ptr_on_device() + begin )
                  , thrust::device_ptr<ValueType>( v.ptr_on_device() + end ) );
    }
};
#endif



//----------------------------------------------------------------------------

template< class ExecSpace >
struct GrowArrayFunctor {

  typedef ExecSpace  execution_space ;

  enum { SHIFT = sizeof(int) == 8 ? 6 : 5 }; // 8 or 4 byte int
  enum { MASK  = ( 1 << SHIFT ) - 1 };

  const Kokkos::View<int*,ExecSpace>  m_search_flags ; // bit flags for values to append
  const Kokkos::View<int*,ExecSpace>  m_search_array ; // array to append values
  const Kokkos::View<int,ExecSpace>   m_search_count ; // offset
  const int m_search_total ;
  const int m_search_team_chunk ;

  GrowArrayFunctor( int array_length , int search_length , int print = 1 )
    : m_search_flags( "flags" , ( search_length + MASK ) >> SHIFT ) // One bit per search entry
    , m_search_array( "array" , array_length )
    , m_search_count( "count" )
    , m_search_total( search_length )
    , m_search_team_chunk( 2048 )
    {}

  KOKKOS_INLINE_FUNCTION
  bool flag_is_set( const int index ) const
    {
      // 64 or 32 bit integer:

      const int j = index >> SHIFT ; // which integer flag
      const int k = 1 << ( index & MASK ); // which bit in that integer
      const int s = ( j < int(m_search_flags.dimension_0()) ) && ( 0 != ( m_search_flags(j) & k ) );

      return s ;
    }

  typedef typename Kokkos::TeamPolicy<ExecSpace>::member_type team_member ;

  KOKKOS_INLINE_FUNCTION
  void operator()( const team_member & member ) const
    {
      enum { LOCAL_BUFFER_LENGTH = 16 };

      int local_buffer[ LOCAL_BUFFER_LENGTH ] ;
      int local_count = 0 ;

      // Each team searches 'm_search_team_chunk' indices.
      // The threads of a team must iterate together because all
      // threads in the team must call 'team_scan' to prevent deadlock in the team.

            int search_team_begin = member.league_rank() * m_search_team_chunk ;
      const int search_team_end   = search_team_begin + m_search_team_chunk ;

      int k = 0 ;

      while ( search_team_begin < search_team_end ) {

        // This iteration searches [ search_team_begin .. search_team_begin + member.team_size() ]
        const int thread_search_index = search_team_begin + member.team_rank();

        // If this thread's search index is in the range
        // and the flag is set, push into this thread's local buffer.
        if ( thread_search_index < m_search_total && flag_is_set(thread_search_index) ) {
          local_buffer[ local_count ] = thread_search_index ;
          ++local_count ;
        }

        // Move the team's search range forward
        search_team_begin += member.team_size(); // Striding team by team size

        // Count number of times a thread's buffer might have grown:
        ++k ;

        // Write buffer if end of search or a thread might have filled its buffer.
        if ( k == LOCAL_BUFFER_LENGTH /* A thread in my team might have filled its buffer */ ||
             ! ( search_team_begin < search_team_end ) /* Team is at the end of its search */ ) {

          // Team's exclusive scan of threads' contributions, with global offset.
          // This thread writes its buffer into [ team_offset .. team_offset + local_count )
          const int team_offset = member.team_scan( local_count , & *m_search_count );

          // Copy locally buffered entries into global array:
          for ( int i = 0 ; i < local_count ; ++i ) {
            m_search_array( team_offset + i ) = local_buffer[i] ;
          }

          k = 0 ;
          local_count = 0 ;
        }
      }
    }
};


template< class ExecSpace >
void grow_array( int array_length , int search_length , int print = 1 )
{
  typedef GrowArrayFunctor< ExecSpace > FunctorType ;

  FunctorType functor( array_length , search_length , print );

  typename Kokkos::View<int,ExecSpace>::HostMirror  count = Kokkos::create_mirror_view( functor.m_search_count );
  typename Kokkos::View<int*,ExecSpace>::HostMirror flags = Kokkos::create_mirror_view( functor.m_search_flags );

  // Set at most 'array_length' random bits over the search length.
  for ( int i = 0 ; i < array_length ; ++i ) {
    // 'lrand48()' generates random number between [0..2^31]
    // index = ( lrand48() * search_length ) / ( 2^31 )
    const long int index = ( lrand48() * search_length ) >> 31 ;
    // set the bit within the flags:
    flags( index >> FunctorType::SHIFT ) |= ( 1 << ( index & FunctorType::MASK ) );
  }

  Kokkos::deep_copy( functor.m_search_flags , flags );

  // Each team works on 'functor.m_search_team_chunk' span of the search_length
  Kokkos::TeamPolicy< ExecSpace >
    work( /* #teams */ ( search_length + functor.m_search_team_chunk - 1 ) / functor.m_search_team_chunk
        , /* threads/team */ Kokkos::TeamPolicy< ExecSpace >::team_size_max( functor ) );

  // Fill array:
  Kokkos::parallel_for( work , functor );

  // How much was filled:
  Kokkos::deep_copy( count , functor.m_search_count );

  // Sort array:
  SortView< ExecSpace >( functor.m_search_array , 0 , *count );

  // Mirror the results:
  typename Kokkos::View<int*,ExecSpace>::HostMirror results = Kokkos::create_mirror_view( functor.m_search_array );
  Kokkos::deep_copy( results , functor.m_search_array );

  // Verify results:
  int result_error_count = 0 ;
  int flags_error_count = 0 ;
  for ( int i = 0 ; i < *count ; ++i ) {
    const int index = results(i);
    const int entry = index >> FunctorType::SHIFT ;
    const int bit   = 1 << ( index & FunctorType::MASK );
    const bool flag = 0 != ( flags( entry ) & bit );
    if ( ! flag ) {
      if ( print ) std::cerr << "result( " << i << " : " << index << " )";
      ++result_error_count ;
    }
    flags( entry ) &= ~bit ; // Clear that verified bit
  }

  for ( int i = 0 ; i < int(flags.dimension_0()) ; ++i ) {
    // If any uncleared bits then an error
    if ( flags(i) ) {
      if ( print ) std::cerr << "flags( " << i << " : " << flags(i) << " )" ;
      ++flags_error_count ;
    }
  }

  if ( result_error_count || flags_error_count ) {
    std::cerr << std::endl << "Example::GrowArrayFunctor( " << array_length
              << " , " << search_length
              << " ) result_error_count( " << result_error_count << " )"
              << " ) flags_error_count( " << flags_error_count << " )"
              << std::endl ;
  }
}


} // namespace Example

//----------------------------------------------------------------------------

#endif /* #ifndef EXAMPLE_GROW_ARRAY */

