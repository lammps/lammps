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

#ifndef KOKKOS_CUDA_TEAM_HPP
#define KOKKOS_CUDA_TEAM_HPP

#include <iostream>
#include <algorithm>
#include <stdio.h>

#include <Kokkos_Macros.hpp>

/* only compile this file if CUDA is enabled for Kokkos */
#if defined( __CUDACC__ ) && defined( KOKKOS_ENABLE_CUDA )

#include <utility>
#include <Kokkos_Parallel.hpp>

#include <Cuda/Kokkos_CudaExec.hpp>
#include <Cuda/Kokkos_Cuda_ReduceScan.hpp>
#include <Cuda/Kokkos_Cuda_Internal.hpp>
#include <Kokkos_Vectorization.hpp>

#if defined(KOKKOS_ENABLE_PROFILING)
#include <impl/Kokkos_Profiling_Interface.hpp>
#include <typeinfo>
#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< typename Type >
struct CudaJoinFunctor {
  typedef Type value_type ;

  KOKKOS_INLINE_FUNCTION
  static void join( volatile value_type & update ,
                    volatile const value_type & input )
    { update += input ; }
};

/**\brief  Team member_type passed to TeamPolicy or TeamTask closures.
 *
 *  Cuda thread blocks for team closures are dimensioned as:
 *    blockDim.x == number of "vector lanes" per "thread"
 *    blockDim.y == number of "threads" per team
 *    blockDim.z == number of teams in a block
 *  where
 *    A set of teams exactly fill a warp OR a team is the whole block
 *      ( 0 == WarpSize % ( blockDim.x * blockDim.y ) )
 *      OR
 *      ( 1 == blockDim.z )
 *
 *  Thus when 1 < blockDim.z the team is warp-synchronous
 *  and __syncthreads should not be called in team collectives.
 *
 *  When multiple teams are mapped onto a single block then the
 *  total available shared memory must be partitioned among teams.
 */
class CudaTeamMember {
private:

  typedef Kokkos::Cuda                           execution_space ;
  typedef execution_space::scratch_memory_space  scratch_memory_space ;

  mutable void        * m_team_reduce ;
  scratch_memory_space  m_team_shared ;
  int                   m_team_reduce_size ;
  int                   m_league_rank ;
  int                   m_league_size ;

public:

  KOKKOS_INLINE_FUNCTION
  const execution_space::scratch_memory_space & team_shmem() const
    { return m_team_shared.set_team_thread_mode(0,1,0) ; }

  KOKKOS_INLINE_FUNCTION
  const execution_space::scratch_memory_space &
    team_scratch(const int& level) const
      { return m_team_shared.set_team_thread_mode(level,1,0) ; }

  KOKKOS_INLINE_FUNCTION
  const execution_space::scratch_memory_space &
    thread_scratch(const int& level) const
      { return m_team_shared.set_team_thread_mode(level,team_size(),team_rank()) ; }

  KOKKOS_INLINE_FUNCTION int league_rank() const { return m_league_rank ; }
  KOKKOS_INLINE_FUNCTION int league_size() const { return m_league_size ; }
  KOKKOS_INLINE_FUNCTION int team_rank() const
    {
      #ifdef __CUDA_ARCH__
      return threadIdx.y ;
      #else
      return 0;
      #endif
    }

  KOKKOS_INLINE_FUNCTION int team_size() const
    {
      #ifdef __CUDA_ARCH__
      return blockDim.y ;
      #else
      return 1;
      #endif
    }

  KOKKOS_INLINE_FUNCTION void team_barrier() const
    {
      #ifdef __CUDA_ARCH__
      if ( 1 == blockDim.z ) __syncthreads();       // team == block
      else                   __threadfence_block(); // team <= warp
      #endif
    }

  //--------------------------------------------------------------------------

  template<class ValueType>
  KOKKOS_INLINE_FUNCTION
  void team_broadcast( ValueType & val, const int& thread_id ) const
    {
      #ifdef __CUDA_ARCH__
      if ( 1 == blockDim.z ) { // team == block
        __syncthreads();
        // Wait for shared data write until all threads arrive here
        if ( threadIdx.x == 0u && threadIdx.y == (uint32_t)thread_id ) {
          *((ValueType*) m_team_reduce) = val ;
        }
        __syncthreads(); // Wait for shared data read until root thread writes
        val = *((ValueType*) m_team_reduce);
      }
      else { // team <= warp
        ValueType tmp( val ); // input might not be a register variable
        cuda_shfl( val, tmp, blockDim.x * thread_id, blockDim.x * blockDim.y );
      }
      #endif
    }
	
  template<class Closure, class ValueType>
  KOKKOS_INLINE_FUNCTION
  void team_broadcast( Closure const & f, ValueType & val, const int& thread_id ) const
    {
      #ifdef __CUDA_ARCH__
      f( val );

      if ( 1 == blockDim.z ) { // team == block
        __syncthreads();
        // Wait for shared data write until all threads arrive here
        if ( threadIdx.x == 0u && threadIdx.y == (uint32_t)thread_id ) {
          *((ValueType*) m_team_reduce) = val ;
        }
        __syncthreads(); // Wait for shared data read until root thread writes
        val = *((ValueType*) m_team_reduce);
      }
      else { // team <= warp
        ValueType tmp( val ); // input might not be a register variable
        cuda_shfl( val, tmp, blockDim.x * thread_id, blockDim.x * blockDim.y );
      }
      #endif
    }

  //--------------------------------------------------------------------------
  /**\brief  Reduction across a team
   *
   *  Mapping of teams onto blocks:
   *    blockDim.x  is "vector lanes"
   *    blockDim.y  is team "threads"
   *    blockDim.z  is number of teams per block
   *
   *  Requires:
   *    blockDim.x is power two
   *    blockDim.x <= CudaTraits::WarpSize
   *    ( 0 == CudaTraits::WarpSize % ( blockDim.x * blockDim.y )
   *      OR
   *    ( 1 == blockDim.z )
   */
  template< typename ReducerType >
  KOKKOS_INLINE_FUNCTION
  typename std::enable_if< is_reducer< ReducerType >::value >::type
  team_reduce( ReducerType const & reducer ) const noexcept
    {
      #ifdef __CUDA_ARCH__
      cuda_intra_block_reduction(reducer,blockDim.y);
      #endif /* #ifdef __CUDA_ARCH__ */
    }

  //--------------------------------------------------------------------------
  /** \brief  Intra-team exclusive prefix sum with team_rank() ordering
   *          with intra-team non-deterministic ordering accumulation.
   *
   *  The global inter-team accumulation value will, at the end of the
   *  league's parallel execution, be the scan's total.
   *  Parallel execution ordering of the league's teams is non-deterministic.
   *  As such the base value for each team's scan operation is similarly
   *  non-deterministic.
   */
  template< typename Type >
  KOKKOS_INLINE_FUNCTION
  Type team_scan( const Type & value , Type * const global_accum ) const
    {
      #ifdef __CUDA_ARCH__
      Type * const base_data = (Type *) m_team_reduce ;

      __syncthreads(); // Don't write in to shared data until all threads have entered this function

      if ( 0 == threadIdx.y ) { base_data[0] = 0 ; }

      base_data[ threadIdx.y + 1 ] = value ;

      Impl::cuda_intra_block_reduce_scan<true,Impl::CudaJoinFunctor<Type>,void>( Impl::CudaJoinFunctor<Type>() , base_data + 1 );

      if ( global_accum ) {
        if ( blockDim.y == threadIdx.y + 1 ) {
          base_data[ blockDim.y ] = atomic_fetch_add( global_accum , base_data[ blockDim.y ] );
        }
        __syncthreads(); // Wait for atomic
        base_data[ threadIdx.y ] += base_data[ blockDim.y ] ;
      }

      return base_data[ threadIdx.y ];
      #else
      return Type();
      #endif
    }

  /** \brief  Intra-team exclusive prefix sum with team_rank() ordering.
   *
   *  The highest rank thread can compute the reduction total as
   *    reduction_total = dev.team_scan( value ) + value ;
   */
  template< typename Type >
  KOKKOS_INLINE_FUNCTION Type team_scan( const Type & value ) const {
    return this->template team_scan<Type>( value , 0 );
  }

  //----------------------------------------

  template< typename ReducerType >
  KOKKOS_INLINE_FUNCTION static
  typename std::enable_if< is_reducer< ReducerType >::value >::type
  vector_reduce( ReducerType const & reducer )
    {

      #ifdef __CUDA_ARCH__
      if(blockDim.x == 1) return;

      // Intra vector lane shuffle reduction:
      typename ReducerType::value_type tmp ( reducer.reference() );

      unsigned mask = blockDim.x==32?0xffffffff:((1<<blockDim.x)-1)<<((threadIdx.y%(32/blockDim.x))*blockDim.x);

      for ( int i = blockDim.x ; ( i >>= 1 ) ; ) {
        cuda_shfl_down( reducer.reference() , tmp , i , blockDim.x , mask );
        if ( (int)threadIdx.x < i ) { reducer.join( tmp , reducer.reference() ); }
      }

      // Broadcast from root lane to all other lanes.
      // Cannot use "butterfly" algorithm to avoid the broadcast
      // because floating point summation is not associative
      // and thus different threads could have different results.

      cuda_shfl( reducer.reference() , tmp , 0 , blockDim.x , mask );
      #endif
    }

  //--------------------------------------------------------------------------
  /**\brief  Global reduction across all blocks
   *
   *  Return !0 if reducer contains the final value
   */
  template< typename ReducerType >
  KOKKOS_INLINE_FUNCTION static
  typename std::enable_if< is_reducer< ReducerType >::value , int >::type
  global_reduce( ReducerType const & reducer
               , int  * const global_scratch_flags
               , void * const global_scratch_space
               , void * const shmem
               , int    const shmem_size
               )
    {
    #ifdef __CUDA_ARCH__

      typedef typename ReducerType::value_type value_type ;
      typedef value_type volatile * pointer_type ;

      // Number of shared memory entries for the reduction:
      const int nsh = shmem_size / sizeof(value_type);

      // Number of CUDA threads in the block, rank within the block
      const int nid = blockDim.x * blockDim.y * blockDim.z ;
      const int tid = threadIdx.x + blockDim.x * (
                      threadIdx.y + blockDim.y * threadIdx.z );

      // Reduces within block using all available shared memory
      // Contributes if it is the root "vector lane"

      // wn == number of warps in the block
      // wx == which lane within the warp
      // wy == which warp within the block

      const int wn = ( nid + CudaTraits::WarpIndexMask ) >> CudaTraits::WarpIndexShift ;
      const int wx = tid &  CudaTraits::WarpIndexMask ;
      const int wy = tid >> CudaTraits::WarpIndexShift ;

      //------------------------
      { // Intra warp shuffle reduction from contributing CUDA threads

        value_type tmp( reducer.reference() );

        for ( int i = CudaTraits::WarpSize ; (int)blockDim.x <= ( i >>= 1 ) ; ) {

          cuda_shfl_down( reducer.reference(), tmp, i, CudaTraits::WarpSize );

          // Root of each vector lane reduces "thread" contribution
          if ( 0 == threadIdx.x && wx < i ) {
            reducer.join( & tmp , reducer.data() );
          }
        }

        // Reduce across warps using shared memory.
        // Number of warps may not be power of two.

        __syncthreads(); // Wait before shared data write

        // Number of shared memory entries for the reduction
        // is at most one per warp
        const int nentry = wn < nsh ? wn : nsh ;

        if ( 0 == wx && wy < nentry ) {
          // Root thread of warp 'wy' has warp's value to contribute
          ((value_type*) shmem)[wy] = tmp ;
        }

        __syncthreads(); // Wait for write to be visible to block

        // When more warps than shared entries
        // then warps must take turns joining their contribution
        // to the designated shared memory entry.
        for ( int i = nentry ; i < wn ; i += nentry ) {

          const int k = wy - i ;

          if ( 0 == wx && i <= wy && k < nentry ) {
            // Root thread of warp 'wy' has warp's value to contribute
            reducer.join( ((value_type*) shmem) + k , & tmp );
          }

          __syncthreads(); // Wait for write to be visible to block
        }

        // One warp performs the inter-warp reduction:

        if ( 0 == wy ) {

          // Start fan-in at power of two covering nentry

          for ( int i = ( 1 << ( 32 - __clz(nentry-1) ) ) ; ( i >>= 1 ) ; ) {
            const int k = wx + i ;
            if ( wx < i && k < nentry ) {
              reducer.join( ((pointer_type)shmem) + wx
                          , ((pointer_type)shmem) + k );
              __threadfence_block(); // Wait for write to be visible to warp
            }
          }
        }
      }
      //------------------------
      { // Write block's value to global_scratch_memory

        int last_block = 0 ;

        if ( 0 == wx ) {
          reducer.copy( ((pointer_type)global_scratch_space)
                        + blockIdx.x * reducer.length()
                      , reducer.data() );

          __threadfence(); // Wait until global write is visible.

          last_block = (int)gridDim.x ==
                       1 + Kokkos::atomic_fetch_add(global_scratch_flags,1);

          // If last block then reset count
          if ( last_block ) *global_scratch_flags = 0 ;
        }

        last_block = __syncthreads_or( last_block );

        if ( ! last_block ) return 0 ;

      }
      //------------------------
      // Last block reads global_scratch_memory into shared memory.

      const int nentry = nid < gridDim.x ?
                       ( nid       < nsh ? nid       : nsh ) :
                       ( gridDim.x < nsh ? gridDim.x : nsh ) ;

      // nentry = min( nid , nsh , gridDim.x )

      // whole block reads global memory into shared memory:

      if ( tid < nentry ) {

        const int offset = tid * reducer.length();

        reducer.copy( ((pointer_type)shmem) + offset
                    , ((pointer_type)global_scratch_space) + offset );

        for ( int i = nentry + tid ; i < (int)gridDim.x ; i += nentry ) {
          reducer.join( ((pointer_type)shmem) + offset
                      , ((pointer_type)global_scratch_space)
                        + i * reducer.length() );
        }
      }

      __syncthreads(); // Wait for writes to be visible to block

      if ( 0 == wy ) {

        // Iterate to reduce shared memory to single warp fan-in size

        const int nreduce = CudaTraits::WarpSize < nentry
                          ? CudaTraits::WarpSize : nentry ;

        // nreduce = min( CudaTraits::WarpSize , nsh , gridDim.x )

        if ( wx < nreduce && nreduce < nentry ) {
          for ( int i = nreduce + wx ; i < nentry ; i += nreduce ) {
            reducer.join( ((pointer_type)shmem) + wx
                        , ((pointer_type)shmem) + i );
          }
          __threadfence_block(); // Wait for writes to be visible to warp
        }

        // Start fan-in at power of two covering nentry

        for ( int i = ( 1 << ( 32 - __clz(nreduce-1) ) ) ; ( i >>= 1 ) ; ) {
          const int k = wx + i ;
          if ( wx < i && k < nreduce ) {
            reducer.join( ((pointer_type)shmem) + wx
                        , ((pointer_type)shmem) + k );
            __threadfence_block(); // Wait for writes to be visible to warp
          }
        }

        if ( 0 == wx ) {
          reducer.copy( reducer.data() , (pointer_type)shmem );
          return 1 ;
        }
      }
      return 0 ;

    #else
      return 0 ;
    #endif
    }

  //----------------------------------------
  // Private for the driver

  KOKKOS_INLINE_FUNCTION
  CudaTeamMember( void * shared
                , const int shared_begin
                , const int shared_size
                , void*     scratch_level_1_ptr
                , const int scratch_level_1_size
                , const int arg_league_rank
                , const int arg_league_size )
    : m_team_reduce( shared )
    , m_team_shared( ((char *)shared) + shared_begin , shared_size,  scratch_level_1_ptr, scratch_level_1_size)
    , m_team_reduce_size( shared_begin )
    , m_league_rank( arg_league_rank )
    , m_league_size( arg_league_size )
    {}

public:
  // Declare to avoid unused private member warnings which are trigger
  // when SFINAE excludes the member function which uses these variables
  // Making another class a friend also surpresses these warnings
  bool impl_avoid_sfinae_warning() const noexcept
  {
    return m_team_reduce_size > 0 && m_team_reduce != nullptr;
  }
};

} // namspace Impl
} // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template<typename iType>
struct TeamThreadRangeBoundariesStruct<iType,CudaTeamMember> {
  typedef iType index_type;
  const CudaTeamMember& member;
  const iType start;
  const iType end;

  KOKKOS_INLINE_FUNCTION
  TeamThreadRangeBoundariesStruct (const CudaTeamMember& thread_, const iType& count)
    : member(thread_)
    , start( 0 )
    , end( count ) {}

  KOKKOS_INLINE_FUNCTION
  TeamThreadRangeBoundariesStruct (const CudaTeamMember& thread_,  const iType& begin_, const iType& end_)
    : member(thread_)
    , start( begin_ )
    , end( end_ ) {}
};



template<typename iType>
struct ThreadVectorRangeBoundariesStruct<iType,CudaTeamMember> {
  typedef iType index_type;
  const index_type start;
  const index_type end;

  KOKKOS_INLINE_FUNCTION
  ThreadVectorRangeBoundariesStruct (const CudaTeamMember, const index_type& count)
    : start( static_cast<index_type>(0) ), end( count ) {}

  KOKKOS_INLINE_FUNCTION
  ThreadVectorRangeBoundariesStruct (const index_type& count)
    : start( static_cast<index_type>(0) ), end( count ) {}

  KOKKOS_INLINE_FUNCTION
  ThreadVectorRangeBoundariesStruct (const CudaTeamMember, const index_type& arg_begin, const index_type& arg_end)
    : start( arg_begin ), end( arg_end ) {}

  KOKKOS_INLINE_FUNCTION
  ThreadVectorRangeBoundariesStruct (const index_type& arg_begin, const index_type& arg_end)
    : start( arg_begin ), end( arg_end ) {}
};

} // namespace Impl

template<typename iType>
KOKKOS_INLINE_FUNCTION
Impl::TeamThreadRangeBoundariesStruct< iType, Impl::CudaTeamMember >
TeamThreadRange( const Impl::CudaTeamMember & thread, const iType & count ) {
  return Impl::TeamThreadRangeBoundariesStruct< iType, Impl::CudaTeamMember >( thread, count );
}

template< typename iType1, typename iType2 >
KOKKOS_INLINE_FUNCTION
Impl::TeamThreadRangeBoundariesStruct< typename std::common_type< iType1, iType2 >::type,
                                       Impl::CudaTeamMember >
TeamThreadRange( const Impl::CudaTeamMember & thread, const iType1 & begin, const iType2 & end ) {
  typedef typename std::common_type< iType1, iType2 >::type iType;
  return Impl::TeamThreadRangeBoundariesStruct< iType, Impl::CudaTeamMember >( thread, iType(begin), iType(end) );
}

template<typename iType>
KOKKOS_INLINE_FUNCTION
Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::CudaTeamMember >
ThreadVectorRange(const Impl::CudaTeamMember& thread, const iType& count) {
  return Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::CudaTeamMember >(thread,count);
}

template<typename iType>
KOKKOS_INLINE_FUNCTION
Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::CudaTeamMember >
ThreadVectorRange(const Impl::CudaTeamMember& thread, const iType& arg_begin, const iType& arg_end) {
  return Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::CudaTeamMember >(thread,arg_begin,arg_end);
}

KOKKOS_INLINE_FUNCTION
Impl::ThreadSingleStruct<Impl::CudaTeamMember> PerTeam(const Impl::CudaTeamMember& thread) {
  return Impl::ThreadSingleStruct<Impl::CudaTeamMember>(thread);
}

KOKKOS_INLINE_FUNCTION
Impl::VectorSingleStruct<Impl::CudaTeamMember> PerThread(const Impl::CudaTeamMember& thread) {
  return Impl::VectorSingleStruct<Impl::CudaTeamMember>(thread);
}

//----------------------------------------------------------------------------

/** \brief  Inter-thread parallel_for.
 *
 *  Executes closure(iType i) for each i=[0..N).
 *
 * The range [0..N) is mapped to all threads of the the calling thread team.
 */
template<typename iType, class Closure >
KOKKOS_INLINE_FUNCTION
void parallel_for
  ( const Impl::TeamThreadRangeBoundariesStruct<iType,Impl::CudaTeamMember>&
      loop_boundaries
  , const Closure & closure
  )
{
  #ifdef __CUDA_ARCH__
  for( iType i = loop_boundaries.start + threadIdx.y
     ; i < loop_boundaries.end
     ; i += blockDim.y )
    closure(i);
  #endif
}

//----------------------------------------------------------------------------

/** \brief  Inter-thread parallel_reduce with a reducer.
 *
 *  Executes closure(iType i, ValueType & val) for each i=[0..N)
 *
 *  The range [0..N) is mapped to all threads of the
 *  calling thread team and a summation of val is
 *  performed and put into result.
 */
template< typename iType, class Closure, class ReducerType >
KOKKOS_INLINE_FUNCTION
typename std::enable_if< Kokkos::is_reducer< ReducerType >::value >::type
parallel_reduce
  ( const Impl::TeamThreadRangeBoundariesStruct<iType,Impl::CudaTeamMember> &
      loop_boundaries
  , const Closure & closure
  , const ReducerType & reducer
  )
{
#ifdef __CUDA_ARCH__

  reducer.init( reducer.reference() );

  for( iType i = loop_boundaries.start + threadIdx.y
     ; i < loop_boundaries.end
     ; i += blockDim.y ) {
    closure(i,reducer.reference());
  }

  loop_boundaries.member.team_reduce( reducer );

#endif
}


/** \brief  Inter-thread parallel_reduce assuming summation.
 *
 *  Executes closure(iType i, ValueType & val) for each i=[0..N)
 *
 *  The range [0..N) is mapped to all threads of the
 *  calling thread team and a summation of val is
 *  performed and put into result.
 */
template< typename iType, class Closure, typename ValueType >
KOKKOS_INLINE_FUNCTION
typename std::enable_if< ! Kokkos::is_reducer< ValueType >::value >::type
parallel_reduce
  ( const Impl::TeamThreadRangeBoundariesStruct<iType,Impl::CudaTeamMember> &
      loop_boundaries
  , const Closure & closure
  , ValueType & result
  )
{
#ifdef __CUDA_ARCH__

  Kokkos::Sum<ValueType> reducer(result);

  reducer.init( reducer.reference() );

  for( iType i = loop_boundaries.start + threadIdx.y
     ; i < loop_boundaries.end
     ; i += blockDim.y ) {
    closure(i,result);
  }

  loop_boundaries.member.team_reduce( reducer );

#endif
}

//----------------------------------------------------------------------------

/** \brief  Intra-thread vector parallel_for.
 *
 *  Executes closure(iType i) for each i=[0..N)
 *
 * The range [0..N) is mapped to all vector lanes of the the calling thread.
 */
template<typename iType, class Closure >
KOKKOS_INLINE_FUNCTION
void parallel_for
  ( const Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::CudaTeamMember>&
    loop_boundaries
  , const Closure & closure
  )
{
#ifdef __CUDA_ARCH__
  for ( iType i = loop_boundaries.start + threadIdx.x
      ; i < loop_boundaries.end
      ; i += blockDim.x ) {
    closure(i);
  }
  #ifdef KOKKOS_IMPL_CUDA_SYNCWARP_NEEDS_MASK
  KOKKOS_IMPL_CUDA_SYNCWARP_MASK(blockDim.x==32?0xffffffff:((1<<blockDim.x)-1)<<(threadIdx.y%(32/blockDim.x))*blockDim.x);
  #else
  KOKKOS_IMPL_CUDA_SYNCWARP;
  #endif
#endif
}

//----------------------------------------------------------------------------

/** \brief  Intra-thread vector parallel_reduce.
 *
 *  Calls closure(iType i, ValueType & val) for each i=[0..N).
 *
 *  The range [0..N) is mapped to all vector lanes of
 *  the calling thread and a reduction of val is performed using +=
 *  and output into result.
 *
 *  The identity value for the += operator is assumed to be the default
 *  constructed value.
 */
template< typename iType, class Closure, class ReducerType >
KOKKOS_INLINE_FUNCTION
typename std::enable_if< is_reducer< ReducerType >::value >::type
parallel_reduce
  ( Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::CudaTeamMember>
      const & loop_boundaries
  , Closure const & closure
  , ReducerType const & reducer )
{
#ifdef __CUDA_ARCH__

  reducer.init( reducer.reference() );

  for ( iType i = loop_boundaries.start + threadIdx.x
      ; i < loop_boundaries.end
      ; i += blockDim.x ) {
    closure(i,reducer.reference());
  }

  Impl::CudaTeamMember::vector_reduce( reducer );

#endif
}

/** \brief  Intra-thread vector parallel_reduce.
 *
 *  Calls closure(iType i, ValueType & val) for each i=[0..N).
 *
 *  The range [0..N) is mapped to all vector lanes of
 *  the calling thread and a reduction of val is performed using +=
 *  and output into result.
 *
 *  The identity value for the += operator is assumed to be the default
 *  constructed value.
 */
template< typename iType, class Closure, typename ValueType >
KOKKOS_INLINE_FUNCTION
typename std::enable_if< ! is_reducer< ValueType >::value >::type
parallel_reduce
  ( Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::CudaTeamMember>
      const & loop_boundaries
  , Closure const & closure
  , ValueType & result )
{
#ifdef __CUDA_ARCH__
  result = ValueType();

  for ( iType i = loop_boundaries.start + threadIdx.x
      ; i < loop_boundaries.end
      ; i += blockDim.x ) {
    closure(i,result);
  }

  Impl::CudaTeamMember::vector_reduce(
    Kokkos::Sum<ValueType>(result ) );

#endif
}

//----------------------------------------------------------------------------

/** \brief  Intra-thread vector parallel exclusive prefix sum.
 *
 *  Executes closure(iType i, ValueType & val, bool final) for each i=[0..N)
 *
 *  The range [0..N) is mapped to all vector lanes in the
 *  thread and a scan operation is performed.
 *  The last call to closure has final == true.
 */
template< typename iType, class Closure >
KOKKOS_INLINE_FUNCTION
void parallel_scan
  ( const Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::CudaTeamMember >&
      loop_boundaries
  , const Closure & closure
  )
{

#ifdef __CUDA_ARCH__

  // Extract value_type from closure

  using value_type =
    typename Kokkos::Impl::FunctorAnalysis
      < Kokkos::Impl::FunctorPatternInterface::SCAN
      , void
      , Closure >::value_type ;

  // Loop through boundaries by vector-length chunks
  // must scan at each iteration

  value_type accum = 0 ;

  // All thread "lanes" must loop the same number of times.
  // Determine an loop end for all thread "lanes."
  // Requires:
  //   blockDim.x is power of two and thus
  //     ( end % blockDim.x ) == ( end & ( blockDim.x - 1 ) )
  //   1 <= blockDim.x <= CudaTraits::WarpSize

  const int mask = blockDim.x - 1 ;
  const unsigned active_mask = blockDim.x==32?0xffffffff:((1<<blockDim.x)-1)<<(threadIdx.y%(32/blockDim.x))*blockDim.x;
  const int rem  = loop_boundaries.end & mask ; // == end % blockDim.x
  const int end  = loop_boundaries.end + ( rem ? blockDim.x - rem : 0 );

  for ( int i = threadIdx.x ; i < end ; i += blockDim.x ) {

    value_type val = 0 ;

    // First acquire per-lane contributions:
    if ( i < loop_boundaries.end ) closure( i , val , false );

    value_type sval = val ;

    // Bottom up inclusive scan in triangular pattern
    // where each CUDA thread is the root of a reduction tree
    // from the zeroth "lane" to itself.
    //  [t] += [t-1] if t >= 1
    //  [t] += [t-2] if t >= 2
    //  [t] += [t-4] if t >= 4
    //  ...

    for ( int j = 1 ; j < (int)blockDim.x ; j <<= 1 ) {
      value_type tmp = 0 ;
      Impl::cuda_shfl_up(tmp, sval , j , blockDim.x, active_mask );
      if ( j <= (int)threadIdx.x ) { sval += tmp ; }
    }

    // Include accumulation and remove value for exclusive scan:
    val = accum + sval - val ;

    // Provide exclusive scan value:
    if ( i < loop_boundaries.end ) closure( i , val , true );

    // Accumulate the last value in the inclusive scan:
    Impl::cuda_shfl( sval , sval , mask , blockDim.x, active_mask );

    accum += sval ;
  }

#endif
}

}

namespace Kokkos {

template<class FunctorType>
KOKKOS_INLINE_FUNCTION
void single(const Impl::VectorSingleStruct<Impl::CudaTeamMember>& , const FunctorType& lambda) {
#ifdef __CUDA_ARCH__
  if(threadIdx.x == 0) lambda();
  #ifdef KOKKOS_IMPL_CUDA_SYNCWARP_NEEDS_MASK
  KOKKOS_IMPL_CUDA_SYNCWARP_MASK(blockDim.x==32?0xffffffff:((1<<blockDim.x)-1)<<(threadIdx.y%(32/blockDim.x))*blockDim.x);
  #else
  KOKKOS_IMPL_CUDA_SYNCWARP;
  #endif
#endif
}

template<class FunctorType>
KOKKOS_INLINE_FUNCTION
void single(const Impl::ThreadSingleStruct<Impl::CudaTeamMember>& , const FunctorType& lambda) {
#ifdef __CUDA_ARCH__
  if(threadIdx.x == 0 && threadIdx.y == 0) lambda();
  #ifdef KOKKOS_IMPL_CUDA_SYNCWARP_NEEDS_MASK
  KOKKOS_IMPL_CUDA_SYNCWARP_MASK(blockDim.x==32?0xffffffff:((1<<blockDim.x)-1)<<(threadIdx.y%(32/blockDim.x))*blockDim.x);
  #else
  KOKKOS_IMPL_CUDA_SYNCWARP;
  #endif
#endif
}

template<class FunctorType, class ValueType>
KOKKOS_INLINE_FUNCTION
void single(const Impl::VectorSingleStruct<Impl::CudaTeamMember>& , const FunctorType& lambda, ValueType& val) {
#ifdef __CUDA_ARCH__
  if(threadIdx.x == 0) lambda(val);
  unsigned mask = blockDim.x==32?0xffffffff:((1<<blockDim.x)-1)<<((threadIdx.y%(32/blockDim.x))*blockDim.x);
  Impl::cuda_shfl(val,val,0,blockDim.x,mask);
#endif
}

template<class FunctorType, class ValueType>
KOKKOS_INLINE_FUNCTION
void single(const Impl::ThreadSingleStruct<Impl::CudaTeamMember>& single_struct, const FunctorType& lambda, ValueType& val) {
#ifdef __CUDA_ARCH__
  if(threadIdx.x == 0 && threadIdx.y == 0) {
    lambda(val);
  }
  single_struct.team_member.team_broadcast(val,0);
#endif
}

} // namespace Kokkos

#endif /* defined( __CUDACC__ ) */

#endif /* #ifndef KOKKOS_CUDA_TEAM_HPP */

