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

#ifndef KOKKOS_IMPL_ROCM_TASK_HPP
#define KOKKOS_IMPL_ROCM_TASK_HPP

#if defined( KOKKOS_ENABLE_TASKDAG )

#include <ROCm/Kokkos_ROCm_Vectorization.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template< class > class TaskExec ; 

template<>
class TaskQueueSpecialization< Kokkos::Experimental::ROCm >
{
public:

  using execution_space = Kokkos::Experimental::ROCm ;
  using queue_type      = Kokkos::Impl::TaskQueue< execution_space > ;
  using task_base_type  = Kokkos::Impl::TaskBase< execution_space , void , void > ;
  using member_type     = TaskExec< execution_space > ;

  // Must specify memory space
  using memory_space = Kokkos::HostSpace ;

  static
  void iff_single_thread_recursive_execute( queue_type * const ) {}

  KOKKOS_INLINE_FUNCTION
  static void driver( queue_type * const, hc::tiled_index<3> );

  // Must provide task queue execution function
  static void execute( queue_type * const );

  // Must provide mechanism to set function pointer in
  // execution space from the host process.
  template< typename FunctorType >
  static
  void proc_set_apply( typename TaskBase< Kokkos::Experimental::ROCm
                               , typename FunctorType::value_type
                               , FunctorType
                               >::function_type * ptr )
    {
      using TaskType = TaskBase< Kokkos::Experimental::ROCm
                               , typename FunctorType::value_type
                               , FunctorType
                               > ;
      hc::extent< 1 > flat_extent( 1 );
      hc::tiled_extent< 1 > team_extent = flat_extent.tile( 1);

      hc::parallel_for_each( team_extent , [&](hc::tiled_index<1> idx) [[hc]]
      {
         *ptr = TaskType::apply ;
      }).wait();
    }
};

/*template<>
KOKKOS_FUNCTION 
void TaskQueue<Kokkos::Experimental::ROCm>::decrement( typename TaskQueue<Kokkos::Experimental::ROCm>::task_root_type *
) {}
*/
extern template class TaskQueue< Kokkos::Experimental::ROCm > ;

//----------------------------------------------------------------------------
/**\brief  Impl::TaskExec<ROCm> is the TaskScheduler<ROCm>::member_type
 *         passed to tasks running in a ROCm space.
 *
 *  ROCm thread blocks for tasking are dimensioned:
 *    idx.tile_dim[0] == vector length
 *    idx.tile_dim[1] == team size
 *    idx.tile_dim[2] == number of teams
 *  where
 *    idx.tile_dim[0] * idx.tile_dim[1] == WavefrontSize
 *
 *  Both single thread and thread team tasks are run by a full ROCm warp.
 *  A single thread task is called by warp lane #0 and the remaining
 *  lanes of the warp are idle.
 */
template<>
class TaskExec< Kokkos::Experimental::ROCm >
{
private:

  TaskExec( TaskExec && ) = delete ;
  TaskExec( TaskExec const & ) = delete ;
  TaskExec & operator = ( TaskExec && ) = delete ;
  TaskExec & operator = ( TaskExec const & ) = delete ;


  friend class Kokkos::Impl::TaskQueue< Kokkos::Experimental::ROCm > ;
  friend class Kokkos::Impl::TaskQueueSpecialization< Kokkos::Experimental::ROCm > ;

  int              m_team_size ;
  hc::tiled_index<3>      m_idx;

//  KOKKOS_INLINE_FUNCTION TaskExec( int arg_team_size )  //TODO: tile_dim[0]
//    : m_team_size( arg_team_size ) {}

  KOKKOS_INLINE_FUNCTION TaskExec( int arg_team_size,
                                   hc::tiled_index<3> tidx)  
    : m_team_size( arg_team_size),
      m_idx( tidx ) {}

public:
//      const auto local = t_idx.local[0];
//      const auto global = t_idx.global[0];
//     const auto tile = t_idx.tile[0];

  hc::tiled_index<3> idx() const { return m_idx;}

#if defined( __HCC_ACCELERATOR__ )
  KOKKOS_INLINE_FUNCTION void team_barrier() { /* __threadfence_block(); */ }
  KOKKOS_INLINE_FUNCTION int  team_rank() const { return m_idx.local[1] ; } // t_idx.tile[0];
  KOKKOS_INLINE_FUNCTION int  team_size() const { return m_team_size ; }
#else
  KOKKOS_INLINE_FUNCTION void team_barrier() {}
  KOKKOS_INLINE_FUNCTION int  team_rank() const { return 0 ; }
  KOKKOS_INLINE_FUNCTION int  team_size() const { return 0 ; }
#endif
};

}} /* namespace Kokkos::Impl */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
namespace Kokkos {

template<typename iType>
KOKKOS_INLINE_FUNCTION
Impl::TeamThreadRangeBoundariesStruct<iType,Impl::TaskExec< Kokkos::Experimental::ROCm > >
TeamThreadRange
  ( Impl::TaskExec< Kokkos::Experimental::ROCm > & thread, const iType & count )
{
  return Impl::TeamThreadRangeBoundariesStruct<iType,Impl::TaskExec< Kokkos::Experimental::ROCm > >(thread,count);
}

template<typename iType1, typename iType2>
KOKKOS_INLINE_FUNCTION
Impl::TeamThreadRangeBoundariesStruct< typename std::common_type< iType1, iType2 >::type,
                                       Impl::TaskExec< Kokkos::Experimental::ROCm > >
TeamThreadRange
  ( Impl:: TaskExec< Kokkos::Experimental::ROCm > & thread, const iType1 & begin, const iType2 & end )
{
  typedef typename std::common_type<iType1, iType2>::type iType;
  return Impl::TeamThreadRangeBoundariesStruct<iType, Impl::TaskExec< Kokkos::Experimental::ROCm > >(thread, begin, end);
}

template<typename iType>
KOKKOS_INLINE_FUNCTION
Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::TaskExec< Kokkos::Experimental::ROCm > >
ThreadVectorRange
  ( Impl::TaskExec< Kokkos::Experimental::ROCm > & thread
  , const iType & count )
{
  return Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::TaskExec< Kokkos::Experimental::ROCm > >(thread,count);
}

/** \brief  Inter-thread parallel_for. Executes lambda(iType i) for each i=0..N-1.
 *
 * The range i=0..N-1 is mapped to all threads of the the calling thread team.
 * This functionality requires C++11 support.
*/
template<typename iType, class Lambda>
KOKKOS_INLINE_FUNCTION
void parallel_for
  ( const Impl::TeamThreadRangeBoundariesStruct<iType,Impl:: TaskExec< Kokkos::Experimental::ROCm > >& loop_boundaries
  , const Lambda& lambda
  )
{
  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment) {
    lambda(i);
  }
}

// reduce across corresponding lanes between team members within workgroup
// assume stride*team_size == workgroup_size
template< typename ValueType >
KOKKOS_INLINE_FUNCTION
void strided_shfl_workgroup_reduction
  (const ValueType& f(),
   ValueType& val,
   int team_size,
   int stride)
{
  for (int lane_delta=(team_size*stride)>>1; lane_delta>=stride; lane_delta>>=1) {
    f(val, Kokkos::shfl_down(val, lane_delta, team_size*stride));
  }
}

template< typename ValueType, class JoinType >
KOKKOS_INLINE_FUNCTION
void strided_shfl_workgroup_reduction
  (const JoinType& join,
   ValueType& val,
   int team_size,
   int stride)
{
  for (int lane_delta=(team_size*stride)>>1; lane_delta>=stride; lane_delta>>=1) {
    join(val, shfl_down(val, lane_delta, team_size*stride));
  }
}

// multiple within-workgroup non-strided reductions
template< typename ValueType, class JoinType >
KOKKOS_INLINE_FUNCTION
void multi_shfl_workgroup_reduction
  (const JoinType& join,
   ValueType& val,
   int vec_length)
{
  for (int lane_delta=vec_length>>1; lane_delta; lane_delta>>=1) {
    join(val, shfl_down(val, lane_delta, vec_length));
  }
}

// broadcast within workgroup
template< class ValueType >
KOKKOS_INLINE_FUNCTION
ValueType shfl_workgroup_broadcast
  (ValueType& val,
   int src_lane,
   int width)
{
  return shfl(val, src_lane, width);
}

// all-reduce across corresponding vector lanes between team members within workgroup
// assume vec_length*team_size == workgroup_size
// blockDim.x == vec_length == stride
// blockDim.y == team_size
// threadIdx.x == position in vec
// threadIdx.y == member number

template<typename iType, class Lambda, typename ValueType>
KOKKOS_INLINE_FUNCTION
void parallel_reduce
  ( const Impl::TeamThreadRangeBoundariesStruct<iType,Impl:: TaskExec< Kokkos::Experimental::ROCm > >& loop_boundaries
  , const Lambda& lambda
  , ValueType& initialized_result)
{
  int team_rank = loop_boundaries.thread.team_rank(); // member num within the team
  ValueType result = initialized_result;
  hc::tiled_index<3> idx = loop_boundaries.thread.idx();

  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment) {
    lambda(i, result);
  }
  initialized_result = result;

  strided_shfl_workgroup_reduction(
                          [&] (ValueType& val1, const ValueType& val2) { val1 += val2; },
                          initialized_result,
                          loop_boundaries.thread.team_size(),
                          idx.tile_dim[0]);
  initialized_result = shfl_workgroup_broadcast<ValueType>( initialized_result, idx.local[0], Impl::ROCmTraits::WavefrontSize );

}

template< typename iType, class Lambda, typename ValueType, class JoinType >
KOKKOS_INLINE_FUNCTION
void parallel_reduce
  (const Impl::TeamThreadRangeBoundariesStruct<iType,Impl::TaskExec< Kokkos::Experimental::ROCm > >& loop_boundaries,
   const Lambda & lambda,
   const JoinType & join,
   ValueType& initialized_result)
{
   hc::tiled_index<3> idx = loop_boundaries.thread.idx();
  int team_rank = loop_boundaries.thread.team_rank(); // member num within the team
  ValueType result = initialized_result;

  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment) {
    lambda(i, result);
  }

  strided_shfl_workgroup_reduction<ValueType, JoinType>(
                          join,
                          initialized_result,
                          loop_boundaries.thread.team_size(),
                          idx.tile_dim[0]);
  initialized_result = shfl_workgroup_broadcast<ValueType>( initialized_result, idx.local[0], Impl::ROCmTraits::WavefrontSize );
}

// placeholder for future function
template< typename iType, class Lambda, typename ValueType >
KOKKOS_INLINE_FUNCTION
void parallel_reduce
  (const Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::TaskExec< Kokkos::Experimental::ROCm > >& loop_boundaries,
   const Lambda & lambda,
   ValueType& initialized_result)
{
  ValueType result = initialized_result;
  hc::tiled_index<3> idx = loop_boundaries.thread.idx();

  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment) {
    lambda(i,result);
  }

  initialized_result = result;

  //initialized_result = multi_shfl_workgroup_reduction(
  multi_shfl_workgroup_reduction(
                          [&] (ValueType& val1, const ValueType& val2) { val1 += val2; },
                          initialized_result,
                          idx.tile_dim[0]);
  initialized_result = shfl_workgroup_broadcast<ValueType>( initialized_result, 0, idx.tile_dim[0] );
}

// placeholder for future function
template< typename iType, class Lambda, typename ValueType, class JoinType >
KOKKOS_INLINE_FUNCTION
void parallel_reduce
  (const Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::TaskExec< Kokkos::Experimental::ROCm > >& loop_boundaries,
   const Lambda & lambda,
   const JoinType & join,
   ValueType& initialized_result)
{
  hc::tiled_index<3> idx = loop_boundaries.thread.idx();
  ValueType result = initialized_result;
  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment) {
    lambda(i,result);
  }
  initialized_result = result;

  multi_shfl_workgroup_reduction<ValueType, JoinType>(join, initialized_result, idx.tile_dim[0]);
  initialized_result = shfl_workgroup_broadcast<ValueType>( initialized_result, 0, idx.tile_dim[0] );
}

template< typename ValueType, typename iType, class Lambda >
KOKKOS_INLINE_FUNCTION
void parallel_scan
  (const Impl::TeamThreadRangeBoundariesStruct<iType,Impl::TaskExec< Kokkos::Experimental::ROCm > >& loop_boundaries,
   const Lambda & lambda)
{
  hc::tiled_index<3> idx = loop_boundaries.thread.idx();
  ValueType accum = 0 ;
  ValueType val, y, local_total;

  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment) {
    val = 0;
    lambda(i,val,false);

    // intra-idx.tile_dim[0] exclusive scan on 'val'
    // accum = accumulated, sum in total for this iteration

    // INCLUSIVE scan
    for( int offset = idx.tile_dim[0] ; offset < Impl::ROCmTraits::WavefrontSize ; offset <<= 1 ) {
      y = shfl_up(val, offset, Impl::ROCmTraits::WavefrontSize);
      if(idx.local[1]*idx.tile_dim[0] >= offset) { val += y; }
    }

    // pass accum to all threads
    local_total = shfl_workgroup_broadcast<ValueType>(val,
                                            idx.local[0]+Impl::ROCmTraits::WavefrontSize-idx.tile_dim[0],
                                            Impl::ROCmTraits::WavefrontSize);

    // make EXCLUSIVE scan by shifting values over one
    val = shfl_up(val, idx.tile_dim[0], Impl::ROCmTraits::WavefrontSize);
    if ( idx.local[1] == 0 ) { val = 0 ; }

    val += accum;
    lambda(i,val,true);
    accum += local_total;
  }
}

// placeholder for future function
template< typename iType, class Lambda, typename ValueType >
KOKKOS_INLINE_FUNCTION
void parallel_scan
  (const Impl::ThreadVectorRangeBoundariesStruct<iType,Impl::TaskExec< Kokkos::Experimental::ROCm > >& loop_boundaries,
   const Lambda & lambda)
{
  hc::tiled_index<3> idx = loop_boundaries.thread.idx();
  ValueType accum = 0 ;
  ValueType val, y, local_total;

  for( iType i = loop_boundaries.start; i < loop_boundaries.end; i+=loop_boundaries.increment) {
    val = 0;
    lambda(i,val,false);

    // intra-idx.tile_dim[0] exclusive scan on 'val'
    // accum = accumulated, sum in total for this iteration

    // INCLUSIVE scan
    for( int offset = 1 ; offset < idx.tile_dim[0] ; offset <<= 1 ) {
      y = shfl_up(val, offset, idx.tile_dim[0]);
      if(idx.local[0] >= offset) { val += y; }
    }

    // pass accum to all threads
    local_total = shfl_workgroup_broadcast<ValueType>(val, idx.tile_dim[0]-1, 
                                                 idx.tile_dim[0]);

    // make EXCLUSIVE scan by shifting values over one
    val = shfl_up(val, 1, idx.tile_dim[0]);
    if ( idx.local[0] == 0 ) { val = 0 ; }

    val += accum;
    lambda(i,val,true);
    accum += local_total;
  }
}


} /* namespace Kokkos */
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #if defined( KOKKOS_ENABLE_TASKDAG ) */
#endif /* #ifndef KOKKOS_IMPL_ROCM_TASK_HPP */

