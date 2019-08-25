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

#if defined( KOKKOS_ENABLE_ROCM ) && defined( KOKKOS_ENABLE_TASKDAG )

#include <impl/Kokkos_TaskQueue_impl.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template class TaskQueue< Kokkos::Experimental::ROCm > ;


//----------------------------------------------------------------------------
KOKKOS_INLINE_FUNCTION
void TaskQueueSpecialization< Kokkos::Experimental::ROCm >::driver
  ( TaskQueueSpecialization< Kokkos::Experimental::ROCm >::queue_type * const queue,
    hc::tiled_index<3> threadIdx )
{
  using Member = TaskExec< Kokkos::Experimental::ROCm > ;
  using Queue  = TaskQueue< Kokkos::Experimental::ROCm > ;
  using task_root_type = TaskBase< void , void , void > ;

  task_root_type * const end = (task_root_type *) task_root_type::EndTag ;

  Member single_exec( 1, threadIdx );
  Member team_exec( threadIdx.tile_dim[0], threadIdx );

  const int wavefront_lane = threadIdx.local[0] + threadIdx.local[1]* threadIdx.tile_dim[0] ;

  union {
    task_root_type * ptr ;
    int              raw[2] ;
  } task ;

  // Loop until all queues are empty and no tasks in flight

  do {

    // Each team lead attempts to acquire either a thread team task
    // or collection of single thread tasks for the team.

    if ( 0 == wavefront_lane ) {

      task.ptr = 0 < *((volatile int *) & queue->m_ready_count) ? end : 0 ;

      // Loop by priority and then type
      for ( int i = 0 ; i < Queue::NumQueue && end == task.ptr ; ++i ) {
        for ( int j = 0 ; j < 2 && end == task.ptr ; ++j ) {
          task.ptr = Queue::pop_ready_task( & queue->m_ready[i][j] );
        }
      }

#if 0
printf("TaskQueue<ROCm>::driver(%d,%d) task(%lx)\n",threadIdx.z,blockIdx.x
      , uintptr_t(task.ptr));
#endif

    }

    // shuffle broadcast

    task.raw[0] = hc::__shfl( task.raw[0] , 0 );
    task.raw[1] = hc::__shfl( task.raw[1] , 0 );

    if ( 0 == task.ptr ) break ; // 0 == queue->m_ready_count

    if ( end != task.ptr ) {
      if ( task_root_type::TaskTeam == task.ptr->m_task_type ) {
        // Thread Team Task
        (*task.ptr->m_apply)( task.ptr , & team_exec );
      }
      else if ( 0 == threadIdx.local[1] ) {
        // Single Thread Task
        (*task.ptr->m_apply)( task.ptr , & single_exec );
      }

      if ( 0 == wavefront_lane ) {
        queue->complete( task.ptr );
      }
    }
  } while(1);
}
#if 0
namespace {
KOKKOS_INLINE_FUNCTION
void rocm_task_queue_execute( TaskQueue< Kokkos::Experimental::ROCm > * queue, 
                              hc::tiled_index<3> threadIdx )
{ TaskQueueSpecialization< Kokkos::Experimental::ROCm >::driver( queue, threadIdx ); }

}
#endif
void TaskQueueSpecialization< Kokkos::Experimental::ROCm >::execute
  ( TaskQueue< Kokkos::Experimental::ROCm > * const queue )
{
  const int workgroups_per_wavefront = 4 ;
  const int wavefront_size = Kokkos::Impl::ROCmTraits::WavefrontSize ;
  const int cu_count = Kokkos::Impl::rocm_internal_cu_count();
//  const dim3 grid( Kokkos::Impl::rocm_internal_cu_count() , 1 , 1 );
//  const dim3 block( 1 , Kokkos::Impl::ROCmTraits::WorkGroupSize , workgroups_per_wavefront );



  // Query the stack size, in bytes:
  // If not large enough then set the stack size, in bytes:

// adapted from the cuda code.  TODO: Not at all sure that this is the proper 
// to map the cuda grid/blocks/3D tiling to HCC
#if 0
  hc::extent< 3 > flat_extent(  cu_count,
                                wavefront_size, workgroups_per_wavefront );
  hc::tiled_extent< 3 > team_extent = flat_extent.tile(1,
                                wavefront_size,workgroups_per_wavefront);

  hc::parallel_for_each( team_extent , [&](hc::tiled_index<3> idx) [[hc]]
  {
    TaskQueueSpecialization< Kokkos::Experimental::ROCm >::driver( queue,idx ); 
  }).wait();
#endif
}


}} /* namespace Kokkos::Impl */

//----------------------------------------------------------------------------

#endif /* #if defined( KOKKOS_ENABLE_ROCM ) && defined( KOKKOS_ENABLE_TASKDAG ) */


