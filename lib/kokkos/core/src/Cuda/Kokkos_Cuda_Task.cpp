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

#include <Kokkos_Macros.hpp>
#if defined( KOKKOS_ENABLE_CUDA ) && defined( KOKKOS_ENABLE_TASKDAG )

#include <Kokkos_Core.hpp>

#include <impl/Kokkos_TaskQueue_impl.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template class TaskQueue< Kokkos::Cuda > ;

//----------------------------------------------------------------------------

__device__
void TaskQueueSpecialization< Kokkos::Cuda >::driver
  ( TaskQueueSpecialization< Kokkos::Cuda >::queue_type * const queue )
{
  using Member = TaskExec< Kokkos::Cuda > ;
  using Queue  = TaskQueue< Kokkos::Cuda > ;
  using task_root_type = TaskBase< Kokkos::Cuda , void , void > ;

  task_root_type * const end = (task_root_type *) task_root_type::EndTag ;

  Member single_exec( 1 );
  Member team_exec( blockDim.y );

  const int warp_lane = threadIdx.x + threadIdx.y * blockDim.x ;

  union {
    task_root_type * ptr ;
    int              raw[2] ;
  } task ;

  // Loop until all queues are empty and no tasks in flight

  do {

    // Each team lead attempts to acquire either a thread team task
    // or collection of single thread tasks for the team.

    if ( 0 == warp_lane ) {

      task.ptr = 0 < *((volatile int *) & queue->m_ready_count) ? end : 0 ;

      // Loop by priority and then type
      for ( int i = 0 ; i < Queue::NumQueue && end == task.ptr ; ++i ) {
        for ( int j = 0 ; j < 2 && end == task.ptr ; ++j ) {
          task.ptr = Queue::pop_ready_task( & queue->m_ready[i][j] );
        }
      }

#if 0
printf("TaskQueue<Cuda>::driver(%d,%d) task(%lx)\n",threadIdx.z,blockIdx.x
      , uintptr_t(task.ptr));
#endif

    }

    // shuffle broadcast

    task.raw[0] = __shfl( task.raw[0] , 0 );
    task.raw[1] = __shfl( task.raw[1] , 0 );

    if ( 0 == task.ptr ) break ; // 0 == queue->m_ready_count

    if ( end != task.ptr ) {
      if ( task_root_type::TaskTeam == task.ptr->m_task_type ) {
        // Thread Team Task
        (*task.ptr->m_apply)( task.ptr , & team_exec );
      }
      else if ( 0 == threadIdx.y ) {
        // Single Thread Task
        (*task.ptr->m_apply)( task.ptr , & single_exec );
      }

      if ( 0 == warp_lane ) {
        queue->complete( task.ptr );
      }
    }
  } while(1);
}

namespace {

__global__
void cuda_task_queue_execute( TaskQueue< Kokkos::Cuda > * queue )
{ TaskQueueSpecialization< Kokkos::Cuda >::driver( queue ); }

}

void TaskQueueSpecialization< Kokkos::Cuda >::execute
  ( TaskQueue< Kokkos::Cuda > * const queue )
{
  const int warps_per_block = 4 ;
  const dim3 grid( Kokkos::Impl::cuda_internal_multiprocessor_count() , 1 , 1 );
  const dim3 block( 1 , Kokkos::Impl::CudaTraits::WarpSize , warps_per_block );
  const int shared = 0 ;
  const cudaStream_t stream = 0 ;

  CUDA_SAFE_CALL( cudaDeviceSynchronize() );

#if 0
printf("cuda_task_queue_execute before\n");
#endif

  // Query the stack size, in bytes:
  //
  // size_t stack_size = 0 ;
  // CUDA_SAFE_CALL( cudaDeviceGetLimit( & stack_size , cudaLimitStackSize ) );
  //
  // If not large enough then set the stack size, in bytes:
  //
  // CUDA_SAFE_CALL( cudaDeviceSetLimit( cudaLimitStackSize , stack_size ) );

  cuda_task_queue_execute<<< grid , block , shared , stream >>>( queue );

  CUDA_SAFE_CALL( cudaGetLastError() );

  CUDA_SAFE_CALL( cudaDeviceSynchronize() );

#if 0
printf("cuda_task_queue_execute after\n");
#endif

}

}} /* namespace Kokkos::Impl */

//----------------------------------------------------------------------------
#else
void KOKKOS_CORE_SRC_CUDA_KOKKOS_CUDA_TASK_PREVENT_LINK_ERROR() {}
#endif /* #if defined( KOKKOS_ENABLE_CUDA ) && defined( KOKKOS_ENABLE_TASKDAG ) */

