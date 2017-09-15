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

#if defined( KOKKOS_DEBUG )

__device__
void verify_warp_convergence( const char * const where )
{
  const unsigned b = __ballot(1);

  if ( b != ~0u ) {

printf(" verify_warp_convergence( %s ) (%d,%d,%d) (%d,%d,%d) failed %x\n"
      , where
      , blockIdx.x
      , blockIdx.y
      , blockIdx.z
      , threadIdx.x
      , threadIdx.y
      , threadIdx.z
      , b );

  }
}

#endif // #if defined( KOKKOS_DEBUG )

//----------------------------------------------------------------------------

__device__
void TaskQueueSpecialization< Kokkos::Cuda >::driver
  ( TaskQueueSpecialization< Kokkos::Cuda >::queue_type * const queue 
  , int32_t shmem_per_warp )
{
  using Member = TaskExec< Kokkos::Cuda > ;
  using Queue  = TaskQueue< Kokkos::Cuda > ;
  using task_root_type = TaskBase< void , void , void > ;

  extern __shared__ int32_t shmem_all[];

  task_root_type * const end = (task_root_type *) task_root_type::EndTag ;

  int32_t * const warp_shmem =
    shmem_all + ( threadIdx.z * shmem_per_warp ) / sizeof(int32_t);

  task_root_type * const task_shmem = (task_root_type *) warp_shmem ;

  const int warp_lane = threadIdx.x + threadIdx.y * blockDim.x ;

  Member single_exec( warp_shmem , 1 );
  Member team_exec( warp_shmem , blockDim.y );

  task_root_type * task_ptr ;

  // Loop until all queues are empty and no tasks in flight

  do {

    // Each team lead attempts to acquire either a thread team task
    // or collection of single thread tasks for the team.

    if ( 0 == warp_lane ) {

      task_ptr = 0 < *((volatile int *) & queue->m_ready_count) ? end : 0 ;

      // Loop by priority and then type
      for ( int i = 0 ; i < Queue::NumQueue && end == task_ptr ; ++i ) {
        for ( int j = 0 ; j < 2 && end == task_ptr ; ++j ) {
          task_ptr = Queue::pop_ready_task( & queue->m_ready[i][j] );
        }
      }

#if 0
printf("TaskQueue<Cuda>::driver(%d,%d) task(%lx)\n",threadIdx.z,blockIdx.x
      , uintptr_t(task_ptr));
#endif

    }

    // shuffle broadcast

    ((int*) & task_ptr )[0] = __shfl( ((int*) & task_ptr )[0] , 0 );
    ((int*) & task_ptr )[1] = __shfl( ((int*) & task_ptr )[1] , 0 );

#if defined( KOKKOS_DEBUG )
    verify_warp_convergence("task_ptr");
#endif

    if ( 0 == task_ptr ) break ; // 0 == queue->m_ready_count

    if ( end != task_ptr ) {

      // Whole warp copy task's closure to/from shared memory.
      // Use all threads of warp for coalesced read/write.

      int32_t const b = sizeof(task_root_type) / sizeof(int32_t);
      int32_t const e = *((int32_t volatile *)( & task_ptr->m_alloc_size )) / sizeof(int32_t);

      int32_t volatile * const task_mem = (int32_t volatile *) task_ptr ;

      // copy global to shared memory:

      for ( int32_t i = warp_lane ; i < e ; i += CudaTraits::WarpSize ) {
        warp_shmem[i] = task_mem[i] ;
      }

      Kokkos::memory_fence();

      // Copy done - use memory fence so that memory writes are visible.
      // For reliable warp convergence on Pascal and Volta an explicit
      // warp level synchronization will also be required.

      if ( task_root_type::TaskTeam == task_shmem->m_task_type ) {
        // Thread Team Task
        (*task_shmem->m_apply)( task_shmem , & team_exec );
      }
      else if ( 0 == threadIdx.y ) {
        // Single Thread Task
        (*task_shmem->m_apply)( task_shmem , & single_exec );
      }

      // copy shared to global memory:

      for ( int32_t i = b + warp_lane ; i < e ; i += CudaTraits::WarpSize ) {
        task_mem[i] = warp_shmem[i] ;
      }

      Kokkos::memory_fence();

#if defined( KOKKOS_DEBUG )
    verify_warp_convergence("apply");
#endif

      // If respawn requested copy respawn data back to main memory

      if ( 0 == warp_lane ) {

        if ( ((task_root_type *) task_root_type::LockTag) != task_shmem->m_next ) {
          ( (volatile task_root_type *) task_ptr )->m_next = task_shmem->m_next ;
          ( (volatile task_root_type *) task_ptr )->m_priority = task_shmem->m_priority ;
        }

        queue->complete( task_ptr );
      }
    }
  } while(1);
}

namespace {

__global__
void cuda_task_queue_execute( TaskQueue< Kokkos::Cuda > * queue 
                            , int32_t shmem_size )
{ TaskQueueSpecialization< Kokkos::Cuda >::driver( queue , shmem_size ); }

}

void TaskQueueSpecialization< Kokkos::Cuda >::execute
  ( TaskQueue< Kokkos::Cuda > * const queue )
{
  const int shared_per_warp = 2048 ;
  const int warps_per_block = 4 ;
  const dim3 grid( Kokkos::Impl::cuda_internal_multiprocessor_count() , 1 , 1 );
  const dim3 block( 1 , Kokkos::Impl::CudaTraits::WarpSize , warps_per_block );
  const int shared_total = shared_per_warp * warps_per_block ;
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

  cuda_task_queue_execute<<< grid , block , shared_total , stream >>>( queue , shared_per_warp );

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

