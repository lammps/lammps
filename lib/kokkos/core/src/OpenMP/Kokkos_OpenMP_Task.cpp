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
#if defined( KOKKOS_ENABLE_OPENMP ) && defined( KOKKOS_ENABLE_TASKDAG )

#include <Kokkos_Core.hpp>

#include <impl/Kokkos_TaskQueue_impl.hpp>
#include <impl/Kokkos_HostThreadTeam.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template class TaskQueue< Kokkos::OpenMP > ;

class HostThreadTeamDataSingleton : private HostThreadTeamData {
private:

  HostThreadTeamDataSingleton() : HostThreadTeamData()
    {
      Kokkos::OpenMP::memory_space space ;
      const size_t num_pool_reduce_bytes  =   32 ;
      const size_t num_team_reduce_bytes  =   32 ;
      const size_t num_team_shared_bytes  = 1024 ;
      const size_t num_thread_local_bytes = 1024 ;
      const size_t alloc_bytes =
        HostThreadTeamData::scratch_size( num_pool_reduce_bytes
                                        , num_team_reduce_bytes
                                        , num_team_shared_bytes
                                        , num_thread_local_bytes );

      HostThreadTeamData::scratch_assign
        ( space.allocate( alloc_bytes )
        , alloc_bytes
        , num_pool_reduce_bytes
        , num_team_reduce_bytes
        , num_team_shared_bytes
        , num_thread_local_bytes );
    }

  ~HostThreadTeamDataSingleton()
    {
      Kokkos::OpenMP::memory_space space ;
      space.deallocate( HostThreadTeamData::scratch_buffer()
                      , HostThreadTeamData::scratch_bytes() );
    }

public:

  static HostThreadTeamData & singleton()
    {
      static HostThreadTeamDataSingleton s ;
      return s ;
    }
};

//----------------------------------------------------------------------------

void TaskQueueSpecialization< Kokkos::OpenMP >::execute
  ( TaskQueue< Kokkos::OpenMP > * const queue )
{
  using execution_space = Kokkos::OpenMP ;
  using queue_type      = TaskQueue< execution_space > ;
  using task_root_type  = TaskBase< void , void , void > ;
  using Member          = Impl::HostThreadTeamMember< execution_space > ;

  static task_root_type * const end =
    (task_root_type *) task_root_type::EndTag ;


  HostThreadTeamData & team_data_single =
    HostThreadTeamDataSingleton::singleton();

  Impl::OpenMPExec * instance = t_openmp_instance;
  const int pool_size = OpenMP::thread_pool_size();

  const int team_size = 1;  // Threads per core
  instance->resize_thread_data( 0 /* global reduce buffer */
                              , 512 * team_size /* team reduce buffer */
                              , 0 /* team shared buffer */
                              , 0 /* thread local buffer */
                              );

  #pragma omp parallel num_threads(pool_size)
  {
    Impl::HostThreadTeamData & self = *(instance->get_thread_data());

    // Organizing threads into a team performs a barrier across the
    // entire pool to insure proper initialization of the team
    // rendezvous mechanism before a team rendezvous can be performed.

    if ( self.organize_team( team_size ) ) {

      Member single_exec( team_data_single );
      Member team_exec( self );

      // Loop until all queues are empty and no tasks in flight

      task_root_type * task = 0 ;

      do {
        // Each team lead attempts to acquire either a thread team task
        // or a single thread task for the team.

        if ( 0 == team_exec.team_rank() ) {

          bool leader_loop = false ;

          do {

            if ( 0 != task && end != task ) {
              // team member #0 completes the previously executed task,
              // completion may delete the task
              queue->complete( task );
            }

            // If 0 == m_ready_count then set task = 0

            task = 0 < *((volatile int *) & queue->m_ready_count) ? end : 0 ;

            // Attempt to acquire a task
            // Loop by priority and then type
            for ( int i = 0 ; i < queue_type::NumQueue && end == task ; ++i ) {
              for ( int j = 0 ; j < 2 && end == task ; ++j ) {
                task = queue_type::pop_ready_task( & queue->m_ready[i][j] );
              }
            }

            // If still tasks are still executing
            // and no task could be acquired
            // then continue this leader loop
            leader_loop = end == task ;

            if ( ( ! leader_loop ) &&
                 ( 0 != task ) &&
                 ( task_root_type::TaskSingle == task->m_task_type ) ) {

              // if a single thread task then execute now

              (*task->m_apply)( task , & single_exec );

              leader_loop = true ;
            }
          } while ( leader_loop );
        }

        // Team lead either found 0 == m_ready_count or a team task
        // Team lead broadcast acquired task:

        team_exec.team_broadcast( task , 0);

        if ( 0 != task ) { // Thread Team Task

          (*task->m_apply)( task , & team_exec );

          // The m_apply function performs a barrier
        }
      } while( 0 != task );
    }
    self.disband_team();
  }
}

void TaskQueueSpecialization< Kokkos::OpenMP >::
  iff_single_thread_recursive_execute
    ( TaskQueue< Kokkos::OpenMP > * const queue )
{
  using execution_space = Kokkos::OpenMP ;
  using queue_type      = TaskQueue< execution_space > ;
  using task_root_type  = TaskBase< void , void , void > ;
  using Member          = Impl::HostThreadTeamMember< execution_space > ;

  if ( 1 == OpenMP::thread_pool_size() ) {

    task_root_type * const end = (task_root_type *) task_root_type::EndTag ;

    HostThreadTeamData & team_data_single =
      HostThreadTeamDataSingleton::singleton();

    Member single_exec( team_data_single );

    task_root_type * task = end ;

    do {

      task = end ;

      // Loop by priority and then type
      for ( int i = 0 ; i < queue_type::NumQueue && end == task ; ++i ) {
        for ( int j = 0 ; j < 2 && end == task ; ++j ) {
          task = queue_type::pop_ready_task( & queue->m_ready[i][j] );
        }
      }

      if ( end == task ) break ;

      (*task->m_apply)( task , & single_exec );

      queue->complete( task );

    } while(1);
  }
}

}} /* namespace Kokkos::Impl */

//----------------------------------------------------------------------------
#else
void KOKKOS_CORE_SRC_OPENMP_KOKKOS_OPENMP_TASK_PREVENT_LINK_ERROR() {}
#endif /* #if defined( KOKKOS_ENABLE_OPENMP ) && defined( KOKKOS_ENABLE_TASKDAG ) */

