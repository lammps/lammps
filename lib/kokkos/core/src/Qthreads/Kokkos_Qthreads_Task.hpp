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

#ifndef KOKKOS_IMPL_QTHREADS_TASK_HPP
#define KOKKOS_IMPL_QTHREADS_TASK_HPP

#include <Kokkos_Macros.hpp>
#if defined( KOKKOS_ENABLE_QTHREADS ) && defined( KOKKOS_ENABLE_TASKPOLICY )

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template<>
class TaskQueueSpecialization< Kokkos::Qthreads >
{
public:

  using execution_space = Kokkos::Qthreads ;
  using queue_type      = Kokkos::Impl::TaskQueue< execution_space > ;
  using task_base_type  = Kokkos::Impl::TaskBase< execution_space, void, void > ;

  // Must specify memory space
  using memory_space = Kokkos::HostSpace ;

  static
  void iff_single_thread_recursive_execute( queue_type * const );

  // Must provide task queue execution function
  static void execute( queue_type * const );

  // Must provide mechanism to set function pointer in
  // execution space from the host process.
  template< typename FunctorType >
  static
  void proc_set_apply( task_base_type::function_type * ptr )
    {
      using TaskType = TaskBase< execution_space,
                                 typename FunctorType::value_type,
                                 FunctorType
                               > ;
       *ptr = TaskType::apply ;
    }
};

extern template class TaskQueue< Kokkos::Qthreads > ;

//----------------------------------------------------------------------------

template<>
class TaskExec< Kokkos::Qthreads >
{
private:

  TaskExec( TaskExec && ) = delete ;
  TaskExec( TaskExec const & ) = delete ;
  TaskExec & operator = ( TaskExec && ) = delete ;
  TaskExec & operator = ( TaskExec const & ) = delete ;


  using PoolExec = Kokkos::Impl::QthreadsExec ;

  friend class Kokkos::Impl::TaskQueue< Kokkos::Qthreads > ;
  friend class Kokkos::Impl::TaskQueueSpecialization< Kokkos::Qthreads > ;

  PoolExec * const m_self_exec ;  ///< This thread's thread pool data structure
  PoolExec * const m_team_exec ;  ///< Team thread's thread pool data structure
  int64_t          m_sync_mask ;
  int64_t mutable  m_sync_value ;
  int     mutable  m_sync_step ;
  int              m_group_rank ; ///< Which "team" subset of thread pool
  int              m_team_rank ;  ///< Which thread within a team
  int              m_team_size ;

  TaskExec();
  TaskExec( PoolExec & arg_exec, int arg_team_size );

public:

#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
  void * team_shared() const
    { return m_team_exec ? m_team_exec->scratch_thread() : (void*) 0 ; }

  int team_shared_size() const
    { return m_team_exec ? m_team_exec->scratch_thread_size() : 0 ; }

  /**\brief  Whole team enters this function call
   *         before any teeam member returns from
   *         this function call.
   */
  void team_barrier() const ;
#else
  KOKKOS_INLINE_FUNCTION void team_barrier() const {}
  KOKKOS_INLINE_FUNCTION void * team_shared() const { return 0 ; }
  KOKKOS_INLINE_FUNCTION int team_shared_size() const { return 0 ; }
#endif

  KOKKOS_INLINE_FUNCTION
  int team_rank() const { return m_team_rank ; }

  KOKKOS_INLINE_FUNCTION
  int team_size() const { return m_team_size ; }
};

}} /* namespace Kokkos::Impl */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #if defined( KOKKOS_ENABLE_TASKPOLICY ) */
#endif /* #ifndef KOKKOS_IMPL_QTHREADS_TASK_HPP */

