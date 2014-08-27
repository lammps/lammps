/*
//@HEADER
// ************************************************************************
// 
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
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

// Experimental unified task-data parallel manycore LDRD

#include <Kokkos_Macros.hpp>

#if defined( KOKKOS_HAVE_QTHREAD )

#include <stdio.h>

#include <stdlib.h>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <string>

#include <Kokkos_Atomic.hpp>
#include <Qthread/Kokkos_Qthread_TaskPolicy.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {
namespace {

TaskManager< Kokkos::Qthread > s_task_manager ;

}

typedef TaskMember<  Kokkos::Qthread > Task ;
typedef TaskManager< Kokkos::Qthread > Mgr ;

Task::TaskMember( const function_type    arg_destroy
                , const function_type    arg_apply
                , const std::type_info & arg_type
                )
  : m_typeid(  arg_type )
  , m_destroy( arg_destroy )
  , m_apply(   arg_apply )
  , m_state( STATE_CONSTRUCTING )
  , m_ref_count(0)
  , m_qfeb(0)
{
  qthread_empty( & m_qfeb ); // Set to full when complete
  for ( int i = 0 ; i < MAX_DEPENDENCE ; ++i ) m_dep[i] = 0 ;
}

Mgr::TaskManager()
{}

void * Mgr::memory_allocate( size_t nbytes )
{
  // Counting on 'malloc' thread safety so lock/unlock not required.
  // However, isolate calls here to mitigate future need to introduce lock/unlock.

  // lock

  void * ptr = malloc( nbytes );

  // unlock

  return ptr ;
}

void Mgr::memory_deallocate( void * ptr )
{
  // Counting on 'free' thread safety so lock/unlock not required.
  // However, isolate calls here to mitigate future need to introduce lock/unlock.

  // lock

  free( ptr );

  // unlock
}

void Mgr::assign( Task ** const lhs , Task * const rhs )
{
  if ( *lhs ) {

    // Must de-assign

    const int count = Kokkos::atomic_fetch_add( & (**lhs).m_ref_count , -1 );

    if ( 1 == count ) {

      // Should only be deallocating a completed task
      // TODO: Support deletion of canceled tasks.

      if ( (**lhs).m_state != Task::STATE_COMPLETE ) {
        throw std::runtime_error(
          std::string("Kokkos::Impl::TaskManager<Kokkos::Qthread>::decrement ERROR: not STATE_COMPLETE") );
      }

      // Get destructor function and apply it
      (**lhs).m_destroy( *lhs );

      memory_deallocate( *lhs );
    }
    else if ( count <= 0 ) {
      throw std::runtime_error(std::string("Kokkos::Impl::TaskManager<Kokkos::Qthread>::assign ERROR: reference counting") );
    }
  }

  if ( rhs ) {
    Kokkos::atomic_fetch_add( & (*rhs).m_ref_count , 1 );
  }

  *lhs = rhs ;
}

void Mgr::verify_set_dependence( Task * t , int n )
{
  // Must be either constructing for original spawn or executing for a respawn.

  if ( Task::STATE_CONSTRUCTING != t->m_state &&
       Task::STATE_EXECUTING    != t->m_state ) {
    throw std::runtime_error(std::string("Kokkos::Impl::Task spawn or respawn state error"));
  }

  if ( MAX_DEPENDENCE <= n ) {
    throw std::runtime_error(std::string("Kokkos::Impl::Task spawn or respawn dependence count error"));
  }
}

void Mgr::schedule( Task * t )
{
  // Is waiting for execution

  // spawn in qthread.  must malloc the precondition array and give to qthread.
  // qthread will eventually free this allocation so memory will not be leaked.

  // concern with thread safety of malloc, does this need to be guarded?
  aligned_t ** qprecon = (aligned_t **) memory_allocate( ( MAX_DEPENDENCE + 1 ) * sizeof(aligned_t *) );

  uintptr_t npre = 0 ;
  for ( ; npre < MAX_DEPENDENCE && t->m_dep[npre] ; ++npre ) {
    qprecon[npre+1] = & t->m_dep[npre]->m_qfeb ; // Qthread precondition flag
  }
  qprecon[0] = reinterpret_cast<aligned_t *>( npre );

  t->m_state = Task::STATE_WAITING ;

  qthread_spawn( & Mgr::qthread_func , t , 0 , NULL
               , npre , qprecon
               , NO_SHEPHERD , QTHREAD_SPAWN_SIMPLE );
}

aligned_t Mgr::qthread_func( void * arg )
{
  Task * const task = reinterpret_cast< Task * >(arg);

  task->m_state = Task::STATE_EXECUTING ;

  (*task->m_apply)( task );

  if ( task->m_state == Task::STATE_EXECUTING ) {
    // Task did not respawn, is complete
    task->m_state = Task::STATE_COMPLETE ;

    // Release dependences before allowing dependent tasks to run.
    // Otherwise their is a thread race condition for removing dependences.
    for ( int i = 0 ; i < MAX_DEPENDENCE ; ++i ) {
      assign( & task->m_dep[i] , 0 );
    }

    // Set qthread FEB to full so that dependent tasks are allowed to execute
    qthread_fill( & task->m_qfeb );
  }

  return 0 ;
}


void Mgr::wait( Task * t )
{
  aligned_t tmp ;
  qthread_readFF( & tmp , & t->m_qfeb );
}

} // namespace Impl
} // namespace Kokkos

namespace Kokkos {

TaskPolicy< Kokkos::Qthread >::TaskPolicy()
  : m_task_manager( Impl::s_task_manager )
{}

} // namespace Kokkos

#endif /* #if defined( KOKKOS_HAVE_QTHREAD ) */

