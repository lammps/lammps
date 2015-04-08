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

#include <Kokkos_Core_fwd.hpp>

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

typedef TaskMember< Kokkos::Qthread , void , void > Task ;

namespace {

inline
unsigned padded_sizeof_derived( unsigned sizeof_derived )
{
  return sizeof_derived +
    ( sizeof_derived % sizeof(Task*) ? sizeof(Task*) - sizeof_derived % sizeof(Task*) : 0 );
}

} // namespace

void Task::deallocate( void * ptr )
{
  // Counting on 'free' thread safety so lock/unlock not required.
  // However, isolate calls here to mitigate future need to introduce lock/unlock.

  // lock

  free( ptr );

  // unlock
}

void * Task::allocate( const unsigned arg_sizeof_derived
                     , const unsigned arg_dependence_capacity )
{
  // Counting on 'malloc' thread safety so lock/unlock not required.
  // However, isolate calls here to mitigate future need to introduce lock/unlock.

  // lock

  void * const ptr = malloc( padded_sizeof_derived( arg_sizeof_derived ) + arg_dependence_capacity * sizeof(Task*) );

  // unlock

  return ptr ;
}

Task::~TaskMember()
{

}


Task::TaskMember( const function_verify_type   arg_verify
                , const function_dealloc_type  arg_dealloc
                , const function_apply_type    arg_apply
                , const unsigned               arg_sizeof_derived
                , const unsigned               arg_dependence_capacity
                )
  : m_dealloc( arg_dealloc )
  , m_verify(  arg_verify )
  , m_apply(   arg_apply )
  , m_dep( (Task **)( ((unsigned char *) this) + padded_sizeof_derived( arg_sizeof_derived ) ) )
  , m_dep_capacity( arg_dependence_capacity )
  , m_dep_size( 0 )
  , m_ref_count( 0 )
  , m_state( Kokkos::TASK_STATE_CONSTRUCTING )
  , m_qfeb(0)
{
  qthread_empty( & m_qfeb ); // Set to full when complete
  for ( unsigned i = 0 ; i < arg_dependence_capacity ; ++i ) m_dep[i] = 0 ;
}

Task::TaskMember( const function_dealloc_type  arg_dealloc
                , const function_apply_type    arg_apply
                , const unsigned               arg_sizeof_derived
                , const unsigned               arg_dependence_capacity
                )
  : m_dealloc( arg_dealloc )
  , m_verify(  & Task::verify_type<void> )
  , m_apply(   arg_apply )
  , m_dep( (Task **)( ((unsigned char *) this) + padded_sizeof_derived( arg_sizeof_derived ) ) )
  , m_dep_capacity( arg_dependence_capacity )
  , m_dep_size( 0 )
  , m_ref_count( 0 )
  , m_state( Kokkos::TASK_STATE_CONSTRUCTING )
  , m_qfeb(0)
{
  qthread_empty( & m_qfeb ); // Set to full when complete
  for ( unsigned i = 0 ; i < arg_dependence_capacity ; ++i ) m_dep[i] = 0 ;
}

//----------------------------------------------------------------------------

void Task::throw_error_add_dependence() const
{
  std::cerr << "TaskMember< Qthread >::add_dependence ERROR"
            << " state(" << m_state << ")"
            << " dep_size(" << m_dep_size << ")"
            << std::endl ;
  throw std::runtime_error("TaskMember< Qthread >::add_dependence ERROR");
}

void Task::throw_error_verify_type()
{
  throw std::runtime_error("TaskMember< Qthread >::verify_type ERROR");
}

//----------------------------------------------------------------------------

#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )
void Task::assign( Task ** const lhs , Task * rhs , const bool no_throw )
{
  static const char msg_error_header[]      = "Kokkos::Impl::TaskManager<Kokkos::Qthread>::assign ERROR" ;
  static const char msg_error_count[]       = ": negative reference count" ;
  static const char msg_error_complete[]    = ": destroy task that is not complete" ;
  static const char msg_error_dependences[] = ": destroy task that has dependences" ;
  static const char msg_error_exception[]   = ": caught internal exception" ;

  const char * msg_error = 0 ;

  try {

    if ( *lhs ) {

      const int count = Kokkos::atomic_fetch_add( & (**lhs).m_ref_count , -1 );

      if ( 1 == count ) {

        // Reference count at zero, delete it

        // Should only be deallocating a completed task
        if ( (**lhs).m_state == Kokkos::TASK_STATE_COMPLETE ) {

          // A completed task should not have dependences...
          for ( int i = 0 ; i < (**lhs).m_dep_size && 0 == msg_error ; ++i ) {
            if ( (**lhs).m_dep[i] ) msg_error = msg_error_dependences ;
          }
        }
        else {
          msg_error = msg_error_complete ;
        }

        if ( 0 == msg_error ) {
          // Get deletion function and apply it
          const Task::function_dealloc_type d = (**lhs).m_dealloc ;

          (*d)( *lhs );
        }
      }
      else if ( count <= 0 ) {
        msg_error = msg_error_count ;
      }
    }

    if ( 0 == msg_error && rhs ) { Kokkos::atomic_fetch_add( & (*rhs).m_ref_count , 1 ); }

    *lhs = rhs ;
  }
  catch( ... ) {
    if ( 0 == msg_error ) msg_error = msg_error_exception ;
  }

  if ( 0 != msg_error ) {
    if ( no_throw ) {
      std::cerr << msg_error_header << msg_error << std::endl ;
      std::cerr.flush();
    }
    else {
      std::string msg(msg_error_header);
      msg.append(msg_error);
      throw std::runtime_error( msg );
    }
  }
}
#endif


//----------------------------------------------------------------------------

aligned_t Task::qthread_func( void * arg )
{
  Task * const task = reinterpret_cast< Task * >(arg);

  task->m_state = Kokkos::TASK_STATE_EXECUTING ;

  (*task->m_apply)( task );

  if ( task->m_state == Kokkos::TASK_STATE_EXECUTING ) {
    // Task did not respawn, is complete
    task->m_state = Kokkos::TASK_STATE_COMPLETE ;

    // Release dependences before allowing dependent tasks to run.
    // Otherwise their is a thread race condition for removing dependences.
    for ( int i = 0 ; i < task->m_dep_size ; ++i ) {
      assign( & task->m_dep[i] , 0 );
    }

    // Set qthread FEB to full so that dependent tasks are allowed to execute
    qthread_fill( & task->m_qfeb );
  }

  return 0 ;
}

void Task::schedule()
{
  // Is waiting for execution

  // spawn in qthread.  must malloc the precondition array and give to qthread.
  // qthread will eventually free this allocation so memory will not be leaked.

  // concern with thread safety of malloc, does this need to be guarded?
  aligned_t ** qprecon = (aligned_t **) malloc( ( m_dep_size + 1 ) * sizeof(aligned_t *) );

  qprecon[0] = reinterpret_cast<aligned_t *>( uintptr_t(m_dep_size) );

  for ( int i = 0 ; i < m_dep_size ; ++i ) {
    qprecon[i+1] = & m_dep[i]->m_qfeb ; // Qthread precondition flag
  }

  m_state = Kokkos::TASK_STATE_WAITING ;

  qthread_spawn( & Task::qthread_func , this , 0 , NULL
               , m_dep_size , qprecon
               , NO_SHEPHERD , QTHREAD_SPAWN_SIMPLE );
}

void Task::wait( const Future< void, Kokkos::Qthread> & f )
{
  if ( f.m_task ) {
    aligned_t tmp ;
    qthread_readFF( & tmp , & f.m_task->m_qfeb );
  }
}

} // namespace Impl
} // namespace Kokkos

#endif /* #if defined( KOKKOS_HAVE_QTHREAD ) */

