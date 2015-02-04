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

#include <impl/Kokkos_Serial_TaskPolicy.hpp>

#if defined( KOKKOS_HAVE_SERIAL )
#include <stdlib.h>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <string>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

typedef TaskMember<  Kokkos::Serial , void , void > Task ;

//----------------------------------------------------------------------------

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
  free( ptr );
}

void * Task::allocate( const unsigned arg_sizeof_derived
                     , const unsigned arg_dependence_capacity )
{
  return malloc( padded_sizeof_derived( arg_sizeof_derived ) + arg_dependence_capacity * sizeof(Task*) );
}

Task::~TaskMember()
{

}

Task::TaskMember( const Task::function_verify_type   arg_verify
                , const Task::function_dealloc_type  arg_dealloc
                , const Task::function_apply_type    arg_apply
                , const unsigned                     arg_sizeof_derived
                , const unsigned                     arg_dependence_capacity
                )
  : m_dealloc( arg_dealloc )
  , m_verify(  arg_verify )
  , m_apply(   arg_apply )
  , m_dep( (Task **)( ((unsigned char *) this) + padded_sizeof_derived( arg_sizeof_derived ) ) )
  , m_wait( 0 )
  , m_next( 0 )
  , m_dep_capacity( arg_dependence_capacity )
  , m_dep_size( 0 )
  , m_ref_count( 0 )
  , m_state( TASK_STATE_CONSTRUCTING )
{
  for ( unsigned i = 0 ; i < arg_dependence_capacity ; ++i ) m_dep[i] = 0 ;
}

Task::TaskMember( const Task::function_dealloc_type  arg_dealloc
                , const Task::function_apply_type    arg_apply
                , const unsigned                     arg_sizeof_derived
                , const unsigned                     arg_dependence_capacity
                )
  : m_dealloc( arg_dealloc )
  , m_verify(  & Task::verify_type<void> )
  , m_apply(   arg_apply )
  , m_dep( (Task **)( ((unsigned char *) this) + padded_sizeof_derived( arg_sizeof_derived ) ) )
  , m_wait( 0 )
  , m_next( 0 )
  , m_dep_capacity( arg_dependence_capacity )
  , m_dep_size( 0 )
  , m_ref_count( 0 )
  , m_state( TASK_STATE_CONSTRUCTING )
{
  for ( unsigned i = 0 ; i < arg_dependence_capacity ; ++i ) m_dep[i] = 0 ;
}

//----------------------------------------------------------------------------

void Task::throw_error_add_dependence() const
{
  std::cerr << "TaskMember< Serial >::add_dependence ERROR"
            << " state(" << m_state << ")"
            << " dep_size(" << m_dep_size << ")"
            << std::endl ;
  throw std::runtime_error("TaskMember< Serial >::add_dependence ERROR");
}

void Task::throw_error_verify_type()
{
  throw std::runtime_error("TaskMember< Serial >::verify_type ERROR");
}

//----------------------------------------------------------------------------

#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )

void Task::assign( Task ** const lhs , Task * rhs , const bool no_throw )
{
  static const char msg_error_header[]      = "Kokkos::Impl::TaskManager<Kokkos::Serial>::assign ERROR" ;
  static const char msg_error_count[]       = ": negative reference count" ;
  static const char msg_error_complete[]    = ": destroy task that is not complete" ;
  static const char msg_error_dependences[] = ": destroy task that has dependences" ;
  static const char msg_error_exception[]   = ": caught internal exception" ;

  const char * msg_error = 0 ;

  try {

    if ( *lhs ) {

      const int count = --((**lhs).m_ref_count);

      if ( 0 == count ) {

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

    if ( 0 == msg_error && rhs ) { ++( rhs->m_ref_count ); }

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

namespace {

Task * s_ready = 0 ;
Task * s_denied = reinterpret_cast<Task*>( ~((unsigned long)0) );

}

void Task::schedule()
{
  // Execute ready tasks in case the task being scheduled
  // is dependent upon a waiting and ready task.

  Task::execute_ready_tasks();

  // spawning   : Constructing -> Waiting
  // respawning : Executing    -> Waiting
  // updating   : Waiting      -> Waiting

  // Must not be in a dependence linked list:  0 == t->m_next

  const bool ok_state = TASK_STATE_COMPLETE != m_state ;
  const bool ok_list  = 0 == m_next ;

  if ( ok_state && ok_list ) {

    // Will be waiting for execution upon return from this function

    m_state = Kokkos::TASK_STATE_WAITING ;

    // Insert this task into another dependence that is not complete

    int i = 0 ;
    for ( ; i < m_dep_size ; ++i ) {
      Task * const y = m_dep[i] ;
      if ( y && s_denied != ( m_next = y->m_wait ) ) {
        y->m_wait = this ; // CAS( & y->m_wait , m_next , this );
        break ;
      }
    }
    if ( i == m_dep_size ) {
      // All dependences are complete, insert into the ready list
      m_next  = s_ready ;
      s_ready = this ; // CAS( & s_ready , m_next = s_ready , this );
    }
  }
  else {
    throw std::runtime_error(std::string("Kokkos::Impl::Task spawn or respawn state error"));
  }
}

void Task::execute_ready_tasks()
{
  while ( s_ready ) {

    // Remove this task from the ready list

    // Task * task ;
    // while ( ! CAS( & s_ready , task = s_ready , s_ready->m_next ) );

    Task * const task = s_ready ;
    s_ready = task->m_next ;

    task->m_next = 0 ;

    // precondition: task->m_state = TASK_STATE_WAITING
    // precondition: task->m_dep[i]->m_state == TASK_STATE_COMPLETE  for all i
    // precondition: does not exist T such that T->m_wait = task
    // precondition: does not exist T such that T->m_next = task

    task->m_state = Kokkos::TASK_STATE_EXECUTING ;

    (*task->m_apply)( task );

    if ( task->m_state == Kokkos::TASK_STATE_EXECUTING ) {
      // task did not respawn itself
      task->m_state = Kokkos::TASK_STATE_COMPLETE ;

      // release dependences:
      for ( int i = 0 ; i < task->m_dep_size ; ++i ) {
        assign( task->m_dep + i , 0 );
      }

      // Stop other tasks from adding themselves to 'task->m_wait' ;

      Task * x ;
      // CAS( & task->m_wait , x = task->m_wait , s_denied );
      x = task->m_wait ; task->m_wait = s_denied ;

      // update tasks waiting on this task
      while ( x ) {
        Task * const next = x->m_next ;

        x->m_next = 0 ;

        x->schedule(); // could happen concurrently

        x = next ;
      }
    }
  }
}

void Task::wait( const Future< void , Kokkos::Serial > & )
{ execute_ready_tasks(); }

} // namespace Impl
} // namespace Kokkos

#endif // defined( KOKKOS_HAVE_SERIAL )
