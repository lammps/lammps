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

#include <Kokkos_Core_fwd.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if defined( KOKKOS_HAVE_PTHREAD )

/* Standard 'C' Linux libraries */

#include <pthread.h>
#include <sched.h>
#include <errno.h>

/* Standard C++ libaries */

#include <cstdlib>
#include <string>
#include <iostream>
#include <stdexcept>

#include <Kokkos_Threads.hpp>

//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {
namespace {

pthread_mutex_t host_internal_pthread_mutex = PTHREAD_MUTEX_INITIALIZER ;

// Pthreads compatible driver.
// Recovery from an exception would require constant intra-thread health
// verification; which would negatively impact runtime.  As such simply
// abort the process.

void * internal_pthread_driver( void * )
{
  try {
    ThreadsExec::driver();
  }
  catch( const std::exception & x ) {
    std::cerr << "Exception thrown from worker thread: " << x.what() << std::endl ;
    std::cerr.flush();
    std::abort();
  }
  catch( ... ) {
    std::cerr << "Exception thrown from worker thread" << std::endl ;
    std::cerr.flush();
    std::abort();
  }
  return NULL ;
}

} // namespace

//----------------------------------------------------------------------------
// Spawn a thread

bool ThreadsExec::spawn()
{
  bool result = false ;

  pthread_attr_t attr ;

  if ( 0 == pthread_attr_init( & attr ) ||
       0 == pthread_attr_setscope(       & attr, PTHREAD_SCOPE_SYSTEM ) ||
       0 == pthread_attr_setdetachstate( & attr, PTHREAD_CREATE_DETACHED ) ) {

    pthread_t pt ;

    result = 0 == pthread_create( & pt, & attr, internal_pthread_driver, 0 );
  }

  pthread_attr_destroy( & attr );

  return result ;
}

//----------------------------------------------------------------------------

bool ThreadsExec::is_process()
{
  static const pthread_t master_pid = pthread_self();

  return pthread_equal( master_pid , pthread_self() );
}

void ThreadsExec::global_lock()
{
  pthread_mutex_lock( & host_internal_pthread_mutex );
}

void ThreadsExec::global_unlock()
{
  pthread_mutex_unlock( & host_internal_pthread_mutex );
}

//----------------------------------------------------------------------------

void ThreadsExec::wait_yield( volatile int & flag , const int value )
{
  while ( value == flag ) { sched_yield(); }
}

} // namespace Impl
} // namespace Kokkos

/* end #if defined( KOKKOS_HAVE_PTHREAD ) */
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#elif defined( KOKKOS_HAVE_WINTHREAD )

/* Windows libraries */
#include <windows.h>
#include <process.h>

/* Standard C++ libaries */

#include <cstdlib>
#include <string>
#include <iostream>
#include <stdexcept>

#include <Kokkos_Threads.hpp>

//----------------------------------------------------------------------------
// Driver for each created pthread

namespace Kokkos {
namespace Impl {
namespace {

unsigned WINAPI internal_winthread_driver( void * arg )
{
  ThreadsExec::driver();

  return 0 ;
}

class ThreadLockWindows {
private:
  CRITICAL_SECTION  m_handle ;

  ~ThreadLockWindows()
  { DeleteCriticalSection( & m_handle ); }

  ThreadLockWindows();
  { InitializeCriticalSection( & m_handle ); }

  ThreadLockWindows( const ThreadLockWindows & );
  ThreadLockWindows & operator = ( const ThreadLockWindows & );

public:

  static ThreadLockWindows & singleton();

  void lock()
  { EnterCriticalSection( & m_handle ); }

  void unlock()
  { LeaveCriticalSection( & m_handle ); }
};

ThreadLockWindows & ThreadLockWindows::singleton()
{ static ThreadLockWindows self ; return self ; }

} // namespace <>
} // namespace Kokkos
} // namespace Impl

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

// Spawn this thread

bool ThreadsExec::spawn()
{
  unsigned Win32ThreadID = 0 ;

  HANDLE handle =
    _beginthreadex(0,0,internal_winthread_driver,0,0, & Win32ThreadID );

  return ! handle ;
}

bool ThreadsExec::is_process() { return true ; }

void ThreadsExec::global_lock()
{ ThreadLockWindows::singleton().lock(); }

void ThreadsExec::global_unlock()
{ ThreadLockWindows::singleton().unlock(); }

void ThreadsExec::wait_yield( volatile int & flag , const int value ) {}
{
  while ( value == flag ) { Sleep(0); }
}

} // namespace Impl
} // namespace Kokkos

#endif /* end #elif defined( KOKKOS_HAVE_WINTHREAD ) */
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------



