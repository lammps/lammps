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

#ifndef KOKKOS_QTHREADS_HPP
#define KOKKOS_QTHREADS_HPP

#include <Kokkos_Macros.hpp>
#if defined( KOKKOS_ENABLE_QTHREADS )

#include <Kokkos_Core_fwd.hpp>

// Defines to enable experimental Qthreads functionality.
#define QTHREAD_LOCAL_PRIORITY
#define CLONED_TASKS

#include <qthread.h>

#include <cstddef>
#include <iosfwd>

#include <Kokkos_HostSpace.hpp>
#include <Kokkos_ScratchSpace.hpp>
#include <Kokkos_Parallel.hpp>
//#include <Kokkos_MemoryTraits.hpp>
//#include <Kokkos_ExecPolicy.hpp>
//#include <Kokkos_TaskScheduler.hpp> // Uncomment when Tasking working.
#include <Kokkos_Layout.hpp>
#include <impl/Kokkos_Tags.hpp>
#include <KokkosExp_MDRangePolicy.hpp>

/*--------------------------------------------------------------------------*/

namespace Kokkos {

namespace Impl {

class QthreadsExec;

} // namespace Impl

} // namespace Kokkos

/*--------------------------------------------------------------------------*/

namespace Kokkos {

/** \brief  Execution space supported by Qthreads */
class Qthreads {
public:
  //! \name Type declarations that all Kokkos devices must provide.
  //@{

  //! Tag this class as an execution space
  typedef Qthreads                 execution_space;
  typedef Kokkos::HostSpace        memory_space;
  //! This execution space preferred device_type
  typedef Kokkos::Device< execution_space, memory_space > device_type;

  typedef Kokkos::LayoutRight      array_layout;
  typedef memory_space::size_type  size_type;

  typedef ScratchMemorySpace< Qthreads > scratch_memory_space;

  //@}
  /*------------------------------------------------------------------------*/

  /** \brief  Initialization will construct one or more instances */
  static Qthreads & instance( int = 0 );

  /** \brief  Set the execution space to a "sleep" state.
   *
   * This function sets the "sleep" state in which it is not ready for work.
   * This may consume less resources than in an "ready" state,
   * but it may also take time to transition to the "ready" state.
   *
   * \return True if enters or is in the "sleep" state.
   *         False if functions are currently executing.
   */
  bool sleep();

  /** \brief  Wake from the sleep state.
   *
   *  \return True if enters or is in the "ready" state.
   *          False if functions are currently executing.
   */
  static bool wake();

  /** \brief Wait until all dispatched functions to complete.
   *
   *  The parallel_for or parallel_reduce dispatch of a functor may
   *  return asynchronously, before the functor completes.  This
   *  method does not return until all dispatched functors on this
   *  device have completed.
   */
  static void fence();

  /*------------------------------------------------------------------------*/

  static int in_parallel();

  static int is_initialized();

  /** \brief  Return maximum amount of concurrency */
  static int concurrency();

  static void initialize( int thread_count );
  static void finalize();

  /** \brief Print configuration information to the given output stream. */
  static void print_configuration( std::ostream &, const bool detail = false );

  int shepherd_size() const;
  int shepherd_worker_size() const;

  static const char* name();
};

} // namespace Kokkos

/*--------------------------------------------------------------------------*/

namespace Kokkos {

namespace Impl {

template<>
struct MemorySpaceAccess
  < Kokkos::Qthreads::memory_space
  , Kokkos::Qthreads::scratch_memory_space
  >
{
  enum { assignable = false };
  enum { accessible = true };
  enum { deepcopy   = false };
};

template<>
struct VerifyExecutionCanAccessMemorySpace
  < Kokkos::Qthreads::memory_space
  , Kokkos::Qthreads::scratch_memory_space
  >
{
  enum { value = true };
  inline static void verify( void ) {}
  inline static void verify( const void * ) {}
};

} // namespace Impl

} // namespace Kokkos

/*--------------------------------------------------------------------------*/

#include <Qthreads/Kokkos_QthreadsExec.hpp>
#include <Qthreads/Kokkos_Qthreads_Parallel.hpp>
//#include <Qthreads/Kokkos_Qthreads_Task.hpp> // Uncomment when Tasking working.
//#include <Qthreads/Kokkos_Qthreads_TaskQueue.hpp> // Uncomment when Tasking working.

#endif // #define KOKKOS_ENABLE_QTHREADS
#endif // #define KOKKOS_QTHREADS_HPP

