/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
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

#ifndef KOKKOS_THREADS_HPP
#define KOKKOS_THREADS_HPP

#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ENABLE_THREADS)

#include <Kokkos_Core_fwd.hpp>

#include <cstddef>
#include <iosfwd>
#include <Kokkos_HostSpace.hpp>
#include <Kokkos_ScratchSpace.hpp>
#include <Kokkos_Layout.hpp>
#include <Kokkos_MemoryTraits.hpp>
#include <impl/Kokkos_Profiling_Interface.hpp>
#include <impl/Kokkos_Tags.hpp>

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {
class ThreadsExec;
}  // namespace Impl
}  // namespace Kokkos

/*--------------------------------------------------------------------------*/

namespace Kokkos {

/** \brief  Execution space for a pool of Pthreads or C11 threads on a CPU. */
class Threads {
 public:
  //! \name Type declarations that all Kokkos devices must provide.
  //@{
  //! Tag this class as a kokkos execution space
  typedef Threads execution_space;
  typedef Kokkos::HostSpace memory_space;

  //! This execution space preferred device_type
  typedef Kokkos::Device<execution_space, memory_space> device_type;

  typedef Kokkos::LayoutRight array_layout;
  typedef memory_space::size_type size_type;

  typedef ScratchMemorySpace<Threads> scratch_memory_space;

  //@}
  /*------------------------------------------------------------------------*/
  //! \name Static functions that all Kokkos devices must implement.
  //@{

  /// \brief True if and only if this method is being called in a
  ///   thread-parallel function.
  static int in_parallel();

  /// \brief Print configuration information to the given output stream.
  static void print_configuration(std::ostream&, const bool detail = false);

  /// \brief Wait until all dispatched functors complete.
  ///
  /// The parallel_for or parallel_reduce dispatch of a functor may
  /// return asynchronously, before the functor completes.  This
  /// method does not return until all dispatched functors on this
  /// device have completed.
  static void impl_static_fence();

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
  static void fence();
#else
  void fence() const;
#endif

  /** \brief  Return the maximum amount of concurrency.  */
  static int concurrency();

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
  static bool sleep();

  static bool wake();

  static void finalize();

  static void initialize(unsigned threads_count             = 0,
                         unsigned use_numa_count            = 0,
                         unsigned use_cores_per_numa        = 0,
                         bool allow_asynchronous_threadpool = false);

  static int is_initialized();

  static Threads& instance(int = 0);

  //----------------------------------------

  static int thread_pool_size(int depth = 0);
#if defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST)
  static int thread_pool_rank();
#else
  KOKKOS_INLINE_FUNCTION static int thread_pool_rank() { return 0; }
#endif

  inline static unsigned max_hardware_threads() { return thread_pool_size(0); }
  KOKKOS_INLINE_FUNCTION static unsigned hardware_thread_id() {
    return thread_pool_rank();
  }
#else
  /// \brief Free any resources being consumed by the device.
  ///
  /// For the Threads device, this terminates spawned worker threads.
  static void impl_finalize();

  //@}
  /*------------------------------------------------------------------------*/
  /*------------------------------------------------------------------------*/
  //! \name Space-specific functions
  //@{

  /** \brief Initialize the device in the "ready to work" state.
   *
   *  The device is initialized in a "ready to work" or "awake" state.
   *  This state reduces latency and thus improves performance when
   *  dispatching work.  However, the "awake" state consumes resources
   *  even when no work is being done.  You may call sleep() to put
   *  the device in a "sleeping" state that does not consume as many
   *  resources, but it will take time (latency) to awaken the device
   *  again (via the wake()) method so that it is ready for work.
   *
   *  Teams of threads are distributed as evenly as possible across
   *  the requested number of numa regions and cores per numa region.
   *  A team will not be split across a numa region.
   *
   *  If the 'use_' arguments are not supplied the hwloc is queried
   *  to use all available cores.
   */
  static void impl_initialize(unsigned threads_count             = 0,
                              unsigned use_numa_count            = 0,
                              unsigned use_cores_per_numa        = 0,
                              bool allow_asynchronous_threadpool = false);

  static int impl_is_initialized();

  static Threads& impl_instance(int = 0);

  //----------------------------------------

  static int impl_thread_pool_size(int depth = 0);
#if defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST)
  static int impl_thread_pool_rank();
#else
  KOKKOS_INLINE_FUNCTION static int impl_thread_pool_rank() { return 0; }
#endif

  inline static unsigned impl_max_hardware_threads() {
    return impl_thread_pool_size(0);
  }
  KOKKOS_INLINE_FUNCTION static unsigned impl_hardware_thread_id() {
    return impl_thread_pool_rank();
  }
#endif

  uint32_t impl_instance_id() const noexcept { return 0; }

  static const char* name();
  //@}
  //----------------------------------------
};

namespace Profiling {
namespace Experimental {
template <>
struct DeviceTypeTraits<Threads> {
  static constexpr DeviceType id = DeviceType::Threads;
};
}  // namespace Experimental
}  // namespace Profiling
}  // namespace Kokkos

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

template <>
struct MemorySpaceAccess<Kokkos::Threads::memory_space,
                         Kokkos::Threads::scratch_memory_space> {
  enum { assignable = false };
  enum { accessible = true };
  enum { deepcopy = false };
};

template <>
struct VerifyExecutionCanAccessMemorySpace<
    Kokkos::Threads::memory_space, Kokkos::Threads::scratch_memory_space> {
  enum { value = true };
  inline static void verify(void) {}
  inline static void verify(const void*) {}
};

}  // namespace Impl
}  // namespace Kokkos

/*--------------------------------------------------------------------------*/

#include <Kokkos_ExecPolicy.hpp>
#include <Kokkos_Parallel.hpp>
#include <Threads/Kokkos_ThreadsExec.hpp>
#include <Threads/Kokkos_ThreadsTeam.hpp>
#include <Threads/Kokkos_Threads_Parallel.hpp>

#include <KokkosExp_MDRangePolicy.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #if defined( KOKKOS_ENABLE_THREADS ) */
#endif /* #define KOKKOS_THREADS_HPP */
