//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#include <Kokkos_Macros.hpp>
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#endif
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
#include <impl/Kokkos_InitializationSettings.hpp>

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {
class ThreadsExec;
enum class fence_is_static { yes, no };
}  // namespace Impl
}  // namespace Kokkos

/*--------------------------------------------------------------------------*/

namespace Kokkos {

/** \brief  Execution space for a pool of C++11 threads on a CPU. */
class Threads {
 public:
  //! \name Type declarations that all Kokkos devices must provide.
  //@{
  //! Tag this class as a kokkos execution space
  using execution_space = Threads;
  using memory_space    = Kokkos::HostSpace;

  //! This execution space preferred device_type
  using device_type = Kokkos::Device<execution_space, memory_space>;

  using array_layout = Kokkos::LayoutRight;
  using size_type    = memory_space::size_type;

  using scratch_memory_space = ScratchMemorySpace<Threads>;

  //@}
  /*------------------------------------------------------------------------*/
  //! \name Static functions that all Kokkos devices must implement.
  //@{

  /// \brief True if and only if this method is being called in a
  ///   thread-parallel function.
  static int in_parallel();

  /// \brief Print configuration information to the given output stream.
  void print_configuration(std::ostream& os, bool verbose = false) const;

  /// \brief Wait until all dispatched functors complete.
  ///
  /// The parallel_for or parallel_reduce dispatch of a functor may
  /// return asynchronously, before the functor completes.  This
  /// method does not return until all dispatched functors on this
  /// device have completed.
  static void impl_static_fence(const std::string& name);

  void fence(const std::string& name =
                 "Kokkos::Threads::fence: Unnamed Instance Fence") const;

  /** \brief  Return the maximum amount of concurrency.  */
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  static int concurrency();
#else
  int concurrency() const;
#endif

  /// \brief Free any resources being consumed by the device.
  ///
  /// For the Threads device, this terminates spawned worker threads.
  static void impl_finalize();

  //@}
  /*------------------------------------------------------------------------*/
  /*------------------------------------------------------------------------*/
  //! \name Space-specific functions
  //@{

  static void impl_initialize(InitializationSettings const&);

  static int impl_is_initialized();

  static Threads& impl_instance(int = 0);

  //----------------------------------------

  static int impl_thread_pool_size(int depth = 0);

  static int impl_thread_pool_rank_host();

  static KOKKOS_FUNCTION int impl_thread_pool_rank() {
    KOKKOS_IF_ON_HOST((return impl_thread_pool_rank_host();))

    KOKKOS_IF_ON_DEVICE((return 0;))
  }

  inline static unsigned impl_max_hardware_threads() {
    return impl_thread_pool_size(0);
  }
  KOKKOS_INLINE_FUNCTION static unsigned impl_hardware_thread_id() {
    return impl_thread_pool_rank();
  }

  uint32_t impl_instance_id() const noexcept { return 1; }

  static const char* name();
  //@}
  //----------------------------------------
 private:
  friend bool operator==(Threads const&, Threads const&) { return true; }
  friend bool operator!=(Threads const&, Threads const&) { return false; }
};

namespace Tools {
namespace Experimental {
template <>
struct DeviceTypeTraits<Threads> {
  static constexpr DeviceType id = DeviceType::Threads;
  static int device_id(const Threads&) { return 0; }
};
}  // namespace Experimental
}  // namespace Tools
}  // namespace Kokkos

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

template <>
struct MemorySpaceAccess<Kokkos::Threads::memory_space,
                         Kokkos::Threads::scratch_memory_space> {
  enum : bool { assignable = false };
  enum : bool { accessible = true };
  enum : bool { deepcopy = false };
};

}  // namespace Impl
}  // namespace Kokkos

/*--------------------------------------------------------------------------*/

#include <Kokkos_ExecPolicy.hpp>
#include <Kokkos_Parallel.hpp>
#include <Threads/Kokkos_ThreadsExec.hpp>
#include <Threads/Kokkos_ThreadsTeam.hpp>
#include <Threads/Kokkos_Threads_Parallel_Range.hpp>
#include <Threads/Kokkos_Threads_Parallel_MDRange.hpp>
#include <Threads/Kokkos_Threads_Parallel_Team.hpp>
#include <Threads/Kokkos_Threads_UniqueToken.hpp>

#include <KokkosExp_MDRangePolicy.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #if defined( KOKKOS_ENABLE_THREADS ) */
#endif /* #define KOKKOS_THREADS_HPP */
