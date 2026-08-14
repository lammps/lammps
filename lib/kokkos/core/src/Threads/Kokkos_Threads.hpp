// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

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
#include <impl/Kokkos_CheckUsage.hpp>
#include <impl/Kokkos_InitializationSettings.hpp>
#include <impl/Kokkos_Profiling_Interface.hpp>

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
  int concurrency() const;

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

  KOKKOS_DEFAULTED_FUNCTION Threads(const Threads&) = default;
  KOKKOS_FUNCTION Threads(Threads&& other) noexcept
      : Threads(static_cast<const Threads&>(other)) {}
  KOKKOS_DEFAULTED_FUNCTION Threads& operator=(const Threads&) = default;
  KOKKOS_FUNCTION Threads& operator=(Threads&& other) noexcept {
    return *this = static_cast<const Threads&>(other);
  }

  ~Threads() { Impl::check_execution_space_destructor_precondition(name()); }

  Threads() { Impl::check_execution_space_constructor_precondition(name()); }

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

#endif /* #if defined( KOKKOS_ENABLE_THREADS ) */
#endif /* #define KOKKOS_THREADS_HPP */
