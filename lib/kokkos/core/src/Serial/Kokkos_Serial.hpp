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

/// \file Kokkos_Serial.hpp
/// \brief Declaration and definition of Kokkos::Serial device.

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#include <Kokkos_Macros.hpp>
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#endif
#ifndef KOKKOS_SERIAL_HPP
#define KOKKOS_SERIAL_HPP

#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ENABLE_SERIAL)

#include <cstddef>
#include <iosfwd>
#include <iterator>
#include <mutex>
#include <thread>
#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_TaskScheduler.hpp>
#include <Kokkos_Layout.hpp>
#include <Kokkos_HostSpace.hpp>
#include <Kokkos_ScratchSpace.hpp>
#include <Kokkos_MemoryTraits.hpp>
#include <impl/Kokkos_HostThreadTeam.hpp>
#include <impl/Kokkos_FunctorAnalysis.hpp>
#include <impl/Kokkos_Tools.hpp>
#include <impl/Kokkos_HostSharedPtr.hpp>
#include <impl/Kokkos_InitializationSettings.hpp>

namespace Kokkos {

namespace Impl {
class SerialInternal {
 public:
  SerialInternal() = default;

  bool is_initialized();

  void initialize();

  void finalize();

  static SerialInternal& singleton();

  std::mutex m_thread_team_data_mutex;

  // Resize thread team data scratch memory
  void resize_thread_team_data(size_t pool_reduce_bytes,
                               size_t team_reduce_bytes,
                               size_t team_shared_bytes,
                               size_t thread_local_bytes);

  HostThreadTeamData m_thread_team_data;
  bool m_is_initialized = false;
};
}  // namespace Impl

struct NewInstance {
  explicit NewInstance() = default;
};

/// \class Serial
/// \brief Kokkos device for non-parallel execution
///
/// A "device" represents a parallel execution model.  It tells Kokkos
/// how to parallelize the execution of kernels in a parallel_for or
/// parallel_reduce.  For example, the Threads device uses
/// C++11 threads on a CPU, the OpenMP device uses the OpenMP language
/// extensions, and the Cuda device uses NVIDIA's CUDA programming
/// model.  The Serial device executes "parallel" kernels
/// sequentially.  This is useful if you really do not want to use
/// threads, or if you want to explore different combinations of MPI
/// and shared-memory parallel programming models.
class Serial {
 public:
  //! \name Type declarations that all Kokkos devices must provide.
  //@{

  //! Tag this class as an execution space:
  using execution_space = Serial;
  //! This device's preferred memory space.
  using memory_space = Kokkos::HostSpace;
  //! The size_type alias best suited for this device.
  using size_type = memory_space::size_type;
  //! This execution space preferred device_type
  using device_type = Kokkos::Device<execution_space, memory_space>;

  //! This device's preferred array layout.
  using array_layout = LayoutRight;

  /// \brief  Scratch memory space
  using scratch_memory_space = ScratchMemorySpace<Kokkos::Serial>;

  //@}

  Serial();

  Serial(NewInstance);

  /// \brief True if and only if this method is being called in a
  ///   thread-parallel function.
  ///
  /// For the Serial device, this method <i>always</i> returns false,
  /// because parallel_for or parallel_reduce with the Serial device
  /// always execute sequentially.

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  KOKKOS_DEPRECATED inline static int in_parallel() { return false; }
#endif

  /// \brief Wait until all dispatched functors complete.
  ///
  /// The parallel_for or parallel_reduce dispatch of a functor may
  /// return asynchronously, before the functor completes.  This
  /// method does not return until all dispatched functors on this
  /// device have completed.
  static void impl_static_fence(const std::string& name) {
    Kokkos::Tools::Experimental::Impl::profile_fence_event<Kokkos::Serial>(
        name,
        Kokkos::Tools::Experimental::SpecialSynchronizationCases::
            GlobalDeviceSynchronization,
        []() {});  // TODO: correct device ID
    Kokkos::memory_fence();
  }

  void fence(const std::string& name =
                 "Kokkos::Serial::fence: Unnamed Instance Fence") const {
    Kokkos::Tools::Experimental::Impl::profile_fence_event<Kokkos::Serial>(
        name, Kokkos::Tools::Experimental::Impl::DirectFenceIDHandle{1},
        []() {});  // TODO: correct device ID
    Kokkos::memory_fence();
  }

  /** \brief  Return the maximum amount of concurrency.  */
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  static int concurrency() { return 1; }
#else
  int concurrency() const { return 1; }
#endif

  //! Print configuration information to the given output stream.
  void print_configuration(std::ostream& os, bool verbose = false) const;

  static void impl_initialize(InitializationSettings const&);

  static bool impl_is_initialized();

  //! Free any resources being consumed by the device.
  static void impl_finalize();

  //--------------------------------------------------------------------------

  inline static int impl_thread_pool_size(int = 0) { return 1; }
  KOKKOS_INLINE_FUNCTION static int impl_thread_pool_rank() { return 0; }

  //--------------------------------------------------------------------------

  KOKKOS_INLINE_FUNCTION static unsigned impl_hardware_thread_id() {
    return impl_thread_pool_rank();
  }
  inline static unsigned impl_max_hardware_threads() {
    return impl_thread_pool_size(0);
  }

  uint32_t impl_instance_id() const noexcept { return 1; }

  static const char* name();

  Impl::SerialInternal* impl_internal_space_instance() const {
    return m_space_instance.get();
  }

 private:
  Kokkos::Impl::HostSharedPtr<Impl::SerialInternal> m_space_instance;
  friend bool operator==(Serial const& lhs, Serial const& rhs) {
    return lhs.impl_internal_space_instance() ==
           rhs.impl_internal_space_instance();
  }
  friend bool operator!=(Serial const& lhs, Serial const& rhs) {
    return !(lhs == rhs);
  }
  //--------------------------------------------------------------------------
};

namespace Tools {
namespace Experimental {
template <>
struct DeviceTypeTraits<Serial> {
  static constexpr DeviceType id = DeviceType::Serial;
  static int device_id(const Serial&) { return 0; }
};
}  // namespace Experimental
}  // namespace Tools
}  // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

template <>
struct MemorySpaceAccess<Kokkos::Serial::memory_space,
                         Kokkos::Serial::scratch_memory_space> {
  enum : bool { assignable = false };
  enum : bool { accessible = true };
  enum : bool { deepcopy = false };
};

}  // namespace Impl
}  // namespace Kokkos

namespace Kokkos::Experimental {

template <class... Args>
std::vector<Serial> partition_space(const Serial&, Args...) {
  static_assert(
      (... && std::is_arithmetic_v<Args>),
      "Kokkos Error: partitioning arguments must be integers or floats");
  std::vector<Serial> instances;
  instances.reserve(sizeof...(Args));
  std::generate_n(std::back_inserter(instances), sizeof...(Args),
                  []() { return Serial{NewInstance{}}; });
  return instances;
}

template <class T>
std::vector<Serial> partition_space(const Serial&,
                                    std::vector<T> const& weights) {
  static_assert(
      std::is_arithmetic<T>::value,
      "Kokkos Error: partitioning arguments must be integers or floats");

  // We only care about the number of instances to create and ignore weights
  // otherwise.
  std::vector<Serial> instances;
  instances.reserve(weights.size());
  std::generate_n(std::back_inserter(instances), weights.size(),
                  []() { return Serial{NewInstance{}}; });
  return instances;
}

}  // namespace Kokkos::Experimental

#include <Serial/Kokkos_Serial_Parallel_Range.hpp>
#include <Serial/Kokkos_Serial_Parallel_MDRange.hpp>
#include <Serial/Kokkos_Serial_Parallel_Team.hpp>
#include <Serial/Kokkos_Serial_Task.hpp>
#include <Serial/Kokkos_Serial_UniqueToken.hpp>

#endif  // defined( KOKKOS_ENABLE_SERIAL )
#endif  /* #define KOKKOS_SERIAL_HPP */
