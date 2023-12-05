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
#ifndef KOKKOS_SYCL_HPP
#define KOKKOS_SYCL_HPP

#include <Kokkos_Macros.hpp>

#ifdef KOKKOS_ENABLE_SYCL
// FIXME_SYCL
#if __has_include(<sycl/sycl.hpp>)
#include <sycl/sycl.hpp>
#else
#include <CL/sycl.hpp>
#endif
#include <SYCL/Kokkos_SYCL_Space.hpp>
#include <Kokkos_Layout.hpp>
#include <Kokkos_ScratchSpace.hpp>
#include <impl/Kokkos_Profiling_Interface.hpp>
#include <impl/Kokkos_HostSharedPtr.hpp>
#include <impl/Kokkos_InitializationSettings.hpp>

namespace Kokkos {
namespace Experimental {
namespace Impl {
class SYCLInternal;
}

/// \class SYCL
/// \brief Kokkos device for multicore processors in the host memory space.
class SYCL {
 public:
  //------------------------------------
  //! \name Type declarations that all Kokkos devices must provide.
  //@{

  //! Tag this class as a kokkos execution space
  using execution_space = SYCL;
  using memory_space    = SYCLDeviceUSMSpace;
  using device_type     = Kokkos::Device<execution_space, memory_space>;

  using array_layout = LayoutLeft;
  using size_type    = memory_space::size_type;

  using scratch_memory_space = ScratchMemorySpace<SYCL>;

  SYCL();
  explicit SYCL(const sycl::queue&);

  uint32_t impl_instance_id() const noexcept {
    return m_space_instance->impl_get_instance_id();
  }

  sycl::queue& sycl_queue() const noexcept {
    return *m_space_instance->m_queue;
  }

  //@}
  //------------------------------------
  //! \name Functions that all Kokkos devices must implement.
  //@{

  KOKKOS_INLINE_FUNCTION static int in_parallel() {
#if defined(__SYCL_DEVICE_ONLY__)
    return true;
#else
    return false;
#endif
  }

  /** \brief  Set the device in a "sleep" state. */
  static bool sleep();

  /** \brief Wake the device from the 'sleep' state. A noop for OpenMP. */
  static bool wake();

  /** \brief Wait until all dispatched functors complete. A noop for OpenMP. */
  static void impl_static_fence(const std::string& name);

  void fence(
      const std::string& name =
          "Kokkos::Experimental::SYCL::fence: Unnamed Instance Fence") const;

  /// \brief Print configuration information to the given output stream.
  void print_configuration(std::ostream& os, bool verbose = false) const;

  /// \brief Free any resources being consumed by the device.
  static void impl_finalize();

  static void impl_initialize(InitializationSettings const&);

  static bool impl_is_initialized();

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  static int concurrency();
#else
  int concurrency() const;
#endif

  static const char* name();

  inline Impl::SYCLInternal* impl_internal_space_instance() const {
    return m_space_instance.get();
  }

 private:
  static std::ostream& impl_sycl_info(std::ostream& os,
                                      const sycl::device& device);

  friend bool operator==(SYCL const& lhs, SYCL const& rhs) {
    return lhs.impl_internal_space_instance() ==
           rhs.impl_internal_space_instance();
  }
  friend bool operator!=(SYCL const& lhs, SYCL const& rhs) {
    return !(lhs == rhs);
  }
  Kokkos::Impl::HostSharedPtr<Impl::SYCLInternal> m_space_instance;
};

}  // namespace Experimental

namespace Tools {
namespace Experimental {
template <>
struct DeviceTypeTraits<Kokkos::Experimental::SYCL> {
  /// \brief An ID to differentiate (for example) Serial from OpenMP in Tooling
  static constexpr DeviceType id = DeviceType::SYCL;
  static int device_id(const Kokkos::Experimental::SYCL& exec) {
    return exec.impl_internal_space_instance()->m_syclDev;
  }
};
}  // namespace Experimental
}  // namespace Tools

namespace Experimental {
template <class... Args>
std::vector<SYCL> partition_space(const SYCL& sycl_space, Args...) {
  static_assert(
      (... && std::is_arithmetic_v<Args>),
      "Kokkos Error: partitioning arguments must be integers or floats");

  sycl::context context = sycl_space.sycl_queue().get_context();
  sycl::device device =
      sycl_space.impl_internal_space_instance()->m_queue->get_device();
  std::vector<SYCL> instances;
  instances.reserve(sizeof...(Args));
  for (unsigned int i = 0; i < sizeof...(Args); ++i)
    instances.emplace_back(
        sycl::queue(context, device, sycl::property::queue::in_order()));
  return instances;
}

template <class T>
std::vector<SYCL> partition_space(const SYCL& sycl_space,
                                  std::vector<T> const& weights) {
  static_assert(
      std::is_arithmetic<T>::value,
      "Kokkos Error: partitioning arguments must be integers or floats");

  sycl::context context = sycl_space.sycl_queue().get_context();
  sycl::device device =
      sycl_space.impl_internal_space_instance()->m_queue->get_device();
  std::vector<SYCL> instances;

  // We only care about the number of instances to create and ignore weights
  // otherwise.
  instances.reserve(weights.size());
  for (unsigned int i = 0; i < weights.size(); ++i)
    instances.emplace_back(
        sycl::queue(context, device, sycl::property::queue::in_order()));
  return instances;
}
}  // namespace Experimental

}  // namespace Kokkos

#endif
#endif
