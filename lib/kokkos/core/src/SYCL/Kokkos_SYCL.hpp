// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

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
    // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
    return *m_space_instance->m_queue;
  }

  //@}
  //------------------------------------
  //! \name Functions that all Kokkos devices must implement.
  //@{

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  KOKKOS_DEPRECATED KOKKOS_INLINE_FUNCTION static int in_parallel() {
#if defined(__SYCL_DEVICE_ONLY__)
    return true;
#else
    return false;
#endif
  }
#endif

  /** \brief Wait until all dispatched functors complete. A noop for OpenMP. */
  static void impl_static_fence(const std::string& name);

  void fence(const std::string& name =
                 "Kokkos::SYCL::fence: Unnamed Instance Fence") const;

  /// \brief Print configuration information to the given output stream.
  void print_configuration(std::ostream& os, bool verbose = false) const;

  /// \brief Free any resources being consumed by the device.
  static void impl_finalize();

  static void impl_initialize(InitializationSettings const&);

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

namespace Tools {
namespace Experimental {
template <>
struct DeviceTypeTraits<Kokkos::SYCL> {
  /// \brief An ID to differentiate (for example) Serial from OpenMP in Tooling
  static constexpr DeviceType id = DeviceType::SYCL;
  static int device_id(const Kokkos::SYCL& exec) {
    return exec.impl_internal_space_instance()->m_syclDev;
  }
};
}  // namespace Experimental
}  // namespace Tools

namespace Experimental::Impl {
// For each space in partition, create new queue on the same device as
// base_instance, ignoring weights
template <class T>
std::vector<SYCL> impl_partition_space(const SYCL& base_instance,
                                       const std::vector<T>& weights) {
  sycl::context context = base_instance.sycl_queue().get_context();
  sycl::device device   = base_instance.sycl_queue().get_device();

  std::vector<SYCL> instances;
  instances.reserve(weights.size());
  std::generate_n(std::back_inserter(instances), weights.size(),
                  [&context, &device]() {
                    return SYCL(sycl::queue(context, device
#ifdef KOKKOS_IMPL_SYCL_USE_IN_ORDER_QUEUES
                                            ,
                                            sycl::property::queue::in_order()
#endif
                                                ));
                  });

  return instances;
}
}  // namespace Experimental::Impl

namespace Impl {
std::vector<sycl::device> get_sycl_devices();
}  // namespace Impl

}  // namespace Kokkos

#endif
#endif
