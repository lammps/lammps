// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_HIP_HPP
#define KOKKOS_HIP_HPP

#include <Kokkos_Core_fwd.hpp>

#include <Kokkos_Layout.hpp>
#include <HIP/Kokkos_HIP_Space.hpp>

#include <hip/hip_runtime_api.h>

namespace Kokkos {
namespace Impl {
class HIPInternal;
enum class ManageStream : bool { no, yes };
}  // namespace Impl
/// \class HIP
/// \brief Kokkos device for multicore processors in the host memory space.
class HIP {
 public:
  //------------------------------------
  //! \name Type declarations that all Kokkos devices must provide.
  //@{

  //! Tag this class as a kokkos execution space
  using execution_space = HIP;
  using memory_space    = HIPSpace;
  using device_type     = Kokkos::Device<execution_space, memory_space>;

  using array_layout = LayoutLeft;
  using size_type    = HIPSpace::size_type;

  using scratch_memory_space = ScratchMemorySpace<HIP>;

  KOKKOS_DEFAULTED_FUNCTION HIP(const HIP&) = default;
  KOKKOS_FUNCTION HIP(HIP&& other) noexcept
      : HIP(static_cast<const HIP&>(other)) {}
  KOKKOS_DEFAULTED_FUNCTION HIP& operator=(const HIP&) = default;
  KOKKOS_FUNCTION HIP& operator=(HIP&& other) noexcept {
    return *this = static_cast<const HIP&>(other);
  }
  ~HIP();
  HIP();

  explicit HIP(hipStream_t stream) : HIP(stream, Impl::ManageStream::no) {}

  HIP(hipStream_t stream, Impl::ManageStream manage_stream);

  //@}
  //------------------------------------
  //! \name Functions that all Kokkos devices must implement.
  //@{

  /** \brief Wait until all dispatched functors complete.
   *
   * The parallel_for or parallel_reduce dispatch of a functor may return
   * asynchronously, before the functor completes. This method does not return
   * until all dispatched functors on this device have completed.
   */
  static void impl_static_fence(const std::string& name);

  void fence(const std::string& name =
                 "Kokkos::HIP::fence(): Unnamed Instance Fence") const;

  hipStream_t hip_stream() const;

  /// \brief Print configuration information to the given output stream.
  void print_configuration(std::ostream& os, bool verbose = false) const;

  /// \brief Free any resources being consumed by the device.
  static void impl_finalize();

  /** \brief  Initialize the device.
   *
   */
  int hip_device() const;
  static hipDeviceProp_t const& hip_device_prop();

  static void impl_initialize(InitializationSettings const&);

  int concurrency() const;

  static const char* name();

  inline Impl::HIPInternal* impl_internal_space_instance() const {
    return m_space_instance.get();
  }

  uint32_t impl_instance_id() const noexcept;

 private:
  friend bool operator==(HIP const& lhs, HIP const& rhs) {
    return lhs.impl_internal_space_instance() ==
           rhs.impl_internal_space_instance();
  }
  friend bool operator!=(HIP const& lhs, HIP const& rhs) {
    return !(lhs == rhs);
  }
  Kokkos::Impl::HostSharedPtr<Impl::HIPInternal> m_space_instance;
};

namespace Impl {
template <>
struct MemorySpaceAccess<HIPSpace, HIP::scratch_memory_space> {
  enum : bool { assignable = false };
  enum : bool { accessible = true };
  enum : bool { deepcopy = false };
};
}  // namespace Impl

namespace Tools {
namespace Experimental {
template <>
struct DeviceTypeTraits<HIP> {
  static constexpr DeviceType id = DeviceType::HIP;
  static int device_id(const HIP& exec) { return exec.hip_device(); }
};
}  // namespace Experimental
}  // namespace Tools
}  // namespace Kokkos

#endif
