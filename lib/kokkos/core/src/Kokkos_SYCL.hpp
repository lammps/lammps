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

#ifndef KOKKOS_SYCL_HPP
#define KOKKOS_SYCL_HPP

#include <Kokkos_Macros.hpp>

#ifdef KOKKOS_ENABLE_SYCL
#include <CL/sycl.hpp>
#include <Kokkos_SYCL_Space.hpp>
#include <Kokkos_Layout.hpp>
#include <Kokkos_ScratchSpace.hpp>
#include <impl/Kokkos_ExecSpaceInitializer.hpp>
#include <impl/Kokkos_Profiling_Interface.hpp>
#include <impl/Kokkos_HostSharedPtr.hpp>

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

  sycl::context sycl_context() const noexcept {
    return m_space_instance->m_queue->get_context();
  };

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
  static void impl_static_fence();
  static void impl_static_fence(const std::string&);
  void fence() const;
  void fence(const std::string&) const;

  /// \brief Print configuration information to the given output stream.
  void print_configuration(std::ostream&, const bool detail = false);

  /// \brief Free any resources being consumed by the device.
  static void impl_finalize();

  /** \brief  Initialize the device.
   *
   */

  struct SYCLDevice {
    SYCLDevice() : SYCLDevice(sycl::default_selector()) {}
    explicit SYCLDevice(sycl::device d);
    explicit SYCLDevice(const sycl::device_selector& selector);
    explicit SYCLDevice(size_t id);

    sycl::device get_device() const;

    friend std::ostream& operator<<(std::ostream& os, const SYCLDevice& that) {
      return SYCL::impl_sycl_info(os, that.m_device);
    }

   private:
    sycl::device m_device;
  };

  static void impl_initialize(SYCLDevice = SYCLDevice());

  int sycl_device() const;

  static bool impl_is_initialized();

  static int concurrency();
  static const char* name();

  inline Impl::SYCLInternal* impl_internal_space_instance() const {
    return m_space_instance.get();
  }

 private:
  static std::ostream& impl_sycl_info(std::ostream& os,
                                      const sycl::device& device);

  Kokkos::Impl::HostSharedPtr<Impl::SYCLInternal> m_space_instance;
};

namespace Impl {

class SYCLSpaceInitializer : public Kokkos::Impl::ExecSpaceInitializerBase {
 public:
  void initialize(const InitArguments& args) final;
  void finalize(const bool) final;
  void fence() final;
  void fence(const std::string&) final;
  void print_configuration(std::ostream& msg, const bool detail) final;
};

}  // namespace Impl
}  // namespace Experimental

namespace Tools {
namespace Experimental {
template <>
struct DeviceTypeTraits<Kokkos::Experimental::SYCL> {
  /// \brief An ID to differentiate (for example) Serial from OpenMP in Tooling
  static constexpr DeviceType id = DeviceType::SYCL;
  static int device_id(const Kokkos::Experimental::SYCL& exec) {
    return exec.sycl_device();
  }
};
}  // namespace Experimental
}  // namespace Tools

namespace Experimental {
template <class... Args>
std::vector<SYCL> partition_space(const SYCL& sycl_space, Args...) {
#ifdef __cpp_fold_expressions
  static_assert(
      (... && std::is_arithmetic_v<Args>),
      "Kokkos Error: partitioning arguments must be integers or floats");
#endif

  sycl::context context = sycl_space.sycl_context();
  sycl::default_selector device_selector;
  std::vector<SYCL> instances;
  instances.reserve(sizeof...(Args));
  for (unsigned int i = 0; i < sizeof...(Args); ++i)
    instances.emplace_back(sycl::queue(context, device_selector));
  return instances;
}

template <class T>
std::vector<SYCL> partition_space(const SYCL& sycl_space,
                                  std::vector<T>& weights) {
  static_assert(
      std::is_arithmetic<T>::value,
      "Kokkos Error: partitioning arguments must be integers or floats");

  sycl::context context = sycl_space.sycl_context();
  sycl::default_selector device_selector;
  std::vector<SYCL> instances;
  instances.reserve(weights.size());
  for (unsigned int i = 0; i < weights.size(); ++i)
    instances.emplace_back(sycl::queue(context, device_selector));
  return instances;
}
}  // namespace Experimental

}  // namespace Kokkos

#endif
#endif
