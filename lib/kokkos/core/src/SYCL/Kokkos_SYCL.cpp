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
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#endif

#include <Kokkos_Concepts.hpp>
#include <SYCL/Kokkos_SYCL_Instance.hpp>
#include <SYCL/Kokkos_SYCL.hpp>
#include <Kokkos_HostSpace.hpp>
#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Error.hpp>
#include <impl/Kokkos_DeviceManagement.hpp>
#include <impl/Kokkos_ExecSpaceManager.hpp>

namespace {
template <typename C>
struct Container {
  explicit Container(const C& c) : container(c) {}

  friend std::ostream& operator<<(std::ostream& os, const Container& that) {
    os << that.container.size();
    for (const auto& v : that.container) {
      os << "\n\t" << v;
    }
    return os;
  }

 private:
  const C& container;
};
}  // namespace

namespace Kokkos {
namespace Experimental {
SYCL::SYCL()
    : m_space_instance(&Impl::SYCLInternal::singleton(),
                       [](Impl::SYCLInternal*) {}) {
  Impl::SYCLInternal::singleton().verify_is_initialized(
      "SYCL instance constructor");
}

SYCL::SYCL(const sycl::queue& stream)
    : m_space_instance(new Impl::SYCLInternal, [](Impl::SYCLInternal* ptr) {
        ptr->finalize();
        delete ptr;
      }) {
  // In principle could be guarded with
  // #ifdef KOKKOS_IMPL_SYCL_USE_IN_ORDER_QUEUES
  // but we chose to require user-provided queues to be in-order
  // unconditionally so that code downstream does not break
  // when the backend setting changes.
  if (!stream.is_in_order())
    Kokkos::abort("User provided sycl::queues must be in-order!");
  Impl::SYCLInternal::singleton().verify_is_initialized(
      "SYCL instance constructor");
  m_space_instance->initialize(stream);
}

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
int SYCL::concurrency() {
  return Impl::SYCLInternal::singleton().m_maxConcurrency;
}
#else
int SYCL::concurrency() const { return m_space_instance->m_maxConcurrency; }
#endif

const char* SYCL::name() { return "SYCL"; }

bool SYCL::impl_is_initialized() {
  return Impl::SYCLInternal::singleton().is_initialized();
}

void SYCL::impl_finalize() { Impl::SYCLInternal::singleton().finalize(); }

void SYCL::print_configuration(std::ostream& os, bool verbose) const {
  os << "Devices:\n";
  os << "  KOKKOS_ENABLE_SYCL: yes\n";

  os << "\nRuntime Configuration:\n";

  os << "macro  KOKKOS_ENABLE_SYCL : defined\n";
#ifdef KOKKOS_IMPL_SYCL_DEVICE_GLOBAL_SUPPORTED
  os << "macro  KOKKOS_IMPL_SYCL_DEVICE_GLOBAL_SUPPORTED : defined\n";
#else
  os << "macro  KOKKOS_IMPL_SYCL_DEVICE_GLOBAL_SUPPORTED : undefined\n";
#endif

#ifdef KOKKOS_IMPL_SYCL_USE_IN_ORDER_QUEUES
  os << "macro  KOKKOS_IMPL_SYCL_USE_IN_ORDER_QUEUES : defined\n";
#else
  os << "macro  KOKKOS_IMPL_SYCL_USE_IN_ORDER_QUEUES : undefined\n";
#endif

  if (verbose)
    SYCL::impl_sycl_info(os, m_space_instance->m_queue->get_device());
}

void SYCL::fence(const std::string& name) const {
  Impl::SYCLInternal::fence(*m_space_instance->m_queue, name,
                            impl_instance_id());
}

void SYCL::impl_static_fence(const std::string& name) {
  Kokkos::Tools::Experimental::Impl::profile_fence_event<
      Kokkos::Experimental::SYCL>(
      name,
      Kokkos::Tools::Experimental::SpecialSynchronizationCases::
          GlobalDeviceSynchronization,
      [&]() {
        // guard accessing all_queues
        std::scoped_lock lock(Impl::SYCLInternal::mutex);
        for (auto& queue : Impl::SYCLInternal::all_queues) {
          try {
            (*queue)->wait_and_throw();
          } catch (sycl::exception const& e) {
            Kokkos::Impl::throw_runtime_exception(
                std::string("There was a synchronous SYCL error:\n") +=
                e.what());
          }
        }
      });
}

void SYCL::impl_initialize(InitializationSettings const& settings) {
  std::vector<sycl::device> gpu_devices =
      sycl::device::get_devices(sycl::info::device_type::gpu);
  // If the device id is not specified and there are no GPUs, sidestep Kokkos
  // device selection and use whatever is available (if no GPU architecture is
  // specified).
#if !defined(KOKKOS_ARCH_INTEL_GPU) && !defined(KOKKOS_IMPL_ARCH_NVIDIA_GPU)
  if (!settings.has_device_id() && gpu_devices.empty()) {
    Impl::SYCLInternal::singleton().initialize(sycl::device());
    Impl::SYCLInternal::m_syclDev = 0;
    return;
  }
#endif
  const auto id = ::Kokkos::Impl::get_gpu(settings);
  Impl::SYCLInternal::singleton().initialize(gpu_devices[id]);
  Impl::SYCLInternal::m_syclDev = id;
}

std::ostream& SYCL::impl_sycl_info(std::ostream& os,
                                   const sycl::device& device) {
  using namespace sycl::info;
  return os << "Name: " << device.get_info<device::name>()
            << "\nDriver Version: " << device.get_info<device::driver_version>()
            << "\nIs CPU: " << device.is_cpu()
            << "\nIs GPU: " << device.is_gpu()
            << "\nIs Accelerator: " << device.is_accelerator()
            << "\nVendor Id: " << device.get_info<device::vendor_id>()
            << "\nMax Compute Units: "
            << device.get_info<device::max_compute_units>()
            << "\nMax Work Item Dimensions: "
            << device.get_info<device::max_work_item_dimensions>()
            << "\nMax Work Group Size: "
            << device.get_info<device::max_work_group_size>()
            << "\nPreferred Vector Width Char: "
            << device.get_info<device::preferred_vector_width_char>()
            << "\nPreferred Vector Width Short: "
            << device.get_info<device::preferred_vector_width_short>()
            << "\nPreferred Vector Width Int: "
            << device.get_info<device::preferred_vector_width_int>()
            << "\nPreferred Vector Width Long: "
            << device.get_info<device::preferred_vector_width_long>()
            << "\nPreferred Vector Width Float: "
            << device.get_info<device::preferred_vector_width_float>()
            << "\nPreferred Vector Width Double: "
            << device.get_info<device::preferred_vector_width_double>()
            << "\nPreferred Vector Width Half: "
            << device.get_info<device::preferred_vector_width_half>()
            << "\nNative Vector Width Char: "
            << device.get_info<device::native_vector_width_char>()
            << "\nNative Vector Width Short: "
            << device.get_info<device::native_vector_width_short>()
            << "\nNative Vector Width Int: "
            << device.get_info<device::native_vector_width_int>()
            << "\nNative Vector Width Long: "
            << device.get_info<device::native_vector_width_long>()
            << "\nNative Vector Width Float: "
            << device.get_info<device::native_vector_width_float>()
            << "\nNative Vector Width Double: "
            << device.get_info<device::native_vector_width_double>()
            << "\nNative Vector Width Half: "
            << device.get_info<device::native_vector_width_half>()
            << "\nAddress Bits: " << device.get_info<device::address_bits>()
            << "\nMax Mem Alloc Size: "
            << device.get_info<device::max_mem_alloc_size>()
            << "\nMax Read Image Args: "
            << device.get_info<device::max_read_image_args>()
            << "\nImage2d Max Width: "
            << device.get_info<device::image2d_max_width>()
            << "\nImage2d Max Height: "
            << device.get_info<device::image2d_max_height>()
            << "\nImage3d Max Width: "
            << device.get_info<device::image3d_max_width>()
            << "\nImage3d Max Height: "
            << device.get_info<device::image3d_max_height>()
            << "\nImage3d Max Depth: "
            << device.get_info<device::image3d_max_depth>()
            << "\nImage Max Buffer Size: "
            << device.get_info<device::image_max_buffer_size>()
            << "\nImage Max Array Size: "
            << device.get_info<device::image_max_array_size>()
            << "\nMax Samplers: " << device.get_info<device::max_samplers>()
            << "\nMax Parameter Size: "
            << device.get_info<device::max_parameter_size>()
            << "\nMem Base Addr Align: "
            << device.get_info<device::mem_base_addr_align>()
            << "\nGlobal Cache Mem Line Size: "
            << device.get_info<device::global_mem_cache_line_size>()
            << "\nGlobal Mem Cache Size: "
            << device.get_info<device::global_mem_cache_size>()
            << "\nGlobal Mem Size: "
            << device.get_info<device::global_mem_size>()
            << "\nLocal Mem Size: " << device.get_info<device::local_mem_size>()
            << "\nError Correction Support: "
            << device.get_info<device::error_correction_support>()
            << "\nProfiling Timer Resolution: "
            << device.get_info<device::profiling_timer_resolution>()
            << "\nIs Available: " << device.get_info<device::is_available>()
            << "\nVendor: " << device.get_info<device::vendor>()
            << "\nVersion: " << device.get_info<device::version>()
            << "\nPartition Max Sub Devices: "
            << device.get_info<device::partition_max_sub_devices>()
            << "\nReference Count: "
            << device.get_info<device::reference_count>() << '\n';
}

namespace Impl {

int g_sycl_space_factory_initialized =
    Kokkos::Impl::initialize_space_factory<SYCL>("170_SYCL");

}
}  // namespace Experimental
}  // namespace Kokkos
