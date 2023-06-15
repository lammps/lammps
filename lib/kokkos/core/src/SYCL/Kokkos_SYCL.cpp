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

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#endif

#include <Kokkos_Concepts.hpp>
#include <SYCL/Kokkos_SYCL_Instance.hpp>
#include <Kokkos_SYCL.hpp>
#include <Kokkos_HostSpace.hpp>
#include <Kokkos_Serial.hpp>
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
  Impl::SYCLInternal::singleton().verify_is_initialized(
      "SYCL instance constructor");
  m_space_instance->initialize(stream);
}

int SYCL::concurrency() {
  return Impl::SYCLInternal::singleton().m_maxConcurrency;
}

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

int SYCL::sycl_device() const {
  return impl_internal_space_instance()->m_syclDev;
}

void SYCL::impl_initialize(InitializationSettings const& settings) {
  std::vector<sycl::device> gpu_devices =
      sycl::device::get_devices(sycl::info::device_type::gpu);
  // If the device id is not specified and there are no GPUs, sidestep Kokkos
  // device selection and use whatever is available (if no GPU architecture is
  // specified).
#if !defined(KOKKOS_ARCH_INTEL_GPU) && !defined(KOKKOS_ARCH_KEPLER) && \
    !defined(KOKKOS_ARCH_MAXWELL) && !defined(KOKKOS_ARCH_PASCAL) &&   \
    !defined(KOKKOS_ARCH_VOLTA) && !defined(KOKKOS_ARCH_TURING75) &&   \
    !defined(KOKKOS_ARCH_AMPERE) && !defined(KOKKOS_ARCH_HOPPER)
  if (!settings.has_device_id() && gpu_devices.empty()) {
    Impl::SYCLInternal::singleton().initialize(sycl::device());
    return;
  }
#endif
  using Kokkos::Impl::get_gpu;
  Impl::SYCLInternal::singleton().initialize(gpu_devices[get_gpu(settings)]);
}

std::ostream& SYCL::impl_sycl_info(std::ostream& os,
                                   const sycl::device& device) {
  using namespace sycl::info;
  return os << "Name: " << device.get_info<device::name>()
            << "\nDriver Version: " << device.get_info<device::driver_version>()
            << "\nIs Host: " << device.is_host()
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
            << "\nImage Support: " << device.get_info<device::image_support>()
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
            << "\nHost Unified Memory: "
            << device.get_info<device::host_unified_memory>()
            << "\nProfiling Timer Resolution: "
            << device.get_info<device::profiling_timer_resolution>()
            << "\nIs Endian Little: "
            << device.get_info<device::is_endian_little>()
            << "\nIs Available: " << device.get_info<device::is_available>()
            << "\nIs Compiler Available: "
            << device.get_info<device::is_compiler_available>()
            << "\nIs Linker Available: "
            << device.get_info<device::is_linker_available>()
            << "\nQueue Profiling: "
            << device.get_info<device::queue_profiling>()
            << "\nVendor: " << device.get_info<device::vendor>()
            << "\nProfile: " << device.get_info<device::profile>()
            << "\nVersion: " << device.get_info<device::version>()
            << "\nPrintf Buffer Size: "
            << device.get_info<device::printf_buffer_size>()
            << "\nPreferred Interop User Sync: "
            << device.get_info<device::preferred_interop_user_sync>()
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
