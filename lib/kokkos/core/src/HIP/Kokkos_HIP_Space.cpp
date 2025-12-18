// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#endif

#include <Kokkos_Macros.hpp>

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif
#include <HIP/Kokkos_HIP_Space.hpp>
#include <HIP/Kokkos_HIP_IsXnack.hpp>

#include <HIP/Kokkos_HIP_DeepCopy.hpp>

#include <impl/Kokkos_Error.hpp>
#include <impl/Kokkos_DeviceManagement.hpp>
#include <impl/Kokkos_ExecSpaceManager.hpp>

#include <hip/hip_runtime_api.h>

#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <atomic>

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace {

static std::atomic<bool> is_first_hip_managed_allocation(true);

}  // namespace

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {

HIPSpace::HIPSpace()
    : m_device(HIP().hip_device()), m_stream(HIP().hip_stream()) {}
HIPSpace::HIPSpace(int device_id, hipStream_t stream)
    : m_device(device_id), m_stream(stream) {}

HIPHostPinnedSpace::HIPHostPinnedSpace()
    : m_device(HIP().hip_device()), m_stream(HIP().hip_stream()) {}
HIPHostPinnedSpace::HIPHostPinnedSpace(int device_id, hipStream_t stream)
    : m_device(device_id), m_stream(stream) {}

HIPManagedSpace::HIPManagedSpace()
    : m_device(HIP().hip_device()), m_stream(HIP().hip_stream()) {}
HIPManagedSpace::HIPManagedSpace(int device_id, hipStream_t stream)
    : m_device(device_id), m_stream(stream) {}

void* HIPSpace::allocate(const HIP& exec_space,
                         const size_t arg_alloc_size) const {
  return allocate(exec_space, "[unlabeled]", arg_alloc_size);
}

void* HIPSpace::allocate(const HIP& exec_space, const char* arg_label,
                         const size_t arg_alloc_size,
                         const size_t arg_logical_size) const {
  return impl_allocate(exec_space.hip_device(), exec_space.hip_stream(),
                       arg_label, arg_alloc_size, arg_logical_size, true);
}

void* HIPSpace::allocate(const size_t arg_alloc_size) const {
  return allocate("[unlabeled]", arg_alloc_size);
}

void* HIPSpace::allocate(const char* arg_label, const size_t arg_alloc_size,
                         const size_t arg_logical_size) const {
  return impl_allocate(m_device, m_stream, arg_label, arg_alloc_size,
                       arg_logical_size, false);
}

void* HIPSpace::impl_allocate(const int device_id,
                              [[maybe_unused]] const hipStream_t stream,
                              const char* arg_label,
                              const size_t arg_alloc_size,
                              const size_t arg_logical_size,
                              [[maybe_unused]] bool stream_sync_only) const {
  void* ptr = nullptr;
  // Instead of trying to allocate zero memory, return early.
  if (arg_alloc_size == 0) return ptr;

  KOKKOS_IMPL_HIP_SAFE_CALL(hipSetDevice(device_id));

#ifdef KOKKOS_ENABLE_IMPL_HIP_MALLOC_ASYNC
  auto const error_code = hipMallocAsync(&ptr, arg_alloc_size, stream);
  if (stream_sync_only) {
    KOKKOS_IMPL_HIP_SAFE_CALL(hipStreamSynchronize(stream));
  } else {
    KOKKOS_IMPL_HIP_SAFE_CALL(hipDeviceSynchronize());
  }
#else
  auto const error_code = hipMalloc(&ptr, arg_alloc_size);
#endif

  if (error_code != hipSuccess) {
    // This is the only way to clear the last error, which we should do here
    // since we're turning it into an exception here
    (void)hipGetLastError();
    Kokkos::Impl::throw_bad_alloc(name(), arg_alloc_size, arg_label);
  }
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    const auto arg_handle = Kokkos::Tools::make_space_handle(name());
    const size_t reported_size =
        (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;
    Kokkos::Profiling::allocateData(arg_handle, arg_label, ptr, reported_size);
  }

  return ptr;
}

void* HIPHostPinnedSpace::allocate(const size_t arg_alloc_size) const {
  return allocate("[unlabeled]", arg_alloc_size);
}
void* HIPHostPinnedSpace::allocate(const char* arg_label,
                                   const size_t arg_alloc_size,
                                   const size_t arg_logical_size) const {
  return impl_allocate(arg_label, arg_alloc_size, arg_logical_size);
}
void* HIPHostPinnedSpace::impl_allocate(
    const char* arg_label, const size_t arg_alloc_size,
    const size_t arg_logical_size,
    const Kokkos::Tools::SpaceHandle arg_handle) const {
  void* ptr = nullptr;

  KOKKOS_IMPL_HIP_SAFE_CALL(hipSetDevice(m_device));
  auto const error_code =
      hipHostMalloc(&ptr, arg_alloc_size, hipHostMallocNonCoherent);
  if (error_code != hipSuccess) {
    // This is the only way to clear the last error, which we should do here
    // since we're turning it into an exception here
    (void)hipGetLastError();
    Kokkos::Impl::throw_bad_alloc(name(), arg_alloc_size, arg_label);
  }
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    const size_t reported_size =
        (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;
    Kokkos::Profiling::allocateData(arg_handle, arg_label, ptr, reported_size);
  }

  return ptr;
}

void* HIPManagedSpace::allocate(const size_t arg_alloc_size) const {
  return allocate("[unlabeled]", arg_alloc_size);
}
void* HIPManagedSpace::allocate(const char* arg_label,
                                const size_t arg_alloc_size,
                                const size_t arg_logical_size) const {
  return impl_allocate(arg_label, arg_alloc_size, arg_logical_size);
}
void* HIPManagedSpace::impl_allocate(
    const char* arg_label, const size_t arg_alloc_size,
    const size_t arg_logical_size,
    const Kokkos::Tools::SpaceHandle arg_handle) const {
  void* ptr = nullptr;

  if (arg_alloc_size > 0) {
    KOKKOS_IMPL_HIP_SAFE_CALL(hipSetDevice(m_device));
    if (is_first_hip_managed_allocation.exchange(false) &&
        Kokkos::show_warnings()) {
      do {  // hack to avoid spamming users with too many warnings
        if (!impl_hip_driver_check_page_migration()) {
          std::cerr << R"warning(
Kokkos::HIP::allocation WARNING: The combination of device and system configuration
                                 does not support page migration between device and host.
                                 HIPManagedSpace might not work as expected.
                                 Please refer to the ROCm documentation on unified/managed memory.)warning"
                    << std::endl;
          break;  // do not warn about HSA_XNACK environement variable
        }

        // check for correct runtime environment
        if (!Kokkos::Impl::xnack_environment_enabled())
          std::cerr << R"warning(
Kokkos::HIP::runtime WARNING: Kokkos was not able to verify that xnack is enabled.
                              Without xnack enabled, Kokkos::HIPManaged might not behave as expected.
                              Set HSA_XNACK=1 in your environment. For further information on HMM support
                              call `Kokkos::print_configuration`, or run with KOKKOS_PRINT_CONFIGURATION=1
                              in your environment.
)warning";
      } while (false);
    }
    auto const error_code = hipMallocManaged(&ptr, arg_alloc_size);
    if (error_code != hipSuccess) {
      // This is the only way to clear the last error, which we should do here
      // since we're turning it into an exception here
      (void)hipGetLastError();
      Kokkos::Impl::throw_bad_alloc(name(), arg_alloc_size, arg_label);
    }
    KOKKOS_IMPL_HIP_SAFE_CALL(hipMemAdvise(
        ptr, arg_alloc_size, hipMemAdviseSetCoarseGrain, m_device));
  }

  if (Kokkos::Profiling::profileLibraryLoaded()) {
    const size_t reported_size =
        (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;
    Kokkos::Profiling::allocateData(arg_handle, arg_label, ptr, reported_size);
  }

  return ptr;
}
bool HIPManagedSpace::impl_hip_driver_check_page_migration() const {
  // check with driver if page migrating memory is available
  // this driver query is copied from the hip documentation
  int hasManagedMemory = 0;  // false by default
  KOKKOS_IMPL_HIP_SAFE_CALL(hipDeviceGetAttribute(
      &hasManagedMemory, hipDeviceAttributeManagedMemory, m_device));
  if (!static_cast<bool>(hasManagedMemory)) return false;
  // next, check pageableMemoryAccess
  int hasPageableMemory = 0;  // false by default
  KOKKOS_IMPL_HIP_SAFE_CALL(hipDeviceGetAttribute(
      &hasPageableMemory, hipDeviceAttributePageableMemoryAccess, m_device));
  return static_cast<bool>(hasPageableMemory);
}

void HIPSpace::deallocate(void* const arg_alloc_ptr,
                          const size_t arg_alloc_size) const {
  deallocate("[unlabeled]", arg_alloc_ptr, arg_alloc_size);
}
void HIPSpace::deallocate(const char* arg_label, void* const arg_alloc_ptr,
                          const size_t arg_alloc_size,
                          const size_t arg_logical_size) const {
  impl_deallocate(arg_label, arg_alloc_ptr, arg_alloc_size, arg_logical_size);
}
void HIPSpace::impl_deallocate(
    const char* arg_label, void* const arg_alloc_ptr,
    const size_t arg_alloc_size, const size_t arg_logical_size,
    const Kokkos::Tools::SpaceHandle arg_handle) const {
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    const size_t reported_size =
        (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;
    Kokkos::Profiling::deallocateData(arg_handle, arg_label, arg_alloc_ptr,
                                      reported_size);
  }
#ifdef KOKKOS_ENABLE_IMPL_HIP_MALLOC_ASYNC
  KOKKOS_IMPL_HIP_SAFE_CALL(hipSetDevice(m_device));
  KOKKOS_IMPL_HIP_SAFE_CALL(hipFreeAsync(arg_alloc_ptr, m_stream));
  KOKKOS_IMPL_HIP_SAFE_CALL(hipDeviceSynchronize());
#else
  KOKKOS_IMPL_HIP_SAFE_CALL(hipSetDevice(m_device));
  KOKKOS_IMPL_HIP_SAFE_CALL(hipFree(arg_alloc_ptr));
#endif
}

void HIPHostPinnedSpace::deallocate(void* const arg_alloc_ptr,
                                    const size_t arg_alloc_size) const {
  deallocate("[unlabeled]", arg_alloc_ptr, arg_alloc_size);
}

void HIPHostPinnedSpace::deallocate(const char* arg_label,
                                    void* const arg_alloc_ptr,
                                    const size_t arg_alloc_size,
                                    const size_t arg_logical_size) const {
  impl_deallocate(arg_label, arg_alloc_ptr, arg_alloc_size, arg_logical_size);
}
void HIPHostPinnedSpace::impl_deallocate(
    const char* arg_label, void* const arg_alloc_ptr,
    const size_t arg_alloc_size, const size_t arg_logical_size,
    const Kokkos::Tools::SpaceHandle arg_handle) const {
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    const size_t reported_size =
        (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;
    Kokkos::Profiling::deallocateData(arg_handle, arg_label, arg_alloc_ptr,
                                      reported_size);
  }
  KOKKOS_IMPL_HIP_SAFE_CALL(hipSetDevice(m_device));
  KOKKOS_IMPL_HIP_SAFE_CALL(hipHostFree(arg_alloc_ptr));
}

void HIPManagedSpace::deallocate(void* const arg_alloc_ptr,
                                 const size_t arg_alloc_size) const {
  deallocate("[unlabeled]", arg_alloc_ptr, arg_alloc_size);
}

void HIPManagedSpace::deallocate(const char* arg_label,
                                 void* const arg_alloc_ptr,
                                 const size_t arg_alloc_size,
                                 const size_t arg_logical_size) const {
  impl_deallocate(arg_label, arg_alloc_ptr, arg_alloc_size, arg_logical_size);
}
void HIPManagedSpace::impl_deallocate(
    const char* arg_label, void* const arg_alloc_ptr,
    const size_t arg_alloc_size, const size_t arg_logical_size,
    const Kokkos::Tools::SpaceHandle arg_handle) const {
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    const size_t reported_size =
        (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;
    Kokkos::Profiling::deallocateData(arg_handle, arg_label, arg_alloc_ptr,
                                      reported_size);
  }
  // We have to unset the CoarseGrain property manually as hipFree does not take
  // care of it. Otherwise, the allocation would continue to linger in the
  // kernel mem page table.
  KOKKOS_IMPL_HIP_SAFE_CALL(hipMemAdvise(
      arg_alloc_ptr, arg_alloc_size, hipMemAdviseUnsetCoarseGrain, m_device));
  KOKKOS_IMPL_HIP_SAFE_CALL(hipSetDevice(m_device));
  KOKKOS_IMPL_HIP_SAFE_CALL(hipFree(arg_alloc_ptr));
}

}  // namespace Kokkos
