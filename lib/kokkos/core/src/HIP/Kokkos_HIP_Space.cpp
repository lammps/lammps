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

#include <Kokkos_Macros.hpp>

#include <Kokkos_Core.hpp>
#include <HIP/Kokkos_HIP_Space.hpp>

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

HIPSpace::HIPSpace() : m_device(HIP().hip_device()) {}

HIPHostPinnedSpace::HIPHostPinnedSpace() {}

HIPManagedSpace::HIPManagedSpace() : m_device(HIP().hip_device()) {}

void* HIPSpace::allocate(const size_t arg_alloc_size) const {
  return allocate("[unlabeled]", arg_alloc_size);
}
void* HIPSpace::allocate(

    const char* arg_label, const size_t arg_alloc_size,
    const size_t arg_logical_size) const {
  return impl_allocate(arg_label, arg_alloc_size, arg_logical_size);
}
void* HIPSpace::impl_allocate(
    const char* arg_label, const size_t arg_alloc_size,
    const size_t arg_logical_size,
    const Kokkos::Tools::SpaceHandle arg_handle) const {
  void* ptr = nullptr;

  auto const error_code = hipMalloc(&ptr, arg_alloc_size);
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
        const char* hsa_xnack = std::getenv("HSA_XNACK");
        if (!hsa_xnack)
          std::cerr << R"warning(
Kokkos::HIP::runtime WARNING: Kokkos did not find an environment variable 'HSA_XNACK'
                              for the current process.
                              Nevertheless, xnack is enabled for all processes if
                              amdgpu.noretry=0 was set in the Linux kernel boot line.
                              Without xnack enabled, Kokkos::HIPManaged might not behave
                              as expected.)warning"
                    << std::endl;
        else if (Kokkos::Impl::strcmp(hsa_xnack, "1") != 0)
          std::cerr
              << "Kokkos::HIP::runtime WARNING: Kokkos detected the "
                 "environement variable "
              << "'HSA_XNACK'=" << hsa_xnack << "\n"
              << "Kokkos advises to set it to '1' to enable it per process."
              << std::endl;
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
  KOKKOS_IMPL_HIP_SAFE_CALL(hipFree(arg_alloc_ptr));
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
  KOKKOS_IMPL_HIP_SAFE_CALL(hipFree(arg_alloc_ptr));
}

}  // namespace Kokkos
