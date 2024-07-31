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

#include <Kokkos_Core.hpp>
#include <HIP/Kokkos_HIP.hpp>
#include <HIP/Kokkos_HIP_Instance.hpp>

#include <impl/Kokkos_DeviceManagement.hpp>
#include <impl/Kokkos_ExecSpaceManager.hpp>

#include <hip/hip_runtime_api.h>

namespace Kokkos {

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
int HIP::concurrency() {
#else
int HIP::concurrency() const {
#endif
  return Impl::HIPInternal::concurrency();
}

int HIP::impl_is_initialized() {
  return Impl::HIPInternal::singleton().is_initialized();
}

void HIP::impl_initialize(InitializationSettings const& settings) {
  const std::vector<int>& visible_devices = Impl::get_visible_devices();
  const int hip_device_id =
      Impl::get_gpu(settings).value_or(visible_devices[0]);

  Impl::HIPInternal::m_hipDev = hip_device_id;
  KOKKOS_IMPL_HIP_SAFE_CALL(
      hipGetDeviceProperties(&Impl::HIPInternal::m_deviceProp, hip_device_id));
  const auto& hipProp = Impl::HIPInternal::m_deviceProp;
  KOKKOS_IMPL_HIP_SAFE_CALL(hipSetDevice(hip_device_id));

  // number of multiprocessors
  Impl::HIPInternal::m_multiProcCount = hipProp.multiProcessorCount;

  //----------------------------------
  // Maximum number of warps,
  // at most one warp per thread in a warp for reduction.
  Impl::HIPInternal::m_maxWarpCount =
      hipProp.maxThreadsPerBlock / Impl::HIPTraits::WarpSize;
  if (Impl::HIPTraits::WarpSize < Impl::HIPInternal::m_maxWarpCount) {
    Impl::HIPInternal::m_maxWarpCount = Impl::HIPTraits::WarpSize;
  }

  //----------------------------------
  // Maximum number of blocks
  Impl::HIPInternal::m_maxBlock[0] = hipProp.maxGridSize[0];
  Impl::HIPInternal::m_maxBlock[1] = hipProp.maxGridSize[1];
  Impl::HIPInternal::m_maxBlock[2] = hipProp.maxGridSize[2];

  // theoretically, we can get 40 WF's / CU, but only can sustain 32 see
  // https://github.com/ROCm-Developer-Tools/HIP/blob/a0b5dfd625d99af7e288629747b40dd057183173/vdi/hip_platform.cpp#L742
  Impl::HIPInternal::m_maxWavesPerCU = 32;
  Impl::HIPInternal::m_shmemPerSM    = hipProp.maxSharedMemoryPerMultiProcessor;
  Impl::HIPInternal::m_maxShmemPerBlock = hipProp.sharedMemPerBlock;
  Impl::HIPInternal::m_maxThreadsPerSM =
      Impl::HIPInternal::m_maxWavesPerCU * Impl::HIPTraits::WarpSize;

  // Init the array for used for arbitrarily sized atomics
  desul::Impl::init_lock_arrays();  // FIXME

  // Allocate a staging buffer for constant mem in pinned host memory
  // and an event to avoid overwriting driver for previous kernel launches
  KOKKOS_IMPL_HIP_SAFE_CALL(
      hipHostMalloc((void**)&Impl::HIPInternal::constantMemHostStaging,
                    Impl::HIPTraits::ConstantMemoryUsage));

  KOKKOS_IMPL_HIP_SAFE_CALL(
      hipEventCreate(&Impl::HIPInternal::constantMemReusable));

  hipStream_t singleton_stream;
  KOKKOS_IMPL_HIP_SAFE_CALL(hipStreamCreate(&singleton_stream));
  Impl::HIPInternal::singleton().initialize(singleton_stream);
}

void HIP::impl_finalize() {
  (void)Impl::hip_global_unique_token_locks(true);

  desul::Impl::finalize_lock_arrays();  // FIXME

  KOKKOS_IMPL_HIP_SAFE_CALL(
      hipEventDestroy(Impl::HIPInternal::constantMemReusable));
  KOKKOS_IMPL_HIP_SAFE_CALL(
      hipHostFree(Impl::HIPInternal::constantMemHostStaging));

  Impl::HIPInternal::singleton().finalize();
  KOKKOS_IMPL_HIP_SAFE_CALL(
      hipStreamDestroy(Impl::HIPInternal::singleton().m_stream));
}

HIP::HIP()
    : m_space_instance(&Impl::HIPInternal::singleton(),
                       [](Impl::HIPInternal*) {}) {
  Impl::HIPInternal::singleton().verify_is_initialized(
      "HIP instance constructor");
}

HIP::HIP(hipStream_t const stream, Impl::ManageStream manage_stream)
    : m_space_instance(
          new Impl::HIPInternal, [manage_stream](Impl::HIPInternal* ptr) {
            ptr->finalize();
            if (static_cast<bool>(manage_stream)) {
              KOKKOS_IMPL_HIP_SAFE_CALL(hipStreamDestroy(ptr->m_stream));
            }
            delete ptr;
          }) {
  Impl::HIPInternal::singleton().verify_is_initialized(
      "HIP instance constructor");
  m_space_instance->initialize(stream);
}

KOKKOS_DEPRECATED HIP::HIP(hipStream_t const stream, bool manage_stream)
    : HIP(stream,
          manage_stream ? Impl::ManageStream::yes : Impl::ManageStream::no) {}

void HIP::print_configuration(std::ostream& os, bool /*verbose*/) const {
  os << "Device Execution Space:\n";
  os << "  KOKKOS_ENABLE_HIP: yes\n";

  os << "HIP Options:\n";
  os << "  KOKKOS_ENABLE_HIP_RELOCATABLE_DEVICE_CODE: ";
#ifdef KOKKOS_ENABLE_HIP_RELOCATABLE_DEVICE_CODE
  os << "yes\n";
#else
  os << "no\n";
#endif
#ifdef KOKKOS_ENABLE_IMPL_HIP_UNIFIED_MEMORY
  os << "  KOKKOS_ENABLE_IMPL_HIP_UNIFIED_MEMORY: ";
  os << "yes\n";
#endif

  os << "\nRuntime Configuration:\n";

  m_space_instance->print_configuration(os);
}

uint32_t HIP::impl_instance_id() const noexcept {
  return m_space_instance->impl_get_instance_id();
}
void HIP::impl_static_fence(const std::string& name) {
  Kokkos::Tools::Experimental::Impl::profile_fence_event<HIP>(
      name,
      Kokkos::Tools::Experimental::SpecialSynchronizationCases::
          GlobalDeviceSynchronization,
      [&]() { KOKKOS_IMPL_HIP_SAFE_CALL(hipDeviceSynchronize()); });
}

void HIP::fence(const std::string& name) const {
  m_space_instance->fence(name);
}

hipStream_t HIP::hip_stream() const { return m_space_instance->m_stream; }

int HIP::hip_device() const { return impl_internal_space_instance()->m_hipDev; }

hipDeviceProp_t const& HIP::hip_device_prop() {
  return Impl::HIPInternal::singleton().m_deviceProp;
}

const char* HIP::name() { return "HIP"; }

namespace Impl {

int g_hip_space_factory_initialized = initialize_space_factory<HIP>("150_HIP");

}  // namespace Impl

}  // namespace Kokkos
