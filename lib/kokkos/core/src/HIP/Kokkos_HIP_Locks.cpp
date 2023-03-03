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

#include <HIP/Kokkos_HIP_Locks.hpp>
#include <HIP/Kokkos_HIP_Error.hpp>
#include <HIP/Kokkos_HIP.hpp>
#include <HIP/Kokkos_HIP_Instance.hpp>

#include <hip/hip_runtime.h>

#include <iostream>

namespace Kokkos {

#ifdef KOKKOS_ENABLE_HIP_RELOCATABLE_DEVICE_CODE
namespace Impl {
__device__ __constant__ HIPLockArrays g_device_hip_lock_arrays = {nullptr, 0};
}
#endif

namespace {

__global__ void init_lock_array_kernel_atomic() {
  unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < KOKKOS_IMPL_HIP_SPACE_ATOMIC_MASK + 1) {
    Kokkos::Impl::g_device_hip_lock_arrays.atomic[i] = 0;
  }
}

}  // namespace

namespace Impl {

HIPLockArrays g_host_hip_lock_arrays = {nullptr, 0};

void initialize_host_hip_lock_arrays() {
#ifdef KOKKOS_ENABLE_IMPL_DESUL_ATOMICS
  desul::Impl::init_lock_arrays();

  DESUL_ENSURE_HIP_LOCK_ARRAYS_ON_DEVICE();
#endif

  if (g_host_hip_lock_arrays.atomic != nullptr) return;
  KOKKOS_IMPL_HIP_SAFE_CALL(hipMalloc(
      &g_host_hip_lock_arrays.atomic,
      sizeof(std::int32_t) * (KOKKOS_IMPL_HIP_SPACE_ATOMIC_MASK + 1)));

  g_host_hip_lock_arrays.n = HIPInternal::concurrency();

  KOKKOS_COPY_HIP_LOCK_ARRAYS_TO_DEVICE();
  init_lock_array_kernel_atomic<<<
      (KOKKOS_IMPL_HIP_SPACE_ATOMIC_MASK + 1 + 255) / 256, 256, 0, nullptr>>>();
}

void finalize_host_hip_lock_arrays() {
#ifdef KOKKOS_ENABLE_IMPL_DESUL_ATOMICS
  desul::Impl::finalize_lock_arrays();
#endif

  if (g_host_hip_lock_arrays.atomic == nullptr) return;
  KOKKOS_IMPL_HIP_SAFE_CALL(hipFree(g_host_hip_lock_arrays.atomic));
  g_host_hip_lock_arrays.atomic = nullptr;
  g_host_hip_lock_arrays.n      = 0;
#ifdef KOKKOS_ENABLE_HIP_RELOCATABLE_DEVICE_CODE
  KOKKOS_COPY_HIP_LOCK_ARRAYS_TO_DEVICE();
#endif
}

}  // namespace Impl

}  // namespace Kokkos
