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

#include <HIP/Kokkos_HIP_DeepCopy.hpp>
#include <HIP/Kokkos_HIP_Error.hpp>  // HIP_SAFE_CALL
#include <HIP/Kokkos_HIP.hpp>

namespace Kokkos {
namespace Impl {
namespace {
hipStream_t get_deep_copy_stream() {
  static hipStream_t s = nullptr;
  if (s == nullptr) {
    KOKKOS_IMPL_HIP_SAFE_CALL(hipStreamCreate(&s));
  }
  return s;
}
}  // namespace

void DeepCopyHIP(void* dst, void const* src, size_t n) {
  KOKKOS_IMPL_HIP_SAFE_CALL(hipMemcpyAsync(dst, src, n, hipMemcpyDefault));
}

void DeepCopyAsyncHIP(const HIP& instance, void* dst, void const* src,
                      size_t n) {
  KOKKOS_IMPL_HIP_SAFE_CALL(
      hipMemcpyAsync(dst, src, n, hipMemcpyDefault, instance.hip_stream()));
}

void DeepCopyAsyncHIP(void* dst, void const* src, size_t n) {
  hipStream_t s = get_deep_copy_stream();
  KOKKOS_IMPL_HIP_SAFE_CALL(hipMemcpyAsync(dst, src, n, hipMemcpyDefault, s));
  Kokkos::Tools::Experimental::Impl::profile_fence_event<HIP>(
      "Kokkos::Impl::DeepCopyAsyncHIP: Post Deep Copy Fence on Deep-Copy "
      "stream",
      Kokkos::Tools::Experimental::SpecialSynchronizationCases::
          DeepCopyResourceSynchronization,
      [&]() { KOKKOS_IMPL_HIP_SAFE_CALL(hipStreamSynchronize(s)); });
}
}  // namespace Impl
}  // namespace Kokkos
