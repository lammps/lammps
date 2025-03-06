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
#ifndef KOKKOS_HIP_ZEROMEMSET_HPP
#define KOKKOS_HIP_ZEROMEMSET_HPP

#include <Kokkos_Macros.hpp>
#include <HIP/Kokkos_HIP.hpp>
#include <impl/Kokkos_ZeroMemset_fwd.hpp>

namespace Kokkos {
namespace Impl {

// hipMemsetAsync sets the first `cnt` bytes of `dst` to the provided value
void zero_with_hip_kernel(const HIP& exec_space, void* dst, size_t cnt);

template <>
struct ZeroMemset<HIP> {
  ZeroMemset(const HIP& exec_space, void* dst, size_t cnt) {
    // in ROCm <= 6.2.0, hipMemsetAsync on a host-allocated pointer
    // returns an invalid value error, but accessing the data via a
    // GPU kernel works.
#if defined(KOKKOS_IMPL_HIP_UNIFIED_MEMORY)
    zero_with_hip_kernel(exec_space, dst, cnt);
#else
    KOKKOS_IMPL_HIP_SAFE_CALL(
        hipMemsetAsync(dst, 0, cnt, exec_space.hip_stream()));
#endif
  }
};

}  // namespace Impl
}  // namespace Kokkos

#endif  // !defined(KOKKOS_HIP_ZEROMEMSET_HPP)
