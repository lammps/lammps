// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOS_HIP_ZEROMEMSET_HPP
#define KOKKOS_HIP_ZEROMEMSET_HPP

#include <Kokkos_Macros.hpp>
#include <HIP/Kokkos_HIP.hpp>
#include <HIP/Kokkos_HIP_Instance.hpp>
#include <impl/Kokkos_ZeroMemset_fwd.hpp>

namespace Kokkos {
namespace Impl {

// hipMemsetAsync sets the first `cnt` bytes of `dst` to the provided value
void zero_with_hip_kernel(const HIP& exec_space, void* dst, size_t cnt);

template <>
struct ZeroMemset<HIP> {
  ZeroMemset(const HIP& exec_space, void* dst, size_t cnt) {
    // We allow user on an AMD APU with unified memory to `malloc` and wrap that
    // in an unmanaged SharedSpace view. In ROCm <= 6.2.1 (and possibly later),
    // hipMemsetAsync on a host-allocated pointer returns an invalid value
    // error, but accessing the data via a GPU kernel works as long as xnack is
    // present and enabled (HSA_XNACK=1)
#if defined(KOKKOS_IMPL_HIP_UNIFIED_MEMORY)
    zero_with_hip_kernel(exec_space, dst, cnt);
#else
    KOKKOS_IMPL_HIP_SAFE_CALL(
        exec_space.impl_internal_space_instance()->hip_memset_async_wrapper(
            dst, 0, cnt));
#endif
  }
};

}  // namespace Impl
}  // namespace Kokkos

#endif  // !defined(KOKKOS_HIP_ZEROMEMSET_HPP)
