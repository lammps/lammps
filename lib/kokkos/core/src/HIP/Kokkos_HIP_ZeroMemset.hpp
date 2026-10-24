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

template <>
struct ZeroMemset<HIP> {
  ZeroMemset(const HIP& exec_space, void* dst, size_t cnt) {
    KOKKOS_IMPL_HIP_SAFE_CALL(
        exec_space.impl_internal_space_instance()->hip_memset_async_wrapper(
            dst, 0, cnt));
  }
};

}  // namespace Impl
}  // namespace Kokkos

#endif  // !defined(KOKKOS_HIP_ZEROMEMSET_HPP)
