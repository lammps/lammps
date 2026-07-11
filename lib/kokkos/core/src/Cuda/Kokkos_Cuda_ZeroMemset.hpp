// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project
#ifndef KOKKOS_CUDA_ZEROMEMSET_HPP
#define KOKKOS_CUDA_ZEROMEMSET_HPP

#include <Kokkos_Macros.hpp>
#include <Cuda/Kokkos_Cuda.hpp>
#include <impl/Kokkos_ZeroMemset_fwd.hpp>

namespace Kokkos {
namespace Impl {

template <>
struct ZeroMemset<Kokkos::Cuda> {
  ZeroMemset(const Kokkos::Cuda& exec_space_instance, void* dst, size_t cnt) {
    KOKKOS_IMPL_CUDA_SAFE_CALL(
        (exec_space_instance.impl_internal_space_instance()
             ->cuda_memset_async_wrapper(dst, 0, cnt)));
  }
};

}  // namespace Impl
}  // namespace Kokkos

#endif  // !defined(KOKKOS_CUDA_ZEROMEMSET_HPP)
