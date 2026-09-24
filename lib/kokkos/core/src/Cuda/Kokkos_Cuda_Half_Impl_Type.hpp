// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_CUDA_HALF_IMPL_TYPE_HPP_
#define KOKKOS_CUDA_HALF_IMPL_TYPE_HPP_

#include <Kokkos_Macros.hpp>

#if !(defined(KOKKOS_ARCH_MAXWELL50) || defined(KOKKOS_ARCH_MAXWELL52))

#include <cuda_fp16.h>
#include <cuda_bf16.h>

namespace Kokkos::Impl {

struct half_impl_t {
  using type = __half;
};
struct bhalf_impl_t {
  using type = __nv_bfloat16;
};

}  // namespace Kokkos::Impl

#endif  // Disables for half_t on cuda: MAXWELL50||MAXWELL52

#endif
