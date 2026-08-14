// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_CUDA_HALF_IMPL_TYPE_HPP_
#define KOKKOS_CUDA_HALF_IMPL_TYPE_HPP_

#include <Kokkos_Macros.hpp>

#if !(defined(KOKKOS_ARCH_MAXWELL50) || defined(KOKKOS_ARCH_MAXWELL52))

#include <cuda_fp16.h>
#include <cuda_bf16.h>

#ifndef KOKKOS_IMPL_HALF_TYPE_DEFINED
// Make sure no one else tries to define half_t
#define KOKKOS_IMPL_HALF_TYPE_DEFINED

namespace Kokkos::Impl {

struct half_impl_t {
  using type = __half;
};
#define KOKKOS_IMPL_BHALF_TYPE_DEFINED
struct bhalf_impl_t {
  using type = __nv_bfloat16;
};

}  // namespace Kokkos::Impl

#endif  // KOKKOS_IMPL_HALF_TYPE_DEFINED
#endif  // Disables for half_t on cuda: MAXWELL50||MAXWELL52

#endif
