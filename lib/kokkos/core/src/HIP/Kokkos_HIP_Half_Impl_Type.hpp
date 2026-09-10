// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_HIP_HALF_IMPL_TYPE_HPP_
#define KOKKOS_HIP_HALF_IMPL_TYPE_HPP_

#include <hip/hip_fp16.h>
#include <hip/hip_bf16.h>

namespace Kokkos {
namespace Impl {
struct half_impl_t {
  using type = __half;
};
}  // namespace Impl
}  // namespace Kokkos

namespace Kokkos::Impl {
struct bhalf_impl_t {
  using type = __hip_bfloat16;
};

}  // namespace Kokkos::Impl

#endif  // KOKKOS_ENABLE_HIP
