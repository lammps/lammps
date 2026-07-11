// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_HIP_HALF_IMPL_TYPE_HPP_
#define KOKKOS_HIP_HALF_IMPL_TYPE_HPP_

#include <hip/hip_fp16.h>
#include <hip/hip_bf16.h>

#ifndef KOKKOS_IMPL_HALF_TYPE_DEFINED
// Make sure no one else tries to define half_t
#define KOKKOS_IMPL_HALF_TYPE_DEFINED

namespace Kokkos {
namespace Impl {
struct half_impl_t {
  using type = __half;
};
}  // namespace Impl
}  // namespace Kokkos
#endif  // KOKKOS_IMPL_HALF_TYPE_DEFINED

#ifndef KOKKOS_IMPL_BHALF_TYPE_DEFINED
// Make sure no one else tries to define bhalf_t
#define KOKKOS_IMPL_BHALF_TYPE_DEFINED
namespace Kokkos::Impl {
struct bhalf_impl_t {
  using type = __hip_bfloat16;
};

}  // namespace Kokkos::Impl

#endif  // KOKKOS_IMPL_BHALF_TYPE_DEFINED

#endif  // KOKKOS_ENABLE_HIP
