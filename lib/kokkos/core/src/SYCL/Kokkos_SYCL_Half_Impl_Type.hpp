// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_SYCL_HALF_IMPL_TYPE_HPP_
#define KOKKOS_SYCL_HALF_IMPL_TYPE_HPP_

#include <Kokkos_Macros.hpp>

#include <sycl/sycl.hpp>

namespace Kokkos::Impl {
struct half_impl_t {
  using type = sycl::half;
};
}  // namespace Kokkos::Impl

// FIXME_SYCL Evaluate when to drop the check
#if __has_include(<sycl/ext/oneapi/bfloat16.hpp>)
namespace Kokkos::Impl {
struct bhalf_impl_t {
  using type = sycl::ext::oneapi::bfloat16;
};
}  // namespace Kokkos::Impl
#elif defined(SYCL_EXT_ONEAPI_BFLOAT16) && defined(KOKKOS_ARCH_INTEL_GPU)
// FIXME_SYCL bfloat16 is only supported for compute capability 8.0 or higher
// on Nvidia GPUs but SYCL_EXT_ONEAPI_BFLOAT16 is defined even for lower compute
// capability.
namespace Kokkos::Impl {
struct bhalf_impl_t {
  using type = sycl::ext::oneapi::experimental::bfloat16;
};
}  // namespace Kokkos::Impl
#endif  // test for bfloat16 support
#endif  // KOKKOS_SYCL_HALF_IMPL_TYPE_HPP_
