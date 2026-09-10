// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_OPENACC_TRAITS_HPP
#define KOKKOS_OPENACC_TRAITS_HPP

#include <openacc.h>

namespace Kokkos::Experimental::Impl {

struct OpenACC_Traits {
#if defined(KOKKOS_IMPL_ARCH_NVIDIA_GPU)
  static constexpr acc_device_t dev_type     = acc_device_nvidia;
  static constexpr bool may_fallback_to_host = false;
#elif defined(KOKKOS_ARCH_AMD_GPU)
  static constexpr acc_device_t dev_type     = acc_device_radeon;
  static constexpr bool may_fallback_to_host = false;
#elif defined(KOKKOS_ENABLE_OPENACC_FORCE_HOST_AS_DEVICE)
  static constexpr acc_device_t dev_type     = acc_device_host;
  static constexpr bool may_fallback_to_host = true;
#else
  static constexpr acc_device_t dev_type     = acc_device_default;
  static constexpr bool may_fallback_to_host = true;
#endif
};

}  // namespace Kokkos::Experimental::Impl

#endif
