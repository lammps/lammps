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

#ifndef KOKKOS_OPENACC_TRAITS_HPP
#define KOKKOS_OPENACC_TRAITS_HPP

#include <openacc.h>

namespace Kokkos::Experimental::Impl {

struct OpenACC_Traits {
#if defined(KOKKOS_IMPL_ARCH_NVIDIA_GPU)
  static constexpr acc_device_t dev_type     = acc_device_nvidia;
  static constexpr bool may_fallback_to_host = false;
#else
  static constexpr acc_device_t dev_type     = acc_device_not_host;
  static constexpr bool may_fallback_to_host = true;
#endif
};

}  // namespace Kokkos::Experimental::Impl

#endif
