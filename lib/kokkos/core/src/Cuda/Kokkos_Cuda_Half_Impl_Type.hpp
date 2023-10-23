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

#ifndef KOKKOS_CUDA_HALF_IMPL_TYPE_HPP_
#define KOKKOS_CUDA_HALF_IMPL_TYPE_HPP_

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_CUDA
#if !(defined(KOKKOS_COMPILER_CLANG) && KOKKOS_COMPILER_CLANG < 900) && \
    !(defined(KOKKOS_ARCH_KEPLER) || defined(KOKKOS_ARCH_MAXWELL50) ||  \
      defined(KOKKOS_ARCH_MAXWELL52))
#include <cuda_fp16.h>
#if (CUDA_VERSION >= 11000)
#include <cuda_bf16.h>
#endif  // CUDA_VERSION >= 11000

#ifndef KOKKOS_IMPL_HALF_TYPE_DEFINED
// Make sure no one else tries to define half_t
#define KOKKOS_IMPL_HALF_TYPE_DEFINED
#define KOKKOS_IMPL_CUDA_HALF_TYPE_DEFINED

namespace Kokkos {
namespace Impl {
struct half_impl_t {
  using type = __half;
};
#if (CUDA_VERSION >= 11000)
#define KOKKOS_IMPL_BHALF_TYPE_DEFINED
struct bhalf_impl_t {
  using type = __nv_bfloat16;
};
#endif  // CUDA_VERSION >= 11000
}  // namespace Impl
}  // namespace Kokkos
#endif  // KOKKOS_IMPL_HALF_TYPE_DEFINED
#endif  // Disables for half_t on cuda:
        // Clang/8||KEPLER30||KEPLER32||KEPLER37||MAXWELL50||MAXWELL52
#endif  // KOKKOS_ENABLE_CUDA
#endif
