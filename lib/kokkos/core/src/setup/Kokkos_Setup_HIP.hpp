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

#ifndef KOKKOS_SETUP_HIP_HPP_
#define KOKKOS_SETUP_HIP_HPP_

#if defined(KOKKOS_ENABLE_HIP)

#define KOKKOS_IMPL_HIP_CLANG_WORKAROUND

#include <hip/hip_runtime.h>
#include <hip/hip_runtime_api.h>

#define KOKKOS_LAMBDA [=] __host__ __device__
#define KOKKOS_CLASS_LAMBDA [ =, *this ] __host__ __device__

#define KOKKOS_IMPL_FORCEINLINE_FUNCTION __device__ __host__ __forceinline__
#define KOKKOS_IMPL_INLINE_FUNCTION __device__ __host__ inline
#define KOKKOS_IMPL_FUNCTION __device__ __host__
#define KOKKOS_IMPL_HOST_FUNCTION __host__
#define KOKKOS_IMPL_DEVICE_FUNCTION __device__

#endif  // #if defined( KOKKOS_ENABLE_HIP )

#endif
