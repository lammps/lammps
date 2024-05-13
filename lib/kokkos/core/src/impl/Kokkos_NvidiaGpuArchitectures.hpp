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

#ifndef KOKKOS_CUDA_NVIDIA_GPU_ARCHITECTURES_HPP
#define KOKKOS_CUDA_NVIDIA_GPU_ARCHITECTURES_HPP

#if defined(KOKKOS_ARCH_KEPLER30)
#define KOKKOS_IMPL_ARCH_NVIDIA_GPU 30
#elif defined(KOKKOS_ARCH_KEPLER32)
#define KOKKOS_IMPL_ARCH_NVIDIA_GPU 32
#elif defined(KOKKOS_ARCH_KEPLER35)
#define KOKKOS_IMPL_ARCH_NVIDIA_GPU 35
#elif defined(KOKKOS_ARCH_KEPLER37)
#define KOKKOS_IMPL_ARCH_NVIDIA_GPU 37
#elif defined(KOKKOS_ARCH_MAXWELL50)
#define KOKKOS_IMPL_ARCH_NVIDIA_GPU 50
#elif defined(KOKKOS_ARCH_MAXWELL52)
#define KOKKOS_IMPL_ARCH_NVIDIA_GPU 52
#elif defined(KOKKOS_ARCH_MAXWELL53)
#define KOKKOS_IMPL_ARCH_NVIDIA_GPU 53
#elif defined(KOKKOS_ARCH_PASCAL60)
#define KOKKOS_IMPL_ARCH_NVIDIA_GPU 60
#elif defined(KOKKOS_ARCH_PASCAL61)
#define KOKKOS_IMPL_ARCH_NVIDIA_GPU 61
#elif defined(KOKKOS_ARCH_VOLTA70)
#define KOKKOS_IMPL_ARCH_NVIDIA_GPU 70
#elif defined(KOKKOS_ARCH_VOLTA72)
#define KOKKOS_IMPL_ARCH_NVIDIA_GPU 72
#elif defined(KOKKOS_ARCH_TURING75)
#define KOKKOS_IMPL_ARCH_NVIDIA_GPU 75
#elif defined(KOKKOS_ARCH_AMPERE80)
#define KOKKOS_IMPL_ARCH_NVIDIA_GPU 80
#elif defined(KOKKOS_ARCH_AMPERE86)
#define KOKKOS_IMPL_ARCH_NVIDIA_GPU 86
#elif defined(KOKKOS_ARCH_ADA89)
#define KOKKOS_IMPL_ARCH_NVIDIA_GPU 89
#elif defined(KOKKOS_ARCH_HOPPER90)
#define KOKKOS_IMPL_ARCH_NVIDIA_GPU 90
#elif defined(KOKKOS_ENABLE_CUDA)
// do not raise an error on other backends that may run on NVIDIA GPUs such as
// OpenACC, OpenMPTarget, or SYCL
#error NVIDIA GPU arch not recognized
#endif

#endif
