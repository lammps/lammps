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

#ifndef KOKKOS_OPENMPTARGET_MACROS_HPP
#define KOKKOS_OPENMPTARGET_MACROS_HPP

// Intel architectures prefer the classical hierarchical parallelism that relies
// on OpenMP.
#if defined(KOKKOS_ARCH_INTEL_GPU)
#define KOKKOS_IMPL_OPENMPTARGET_HIERARCHICAL_INTEL_GPU
#endif

// Define a macro for llvm compiler greater than version 17 and on NVIDIA and
// AMD GPUs. This would be useful in cases where non-OpenMP standard llvm
// extensions can be used.
#if defined(KOKKOS_COMPILER_CLANG) && (KOKKOS_COMPILER_CLANG >= 1700) && \
    (defined(KOKKOS_ARCH_AMD_GPU) || defined(KOKKOS_IMPL_ARCH_NVIDIA_GPU))
#define KOKKOS_IMPL_OPENMPTARGET_LLVM_EXTENSIONS
#endif

#define KOKKOS_IMPL_OPENMPTARGET_PRAGMA_HELPER(x) _Pragma(#x)
#define KOKKOS_IMPL_OMPTARGET_PRAGMA(x) \
  KOKKOS_IMPL_OPENMPTARGET_PRAGMA_HELPER(omp target x)

// Use scratch memory extensions to request dynamic shared memory for the
// right compiler/architecture combination.
#ifdef KOKKOS_IMPL_OPENMPTARGET_LLVM_EXTENSIONS
#define KOKKOS_IMPL_OMPX_DYN_CGROUP_MEM(N) ompx_dyn_cgroup_mem(N)
#else
#define KOKKOS_IMPL_OMPX_DYN_CGROUP_MEM(N)
#endif

#endif  // KOKKOS_OPENMPTARGET_MACROS_HPP
