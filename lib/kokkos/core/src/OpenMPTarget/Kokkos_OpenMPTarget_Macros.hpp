// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_OPENMPTARGET_MACROS_HPP
#define KOKKOS_OPENMPTARGET_MACROS_HPP

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
