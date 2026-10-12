// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_NEXTSILICON_DEEP_COPY_HPP
#define KOKKOS_NEXTSILICON_DEEP_COPY_HPP

#include <nextapi/memory.h>

#include <NextSilicon/Kokkos_NextSilicon.hpp>
#include <NextSilicon/Kokkos_NextSiliconSpace.hpp>

#include <Kokkos_Concepts.hpp>

namespace Kokkos {
namespace Impl {

template <>
struct DeepCopy<Kokkos::Experimental::NextSiliconSharedSpace,
                Kokkos::Experimental::NextSiliconSharedSpace,
                Kokkos::Experimental::NextSilicon> {
  DeepCopy(void* dst, const void* src, size_t n) {
    // Let NextSilicon runtime select the correct implementation.
    nextapi_memory_copy(dst, src, n);
  }
  DeepCopy(const Kokkos::Experimental::NextSilicon&, void* dst, const void* src,
           size_t n) {
    // Let NextSilicon runtime select the correct implementation.
    nextapi_memory_copy(dst, src, n);
  }
};

template <class ExecutionSpace>
struct DeepCopy<Kokkos::Experimental::NextSiliconSharedSpace,
                Kokkos::Experimental::NextSiliconSharedSpace, ExecutionSpace> {
  DeepCopy(void* dst, const void* src, size_t n) {
    // Let NextSilicon runtime select the correct implementation.
    nextapi_memory_copy(dst, src, n);
  }
  DeepCopy(const ExecutionSpace& exec, void* dst, const void* src, size_t n) {
    exec.fence(
        "Kokkos::Impl::DeepCopy<NextSiliconSharedSpace, "
        "NextSiliconSharedSpace, "
        "ExecutionSpace>::DeepCopy: fence before copy");
    // Let NextSilicon runtime select the correct implementation.
    nextapi_memory_copy(dst, src, n);
  }
};

template <>
struct DeepCopy<Kokkos::Experimental::NextSiliconSharedSpace, Kokkos::HostSpace,
                Kokkos::Experimental::NextSilicon> {
  DeepCopy(void* dst, const void* src, size_t n) {
    // Let NextSilicon runtime select the correct implementation.
    nextapi_memory_copy(dst, src, n);
  }
  DeepCopy(const Kokkos::Experimental::NextSilicon&, void* dst, const void* src,
           size_t n) {
    // Let NextSilicon runtime select the correct implementation.
    nextapi_memory_copy(dst, src, n);
  }
};

template <class ExecutionSpace>
struct DeepCopy<Kokkos::Experimental::NextSiliconSharedSpace, Kokkos::HostSpace,
                ExecutionSpace> {
  DeepCopy(void* dst, const void* src, size_t n) {
    // Let NextSilicon runtime select the correct implementation.
    nextapi_memory_copy(dst, src, n);
  }
  DeepCopy(const ExecutionSpace& exec, void* dst, const void* src, size_t n) {
    exec.fence(
        "Kokkos::Impl::DeepCopy<NextSiliconSharedSpace, HostSpace, "
        "ExecutionSpace>::DeepCopy: fence before copy");
    // Let NextSilicon runtime select the correct implementation.
    nextapi_memory_copy(dst, src, n);
  }
};

template <>
struct DeepCopy<Kokkos::HostSpace, Kokkos::Experimental::NextSiliconSharedSpace,
                Kokkos::Experimental::NextSilicon> {
  DeepCopy(void* dst, const void* src, size_t n) {
    // Let NextSilicon runtime select the correct implementation.
    nextapi_memory_copy(dst, src, n);
  }
  DeepCopy(const Kokkos::Experimental::NextSilicon&, void* dst, const void* src,
           size_t n) {
    // Let NextSilicon runtime select the correct implementation.
    nextapi_memory_copy(dst, src, n);
  }
};

template <class ExecutionSpace>
struct DeepCopy<Kokkos::HostSpace, Kokkos::Experimental::NextSiliconSharedSpace,
                ExecutionSpace> {
  DeepCopy(void* dst, const void* src, size_t n) {
    // Let NextSilicon runtime select the correct implementation.
    nextapi_memory_copy(dst, src, n);
  }
  DeepCopy(const ExecutionSpace& exec, void* dst, const void* src, size_t n) {
    exec.fence(
        "Kokkos::Impl::DeepCopy<HostSpace, NextSiliconSharedSpace, "
        "ExecutionSpace>::DeepCopy: fence before copy");
    // Let NextSilicon runtime select the correct implementation.
    nextapi_memory_copy(dst, src, n);
  }
};

}  // namespace Impl
}  // namespace Kokkos

#endif
