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

#ifndef KOKKOS_HIP_DEEP_COPY_HPP
#define KOKKOS_HIP_DEEP_COPY_HPP

#include <HIP/Kokkos_HIP_Space.hpp>
#include <HIP/Kokkos_HIP_Error.hpp>  // HIP_SAFE_CALL

#include <hip/hip_runtime_api.h>

namespace Kokkos {
namespace Impl {

void DeepCopyHIP(void* dst, const void* src, size_t n);
void DeepCopyAsyncHIP(const HIP& instance, void* dst, const void* src,
                      size_t n);
void DeepCopyAsyncHIP(void* dst, const void* src, size_t n);

template <class MemSpace>
struct DeepCopy<MemSpace, HostSpace, HIP,
                std::enable_if_t<is_hip_type_space<MemSpace>::value>> {
  DeepCopy(void* dst, const void* src, size_t n) { DeepCopyHIP(dst, src, n); }
  DeepCopy(const HIP& instance, void* dst, const void* src, size_t n) {
    DeepCopyAsyncHIP(instance, dst, src, n);
  }
};

template <class MemSpace>
struct DeepCopy<HostSpace, MemSpace, HIP,
                std::enable_if_t<is_hip_type_space<MemSpace>::value>> {
  DeepCopy(void* dst, const void* src, size_t n) { DeepCopyHIP(dst, src, n); }
  DeepCopy(const HIP& instance, void* dst, const void* src, size_t n) {
    DeepCopyAsyncHIP(instance, dst, src, n);
  }
};

template <class MemSpace1, class MemSpace2>
struct DeepCopy<MemSpace1, MemSpace2, HIP,
                std::enable_if_t<is_hip_type_space<MemSpace1>::value &&
                                 is_hip_type_space<MemSpace2>::value>> {
  DeepCopy(void* dst, const void* src, size_t n) { DeepCopyHIP(dst, src, n); }
  DeepCopy(const HIP& instance, void* dst, const void* src, size_t n) {
    DeepCopyAsyncHIP(instance, dst, src, n);
  }
};

template <class MemSpace1, class MemSpace2, class ExecutionSpace>
struct DeepCopy<MemSpace1, MemSpace2, ExecutionSpace,
                std::enable_if_t<is_hip_type_space<MemSpace1>::value &&
                                 is_hip_type_space<MemSpace2>::value &&
                                 !std::is_same<ExecutionSpace, HIP>::value>> {
  inline DeepCopy(void* dst, const void* src, size_t n) {
    DeepCopyHIP(dst, src, n);
  }

  inline DeepCopy(const ExecutionSpace& exec, void* dst, const void* src,
                  size_t n) {
    exec.fence(fence_string());
    DeepCopyAsyncHIP(dst, src, n);
  }

 private:
  static const std::string& fence_string() {
    static const std::string string =
        std::string("Kokkos::Impl::DeepCopy<") + MemSpace1::name() + "Space, " +
        MemSpace2::name() +
        "Space, ExecutionSpace>::DeepCopy: fence before copy";
    return string;
  }
};

template <class MemSpace, class ExecutionSpace>
struct DeepCopy<MemSpace, HostSpace, ExecutionSpace,
                std::enable_if_t<is_hip_type_space<MemSpace>::value &&
                                 !std::is_same<ExecutionSpace, HIP>::value>> {
  inline DeepCopy(void* dst, const void* src, size_t n) {
    DeepCopyHIP(dst, src, n);
  }

  inline DeepCopy(const ExecutionSpace& exec, void* dst, const void* src,
                  size_t n) {
    exec.fence(fence_string());
    DeepCopyAsyncHIP(dst, src, n);
  }

 private:
  static const std::string& fence_string() {
    static const std::string string =
        std::string("Kokkos::Impl::DeepCopy<") + MemSpace::name() +
        "Space, HostSpace, ExecutionSpace>::DeepCopy: fence before copy";
    return string;
  }
};

template <class MemSpace, class ExecutionSpace>
struct DeepCopy<HostSpace, MemSpace, ExecutionSpace,
                std::enable_if_t<is_hip_type_space<MemSpace>::value &&
                                 !std::is_same<ExecutionSpace, HIP>::value>> {
  inline DeepCopy(void* dst, const void* src, size_t n) {
    DeepCopyHIP(dst, src, n);
  }

  inline DeepCopy(const ExecutionSpace& exec, void* dst, const void* src,
                  size_t n) {
    exec.fence(fence_string());
    DeepCopyAsyncHIP(dst, src, n);
  }

 private:
  static const std::string& fence_string() {
    static const std::string string =
        std::string("Kokkos::Impl::DeepCopy<HostSpace, ") + MemSpace::name() +
        "Space, ExecutionSpace>::DeepCopy: fence before copy";
    return string;
  }
};
}  // namespace Impl
}  // namespace Kokkos

#endif
