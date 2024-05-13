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

#ifndef KOKKOS_SYCLDEEPCOPY_HPP
#define KOKKOS_SYCLDEEPCOPY_HPP

#include <Kokkos_Core_fwd.hpp>
#include <SYCL/Kokkos_SYCL.hpp>

#include <vector>

#ifdef KOKKOS_ENABLE_SYCL

namespace Kokkos {
namespace Impl {

void DeepCopySYCL(void* dst, const void* src, size_t n);
void DeepCopyAsyncSYCL(const Kokkos::Experimental::SYCL& instance, void* dst,
                       const void* src, size_t n);
void DeepCopyAsyncSYCL(void* dst, const void* src, size_t n);

template <class MemSpace>
struct DeepCopy<MemSpace, HostSpace, Kokkos::Experimental::SYCL,
                std::enable_if_t<is_sycl_type_space<MemSpace>::value>> {
  DeepCopy(void* dst, const void* src, size_t n) { DeepCopySYCL(dst, src, n); }
  DeepCopy(const Kokkos::Experimental::SYCL& instance, void* dst,
           const void* src, size_t n) {
    DeepCopyAsyncSYCL(instance, dst, src, n);
  }
};

template <class MemSpace>
struct DeepCopy<HostSpace, MemSpace, Kokkos::Experimental::SYCL,
                std::enable_if_t<is_sycl_type_space<MemSpace>::value>> {
  DeepCopy(void* dst, const void* src, size_t n) { DeepCopySYCL(dst, src, n); }
  DeepCopy(const Kokkos::Experimental::SYCL& instance, void* dst,
           const void* src, size_t n) {
    DeepCopyAsyncSYCL(instance, dst, src, n);
  }
};

template <class MemSpace1, class MemSpace2>
struct DeepCopy<MemSpace1, MemSpace2, Kokkos::Experimental::SYCL,
                std::enable_if_t<is_sycl_type_space<MemSpace1>::value &&
                                 is_sycl_type_space<MemSpace2>::value>> {
  DeepCopy(void* dst, const void* src, size_t n) { DeepCopySYCL(dst, src, n); }
  DeepCopy(const Kokkos::Experimental::SYCL& instance, void* dst,
           const void* src, size_t n) {
    DeepCopyAsyncSYCL(instance, dst, src, n);
  }
};

template <class MemSpace1, class MemSpace2, class ExecutionSpace>
struct DeepCopy<
    MemSpace1, MemSpace2, ExecutionSpace,
    std::enable_if_t<
        is_sycl_type_space<MemSpace1>::value &&
        is_sycl_type_space<MemSpace2>::value &&
        !std::is_same<ExecutionSpace, Kokkos::Experimental::SYCL>::value>> {
  inline DeepCopy(void* dst, const void* src, size_t n) {
    DeepCopySYCL(dst, src, n);
  }

  inline DeepCopy(const ExecutionSpace& exec, void* dst, const void* src,
                  size_t n) {
    exec.fence(fence_string());
    DeepCopyAsyncSYCL(dst, src, n);
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
struct DeepCopy<
    MemSpace, HostSpace, ExecutionSpace,
    std::enable_if_t<
        is_sycl_type_space<MemSpace>::value &&
        !std::is_same<ExecutionSpace, Kokkos::Experimental::SYCL>::value>> {
  inline DeepCopy(void* dst, const void* src, size_t n) {
    DeepCopySYCL(dst, src, n);
  }

  inline DeepCopy(const ExecutionSpace& exec, void* dst, const void* src,
                  size_t n) {
    exec.fence(fence_string());
    DeepCopyAsyncSYCL(dst, src, n);
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
struct DeepCopy<
    HostSpace, MemSpace, ExecutionSpace,
    std::enable_if_t<
        is_sycl_type_space<MemSpace>::value &&
        !std::is_same<ExecutionSpace, Kokkos::Experimental::SYCL>::value>> {
  inline DeepCopy(void* dst, const void* src, size_t n) {
    DeepCopySYCL(dst, src, n);
  }

  inline DeepCopy(const ExecutionSpace& exec, void* dst, const void* src,
                  size_t n) {
    exec.fence(fence_string());
    DeepCopyAsyncSYCL(dst, src, n);
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
#endif
