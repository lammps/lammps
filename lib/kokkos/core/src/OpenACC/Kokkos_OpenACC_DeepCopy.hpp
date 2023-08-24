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

#ifndef KOKKOS_OPENACC_DEEP_COPY_HPP
#define KOKKOS_OPENACC_DEEP_COPY_HPP

#include <OpenACC/Kokkos_OpenACC.hpp>
#include <OpenACC/Kokkos_OpenACCSpace.hpp>

#include <Kokkos_Concepts.hpp>

#include <openacc.h>

template <>
struct Kokkos::Impl::DeepCopy<Kokkos::Experimental::OpenACCSpace,
                              Kokkos::Experimental::OpenACCSpace,
                              Kokkos::Experimental::OpenACC> {
  DeepCopy(void* dst, const void* src, size_t n) {
    // The behavior of acc_memcpy_device when bytes argument is zero is
    // clarified only in the latest OpenACC specification (V3.2), and thus the
    // value checking is added as a safeguard. (The current NVHPC (V22.5)
    // supports OpenACC V2.7.)
    if (n > 0) {
      acc_memcpy_device(dst, const_cast<void*>(src), n);
    }
  }
  DeepCopy(const Kokkos::Experimental::OpenACC& exec, void* dst,
           const void* src, size_t n) {
    if (n > 0) {
      acc_memcpy_device_async(dst, const_cast<void*>(src), n,
                              exec.acc_async_queue());
    }
  }
};

template <class ExecutionSpace>
struct Kokkos::Impl::DeepCopy<Kokkos::Experimental::OpenACCSpace,
                              Kokkos::Experimental::OpenACCSpace,
                              ExecutionSpace> {
  DeepCopy(void* dst, const void* src, size_t n) {
    if (n > 0) {
      acc_memcpy_device(dst, const_cast<void*>(src), n);
    }
  }
  DeepCopy(const ExecutionSpace& exec, void* dst, const void* src, size_t n) {
    exec.fence(
        "Kokkos::Impl::DeepCopy<OpenACCSpace, OpenACCSpace, "
        "ExecutionSpace>::DeepCopy: fence before copy");
    if (n > 0) {
      acc_memcpy_device(dst, const_cast<void*>(src), n);
    }
  }
};

template <>
struct Kokkos::Impl::DeepCopy<Kokkos::Experimental::OpenACCSpace,
                              Kokkos::HostSpace,
                              Kokkos::Experimental::OpenACC> {
  DeepCopy(void* dst, const void* src, size_t n) {
    if (n > 0) acc_memcpy_to_device(dst, const_cast<void*>(src), n);
  }
  DeepCopy(const Kokkos::Experimental::OpenACC& exec, void* dst,
           const void* src, size_t n) {
    if (n > 0)
      acc_memcpy_to_device_async(dst, const_cast<void*>(src), n,
                                 exec.acc_async_queue());
  }
};

template <class ExecutionSpace>
struct Kokkos::Impl::DeepCopy<Kokkos::Experimental::OpenACCSpace,
                              Kokkos::HostSpace, ExecutionSpace> {
  DeepCopy(void* dst, const void* src, size_t n) {
    if (n > 0) {
      acc_memcpy_to_device(dst, const_cast<void*>(src), n);
    }
  }
  DeepCopy(const ExecutionSpace& exec, void* dst, const void* src, size_t n) {
    exec.fence(
        "Kokkos::Impl::DeepCopy<OpenACCSpace, HostSpace, "
        "ExecutionSpace>::DeepCopy: fence before copy");
    if (n > 0) {
      acc_memcpy_to_device(dst, const_cast<void*>(src), n);
    }
  }
};

template <>
struct Kokkos::Impl::DeepCopy<Kokkos::HostSpace,
                              Kokkos::Experimental::OpenACCSpace,
                              Kokkos::Experimental::OpenACC> {
  DeepCopy(void* dst, const void* src, size_t n) {
    if (n > 0) {
      acc_memcpy_from_device(dst, const_cast<void*>(src), n);
    }
  }
  DeepCopy(const Kokkos::Experimental::OpenACC& exec, void* dst,
           const void* src, size_t n) {
    if (n > 0) {
      acc_memcpy_from_device_async(dst, const_cast<void*>(src), n,
                                   exec.acc_async_queue());
    }
  }
};

template <class ExecutionSpace>
struct Kokkos::Impl::DeepCopy<
    Kokkos::HostSpace, Kokkos::Experimental::OpenACCSpace, ExecutionSpace> {
  DeepCopy(void* dst, const void* src, size_t n) {
    if (n > 0) acc_memcpy_from_device(dst, const_cast<void*>(src), n);
  }
  DeepCopy(const ExecutionSpace& exec, void* dst, const void* src, size_t n) {
    exec.fence(
        "Kokkos::Impl::DeepCopy<HostSpace, OpenACCSpace, "
        "ExecutionSpace>::DeepCopy: fence before copy");
    if (n > 0) {
      acc_memcpy_from_device(dst, const_cast<void*>(src), n);
    }
  }
};

#endif
