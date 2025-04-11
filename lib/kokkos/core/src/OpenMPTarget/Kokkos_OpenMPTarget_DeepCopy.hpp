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

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#include <Kokkos_Macros.hpp>
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#endif
#ifndef KOKKOS_OPENMPTARGET_DEEP_COPY_HPP
#define KOKKOS_OPENMPTARGET_DEEP_COPY_HPP

#include <OpenMPTarget/Kokkos_OpenMPTarget_Error.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

// TODO: implement all possible deep_copies
template <class ExecutionSpace>
struct DeepCopy<Kokkos::Experimental::OpenMPTargetSpace,
                Kokkos::Experimental::OpenMPTargetSpace, ExecutionSpace> {
  DeepCopy(void* dst, const void* src, size_t n) {
    // In the Release and RelWithDebInfo builds, the size of the memcpy should
    // be greater than zero to avoid error. omp_target_memcpy returns zero on
    // success.
    if (n > 0)
      KOKKOS_IMPL_OMPT_SAFE_CALL(omp_target_memcpy(
          dst, const_cast<void*>(src), n, 0, 0, omp_get_default_device(),
          omp_get_default_device()));
  }
  DeepCopy(const ExecutionSpace& exec, void* dst, const void* src, size_t n) {
    exec.fence(
        "Kokkos::Impl::DeepCopy<OpenMPTargetSpace, OpenMPTargetSpace>: fence "
        "before "
        "copy");
    if (n > 0)
      KOKKOS_IMPL_OMPT_SAFE_CALL(omp_target_memcpy(
          dst, const_cast<void*>(src), n, 0, 0, omp_get_default_device(),
          omp_get_default_device()));
  }
};

template <class ExecutionSpace>
struct DeepCopy<Kokkos::Experimental::OpenMPTargetSpace, HostSpace,
                ExecutionSpace> {
  DeepCopy(void* dst, const void* src, size_t n) {
    if (n > 0)
      KOKKOS_IMPL_OMPT_SAFE_CALL(omp_target_memcpy(
          dst, const_cast<void*>(src), n, 0, 0, omp_get_default_device(),
          omp_get_initial_device()));
  }
  DeepCopy(const ExecutionSpace& exec, void* dst, const void* src, size_t n) {
    exec.fence(
        "Kokkos::Impl::DeepCopy<OpenMPTargetSpace, HostSpace>: fence before "
        "copy");
    if (n > 0)
      KOKKOS_IMPL_OMPT_SAFE_CALL(omp_target_memcpy(
          dst, const_cast<void*>(src), n, 0, 0, omp_get_default_device(),
          omp_get_initial_device()));
  }
};

template <class ExecutionSpace>
struct DeepCopy<HostSpace, Kokkos::Experimental::OpenMPTargetSpace,
                ExecutionSpace> {
  DeepCopy(void* dst, const void* src, size_t n) {
    if (n > 0)
      KOKKOS_IMPL_OMPT_SAFE_CALL(omp_target_memcpy(
          dst, const_cast<void*>(src), n, 0, 0, omp_get_initial_device(),
          omp_get_default_device()));
  }
  DeepCopy(const ExecutionSpace& exec, void* dst, const void* src, size_t n) {
    exec.fence(
        "Kokkos::Impl::DeepCopy<HostSpace, OpenMPTargetSpace>: fence before "
        "copy");
    if (n > 0)
      KOKKOS_IMPL_OMPT_SAFE_CALL(omp_target_memcpy(
          dst, const_cast<void*>(src), n, 0, 0, omp_get_initial_device(),
          omp_get_default_device()));
  }
};

}  // namespace Impl
}  // namespace Kokkos

#endif  // KOKKOS_OPENMPTARGET_DEEP_COPY_HPP
