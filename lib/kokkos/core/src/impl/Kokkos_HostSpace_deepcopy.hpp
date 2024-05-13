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

#ifndef KOKKOS_IMPL_HOSTSPACE_DEEPCOPY_HPP
#define KOKKOS_IMPL_HOSTSPACE_DEEPCOPY_HPP

#include <cstdint>

namespace Kokkos {

namespace Impl {

void hostspace_fence(const DefaultHostExecutionSpace& exec);

void hostspace_parallel_deepcopy(void* dst, const void* src, ptrdiff_t n);
// DeepCopy called with an execution space that can't access HostSpace
void hostspace_parallel_deepcopy_async(void* dst, const void* src, ptrdiff_t n);
void hostspace_parallel_deepcopy_async(const DefaultHostExecutionSpace& exec,
                                       void* dst, const void* src, ptrdiff_t n);

}  // namespace Impl

}  // namespace Kokkos

#endif  // KOKKOS_IMPL_HOSTSPACE_DEEPCOPY_HPP
