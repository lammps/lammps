// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_IMPL_HOSTSPACE_DEEPCOPY_HPP
#define KOKKOS_IMPL_HOSTSPACE_DEEPCOPY_HPP

#include <cstdint>

namespace Kokkos {

namespace Impl {

void hostspace_fence(const DefaultHostExecutionSpace& exec);

void hostspace_parallel_deepcopy(void* dst, const void* src, ptrdiff_t n);
// DeepCopy called with an execution space that can't access HostSpace
void hostspace_parallel_deepcopy_async(void* dst, const void* src, ptrdiff_t n);
template <typename ExecutionSpace>
void hostspace_parallel_deepcopy_async(const ExecutionSpace& exec, void* dst,
                                       const void* src, ptrdiff_t n);
template <typename ExecutionSpace>
void hostspace_parallel_zeromemset(const ExecutionSpace& exec, void* dst,
                                   size_t n);
}  // namespace Impl

}  // namespace Kokkos

#endif  // KOKKOS_IMPL_HOSTSPACE_DEEPCOPY_HPP
