// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_OPENMP_ZEROMEMSET_HPP
#define KOKKOS_OPENMP_ZEROMEMSET_HPP

#include <OpenMP/Kokkos_OpenMP.hpp>
#include <impl/Kokkos_ZeroMemset_fwd.hpp>

#include <cstring>

namespace Kokkos {
namespace Impl {

template <>
struct ZeroMemset<OpenMP> {
  ZeroMemset(const OpenMP& exec_space, void* dst, size_t cnt) {
    // Threshold chosen based on the ViewFirstTouch_ParallelFor benchmark,
    // run on AMD EPYC Genoa and Intel Xeon Cascade Lake architectures,
    // which have 8 and 2 NUMA nodes respectively.
    constexpr size_t host_memset_limit = 1lu << 17;
    if (cnt < host_memset_limit || exec_space.concurrency() < 4) {
      std::memset(dst, 0, cnt);
    } else {
      hostspace_parallel_zeromemset(exec_space, dst, cnt);
    }
  }
};

}  // namespace Impl
}  // namespace Kokkos

#endif  // KOKKOS_OPENMP_ZEROMEMSET_HPP
