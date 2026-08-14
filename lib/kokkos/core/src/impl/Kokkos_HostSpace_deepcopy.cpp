// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#endif

#include "Kokkos_Core.hpp"
#include "Kokkos_HostSpace_deepcopy.hpp"

namespace Kokkos {

namespace Impl {

void hostspace_fence(const DefaultHostExecutionSpace& exec) {
  exec.fence("HostSpace fence");
}

void hostspace_parallel_deepcopy(void* dst, const void* src, ptrdiff_t n) {
  Kokkos::DefaultHostExecutionSpace exec;
  hostspace_parallel_deepcopy_async(exec, dst, src, n);
}

// DeepCopy called with an execution space that can't access HostSpace
void hostspace_parallel_deepcopy_async(void* dst, const void* src,
                                       ptrdiff_t n) {
  Kokkos::DefaultHostExecutionSpace exec;
  hostspace_parallel_deepcopy_async(exec, dst, src, n);
  exec.fence(
      "Kokkos::Impl::hostspace_parallel_deepcopy_async: fence after copy");
}

template <typename ExecutionSpace>
void hostspace_parallel_deepcopy_async(const ExecutionSpace& exec, void* dst,
                                       const void* src, ptrdiff_t n) {
  using policy_t = Kokkos::RangePolicy<ExecutionSpace>;

  // If the asynchronous HPX backend is enabled, do *not* copy anything
  // synchronously. The deep copy must be correctly sequenced with respect to
  // other kernels submitted to the same instance, so we only use the fallback
  // parallel_for version in this case.
#if !(defined(KOKKOS_ENABLE_HPX) && \
      defined(KOKKOS_ENABLE_IMPL_HPX_ASYNC_DISPATCH))
  constexpr int host_deep_copy_serial_limit = 10 * 8192;
  if ((n < host_deep_copy_serial_limit) || (exec.concurrency() == 1)) {
    if (0 < n) std::memcpy(dst, src, n);
    return;
  }

  // Both src and dst are aligned the same way with respect to 8 byte words
  if (reinterpret_cast<ptrdiff_t>(src) % 8 ==
      reinterpret_cast<ptrdiff_t>(dst) % 8) {
    char* dst_c       = reinterpret_cast<char*>(dst);
    const char* src_c = reinterpret_cast<const char*>(src);
    int count         = 0;
    // get initial bytes copied
    while (reinterpret_cast<ptrdiff_t>(dst_c) % 8 != 0) {
      *dst_c = *src_c;
      dst_c++;
      src_c++;
      count++;
    }

    // copy the bulk of the data
    double* dst_p       = reinterpret_cast<double*>(dst_c);
    const double* src_p = reinterpret_cast<const double*>(src_c);
    Kokkos::parallel_for("Kokkos::Impl::host_space_deepcopy_double",
                         policy_t(exec, 0, (n - count) / 8),
                         [=](const ptrdiff_t i) { dst_p[i] = src_p[i]; });

    // get final data copied
    dst_c += ((n - count) / 8) * 8;
    src_c += ((n - count) / 8) * 8;
    char* dst_end = reinterpret_cast<char*>(dst) + n;
    while (dst_c != dst_end) {
      *dst_c = *src_c;
      dst_c++;
      src_c++;
    }
    return;
  }

  // Both src and dst are aligned the same way with respect to 4 byte words
  if (reinterpret_cast<ptrdiff_t>(src) % 4 ==
      reinterpret_cast<ptrdiff_t>(dst) % 4) {
    char* dst_c       = reinterpret_cast<char*>(dst);
    const char* src_c = reinterpret_cast<const char*>(src);
    int count         = 0;
    // get initial bytes copied
    while (reinterpret_cast<ptrdiff_t>(dst_c) % 4 != 0) {
      *dst_c = *src_c;
      dst_c++;
      src_c++;
      count++;
    }

    // copy the bulk of the data
    int32_t* dst_p       = reinterpret_cast<int32_t*>(dst_c);
    const int32_t* src_p = reinterpret_cast<const int32_t*>(src_c);
    Kokkos::parallel_for("Kokkos::Impl::host_space_deepcopy_int",
                         policy_t(exec, 0, (n - count) / 4),
                         [=](const ptrdiff_t i) { dst_p[i] = src_p[i]; });

    // get final data copied
    dst_c += ((n - count) / 4) * 4;
    src_c += ((n - count) / 4) * 4;
    char* dst_end = reinterpret_cast<char*>(dst) + n;
    while (dst_c != dst_end) {
      *dst_c = *src_c;
      dst_c++;
      src_c++;
    }
    return;
  }
#endif

  // Src and dst are not aligned the same way, we can only to byte wise copy.
  {
    char* dst_p       = reinterpret_cast<char*>(dst);
    const char* src_p = reinterpret_cast<const char*>(src);
    Kokkos::parallel_for("Kokkos::Impl::host_space_deepcopy_char",
                         policy_t(exec, 0, n),
                         [=](const ptrdiff_t i) { dst_p[i] = src_p[i]; });
  }
}

template <typename ExecutionSpace>
void hostspace_parallel_zeromemset(const ExecutionSpace& exec, void* dst,
                                   size_t n) {
  constexpr uint8_t z_u8 = 0x00;
  using policy_t =
      Kokkos::RangePolicy<ExecutionSpace, Kokkos::IndexType<size_t>>;

#if !(defined(KOKKOS_ENABLE_HPX) && \
      defined(KOKKOS_ENABLE_IMPL_HPX_ASYNC_DISPATCH))
  constexpr uint64_t z_u64 = 0x0000000000000000;

  // Align initial bytes to 8-byte boundary
  uint8_t* dst_c = reinterpret_cast<uint8_t*>(dst);
  size_t count   = 0;
  while (reinterpret_cast<size_t>(dst_c) % 8 != 0) {
    *dst_c = z_u8;
    dst_c++;
    count++;
  }

  // Zero-fill with 8-byte words in parallel
  uint64_t* dst_p  = reinterpret_cast<uint64_t*>(dst_c);
  const size_t cnt = (n - count) / 8;
  Kokkos::parallel_for("Kokkos::Impl::hostspace_parallel_zeromemset",
                       policy_t(exec, 0, cnt),
                       [=](const size_t i) { dst_p[i] = z_u64; });

  // Handle any remaining bytes that don't fit in 8-byte words
  dst_c += cnt * 8;
  uint8_t* dst_end = reinterpret_cast<uint8_t*>(dst) + n;
  while (dst_c != dst_end) {
    *dst_c = z_u8;
    dst_c++;
  }
#else
  uint8_t* dst_p = reinterpret_cast<uint8_t*>(dst);
  Kokkos::parallel_for("Kokkos::Impl::hostspace_parallel_zeromemset",
                       policy_t(exec, 0, n),
                       [=](const size_t i) { dst_p[i] = z_u8; });
#endif
}

// Explicit instantiation
template void hostspace_parallel_deepcopy_async<DefaultHostExecutionSpace>(
    const DefaultHostExecutionSpace&, void*, const void*, ptrdiff_t);

template void hostspace_parallel_zeromemset<DefaultHostExecutionSpace>(
    const DefaultHostExecutionSpace&, void*, size_t);

#if defined(KOKKOS_ENABLE_SERIAL) &&                                    \
    (defined(KOKKOS_ENABLE_OPENMP) || defined(KOKKOS_ENABLE_THREADS) || \
     defined(KOKKOS_ENABLE_HPX))
// Instantiate only if both the Serial backend and some other host parallel
// backend are enabled
template void hostspace_parallel_deepcopy_async<Kokkos::Serial>(
    const Kokkos::Serial&, void*, const void*, ptrdiff_t);
#endif
}  // namespace Impl

}  // namespace Kokkos
