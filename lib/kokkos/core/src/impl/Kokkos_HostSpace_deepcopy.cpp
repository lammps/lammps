/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include "Kokkos_Core.hpp"
#include "Kokkos_HostSpace_deepcopy.hpp"

namespace Kokkos {

namespace Impl {

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

void hostspace_parallel_deepcopy_async(const DefaultHostExecutionSpace& exec,
                                       void* dst, const void* src,
                                       ptrdiff_t n) {
  using policy_t = Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>;
  constexpr int host_deep_copy_serial_limit = 10 * 8192;

  // If the asynchronous HPX backend is enabled, do *not* copy anything
  // synchronously. The deep copy must be correctly sequenced with respect to
  // other kernels submitted to the same instance, so we only use the fallback
  // parallel_for version in this case.
#if !(defined(KOKKOS_ENABLE_HPX) && defined(KOKKOS_ENABLE_HPX_ASYNC_DISPATCH))
  if ((n < host_deep_copy_serial_limit) ||
      (DefaultHostExecutionSpace().concurrency() == 1)) {
    std::memcpy(dst, src, n);
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

}  // namespace Impl

}  // namespace Kokkos
