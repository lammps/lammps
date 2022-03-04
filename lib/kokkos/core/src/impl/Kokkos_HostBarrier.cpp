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

#include <Kokkos_Macros.hpp>

#include <impl/Kokkos_HostBarrier.hpp>
#include <impl/Kokkos_BitOps.hpp>

#include <impl/Kokkos_HostBarrier.hpp>

#include <thread>
#if defined(_WIN32)
#include <process.h>
#include <winsock2.h>
#include <windows.h>
#endif

namespace Kokkos {
namespace Impl {

void HostBarrier::impl_backoff_wait_until_equal(
    int* ptr, const int v, const bool active_wait) noexcept {
  unsigned count = 0u;

  while (!test_equal(ptr, v)) {
    const int c = ::Kokkos::log2(++count);
    if (!active_wait || c > log2_iterations_till_sleep) {
      std::this_thread::sleep_for(
          std::chrono::nanoseconds(c < 16 ? 256 * c : 4096));
    } else if (c > log2_iterations_till_yield) {
      std::this_thread::yield();
    }
#if defined(KOKKOS_ENABLE_ASM)
#if defined(__PPC64__)
    for (int j = 0; j < num_nops; ++j) {
      asm volatile("nop\n");
    }
    asm volatile("or 27, 27, 27" ::: "memory");
#elif defined(__amd64) || defined(__amd64__) || defined(__x86_64) || \
    defined(__x86_64__)
    for (int j = 0; j < num_nops; ++j) {
      asm volatile("nop\n");
    }
    asm volatile("pause\n" ::: "memory");
#endif
#endif
  }
}
}  // namespace Impl
}  // namespace Kokkos
