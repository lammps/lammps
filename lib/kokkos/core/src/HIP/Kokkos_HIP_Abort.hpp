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

#ifndef KOKKOS_HIP_ABORT_HPP
#define KOKKOS_HIP_ABORT_HPP

#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ENABLE_HIP)

#include <hip/hip_runtime.h>

// FIXME_HIP ROCm 4.5 version header include would be <rocm/rocm_version.h>
#if __has_include(<rocm_version.h>)
#include <rocm_version.h>
#define KOKKOS_IMPL_ROCM_VERSION \
  ROCM_VERSION_MAJOR * 10000 + ROCM_VERSION_MINOR * 100 + ROCM_VERSION_PATCH
#endif

// FIXME_HIP workaround for ROCm version less than 5.0.2
#if KOKKOS_IMPL_ROCM_VERSION < 50002
#define KOKKOS_IMPL_HIP_ABORT_DOES_NOT_PRINT_MESSAGE
#endif

namespace Kokkos {
namespace Impl {

// The two keywords below are not contradictory. `noinline` is a
// directive to the optimizer.
[[noreturn]] __device__ __attribute__((noinline)) inline void hip_abort(
    char const *msg) {
  const char empty[] = "";
  __assert_fail(msg, empty, 0, empty);
  // This loop is never executed. It's intended to suppress warnings that the
  // function returns, even though it does not. This is necessary because
  // abort() is not marked as [[noreturn]], even though it does not return.
  while (true)
    ;
}

}  // namespace Impl
}  // namespace Kokkos

#endif
#endif
