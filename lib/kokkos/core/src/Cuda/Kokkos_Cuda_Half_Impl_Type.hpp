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

#ifndef KOKKOS_CUDA_HALF_IMPL_TYPE_HPP_
#define KOKKOS_CUDA_HALF_IMPL_TYPE_HPP_

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_CUDA
#if !(defined(KOKKOS_COMPILER_CLANG) && KOKKOS_COMPILER_CLANG < 900) && \
    !(defined(KOKKOS_ARCH_KEPLER) || defined(KOKKOS_ARCH_MAXWELL50) ||  \
      defined(KOKKOS_ARCH_MAXWELL52))
#include <cuda_fp16.h>
#if (CUDA_VERSION >= 11000)
#include <cuda_bf16.h>
#endif  // CUDA_VERSION >= 11000

#ifndef KOKKOS_IMPL_HALF_TYPE_DEFINED
// Make sure no one else tries to define half_t
#define KOKKOS_IMPL_HALF_TYPE_DEFINED
#define KOKKOS_IMPL_CUDA_HALF_TYPE_DEFINED

namespace Kokkos {
namespace Impl {
struct half_impl_t {
  using type = __half;
};
#if (CUDA_VERSION >= 11000)
#define KOKKOS_IMPL_BHALF_TYPE_DEFINED
struct bhalf_impl_t {
  using type = __nv_bfloat16;
};
#endif  // CUDA_VERSION >= 11000
}  // namespace Impl
}  // namespace Kokkos
#endif  // KOKKOS_IMPL_HALF_TYPE_DEFINED
#endif  // Disables for half_t on cuda:
        // Clang/8||KEPLER30||KEPLER32||KEPLER37||MAXWELL50||MAXWELL52
#endif  // KOKKOS_ENABLE_CUDA
#endif
