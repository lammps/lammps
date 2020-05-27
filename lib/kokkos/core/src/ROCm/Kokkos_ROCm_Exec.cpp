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
#ifndef KOKKOS_ROCMEXEC_HPP
#define KOKKOS_ROCMEXEC_HPP

#include <algorithm>
#include <typeinfo>
#include <Kokkos_Macros.hpp>
//#include <ROCm/Kokkos_ROCmExec.hpp>
#include <hc.hpp>

#define ROCM_SPACE_ATOMIC_MASK 0x1FFFF
#define ROCM_SPACE_ATOMIC_XOR_MASK 0x15A39
#define ROCM_CONCURRENCY 20480
//#define ROCM_CONCURRENCY 81920  # for fiji

namespace Kokkos {
static int rocm_space_atomic_locks[ROCM_SPACE_ATOMIC_MASK + 1];
static int rocm_space_scratch_locks[ROCM_CONCURRENCY];
static int rocm_space_threadid_locks[ROCM_CONCURRENCY];
namespace Impl {
// TODO: mimic cuda implemtation, add dgpu capability

void init_rocm_atomic_lock_array() {
  static int is_initialized = 0;
  if (!is_initialized) {
    for (int i = 0; i < ROCM_SPACE_ATOMIC_MASK + 1; i++)
      rocm_space_atomic_locks[i] = 0;
    is_initialized = 1;
  }
}

void init_rocm_scratch_lock_array() {
  static int is_initialized = 0;
  if (!is_initialized) {
    for (int i = 0; i < ROCM_CONCURRENCY; i++) rocm_space_scratch_locks[i] = 0;
    is_initialized = 1;
  }
}

void init_rocm_threadid_lock_array() {
  static int is_initialized = 0;
  if (!is_initialized) {
    for (int i = 0; i < ROCM_CONCURRENCY; i++) rocm_space_threadid_locks[i] = 0;
    is_initialized = 1;
  }
}

void init_lock_arrays_rocm_space() {
  init_rocm_atomic_lock_array();
  //     init_rocm_scratch_lock_array();
  //     init_rocm_threadid_lock_array();
}
}  // namespace Impl

}  // namespace Kokkos
#if 0
namespace Kokkos {
namespace Impl {
KOKKOS_INLINE_FUNCTION
bool lock_address_rocm_space(void* ptr) {
#if 0
return(Kokkos::Impl::lock_address_host_space(ptr));
#else
  size_t offset = size_t(ptr);
  offset = offset >> 2;
  offset = offset & ROCM_SPACE_ATOMIC_MASK;
  return (0 == hc::atomic_compare_exchange(&rocm_space_atomic_locks[offset],0,1));
#endif
}

KOKKOS_INLINE_FUNCTION
void unlock_address_rocm_space(void* ptr) {
#if 0
Kokkos::Impl::unlock_address_host_space(ptr) ;
#else
  size_t offset = size_t(ptr);
  offset = offset >> 2;
  offset = offset & ROCM_SPACE_ATOMIC_MASK;
  hc::atomic_exchange( &rocm_space_atomic_locks[ offset ], 0);
#endif
}

}
} // namespace Kokkos
#endif

#endif /* #ifndef KOKKOS_ROCMEXEC_HPP */
