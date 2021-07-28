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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ATOMIC_HPP) && !defined(KOKKOS_MEMORY_FENCE_HPP)
#define KOKKOS_MEMORY_FENCE_HPP
namespace Kokkos {

//----------------------------------------------------------------------------

KOKKOS_FORCEINLINE_FUNCTION
void memory_fence() {
#if defined(__CUDA_ARCH__)
  __threadfence();
#elif defined(KOKKOS_ENABLE_OPENMPTARGET)
#pragma omp flush
#elif defined(__HIP_DEVICE_COMPILE__)
  __threadfence();
#elif defined(KOKKOS_ENABLE_SYCL) && defined(__SYCL_DEVICE_ONLY__)
  sycl::ONEAPI::atomic_fence(sycl::ONEAPI::memory_order::acq_rel,
                             sycl::ONEAPI::memory_scope::device);
#elif defined(KOKKOS_ENABLE_ASM) && defined(KOKKOS_ENABLE_ISA_X86_64)
  asm volatile("mfence" ::: "memory");
#elif defined(KOKKOS_ENABLE_GNU_ATOMICS) || \
    (defined(KOKKOS_COMPILER_NVCC) && defined(KOKKOS_ENABLE_INTEL_ATOMICS))
  __sync_synchronize();
#elif defined(KOKKOS_ENABLE_INTEL_ATOMICS)
  _mm_mfence();
#elif defined(KOKKOS_ENABLE_OPENMP_ATOMICS)
#pragma omp flush
#elif defined(KOKKOS_ENABLE_WINDOWS_ATOMICS)
  MemoryBarrier();
#elif !defined(KOKKOS_ENABLE_SERIAL_ATOMICS)
#error "Error: memory_fence() not defined"
#endif
}

//////////////////////////////////////////////////////
// store_fence()
//
// If possible use a store fence on the architecture, if not run a full memory
// fence

KOKKOS_FORCEINLINE_FUNCTION
void store_fence() {
#if defined(KOKKOS_ENABLE_ASM) && defined(KOKKOS_ENABLE_ISA_X86_64)
  asm volatile("sfence" ::: "memory");
#else
  memory_fence();
#endif
}

//////////////////////////////////////////////////////
// load_fence()
//
// If possible use a load fence on the architecture, if not run a full memory
// fence

KOKKOS_FORCEINLINE_FUNCTION
void load_fence() {
#if defined(KOKKOS_ENABLE_ASM) && defined(KOKKOS_ENABLE_ISA_X86_64)
  asm volatile("lfence" ::: "memory");
#else
  memory_fence();
#endif
}

}  // namespace Kokkos

#endif
