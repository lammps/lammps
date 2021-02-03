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

#include <HIP/Kokkos_HIP_Locks.hpp>
#include <HIP/Kokkos_HIP_Error.hpp>
#include <Kokkos_HIP_Space.hpp>

#include <hip/hip_runtime.h>

#include <iostream>

namespace Kokkos {

#ifdef KOKKOS_ENABLE_HIP_RELOCATABLE_DEVICE_CODE
namespace Impl {
__device__ __constant__ HIPLockArrays g_device_hip_lock_arrays = {nullptr,
                                                                  nullptr, 0};
}
#endif

namespace {

__global__ void init_lock_array_kernel_atomic() {
  unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < KOKKOS_IMPL_HIP_SPACE_ATOMIC_MASK + 1) {
    Kokkos::Impl::g_device_hip_lock_arrays.atomic[i] = 0;
  }
}

__global__ void init_lock_array_kernel_threadid(int N) {
  unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < static_cast<unsigned>(N)) {
    Kokkos::Impl::g_device_hip_lock_arrays.scratch[i] = 0;
  }
}

}  // namespace

namespace Impl {

HIPLockArrays g_host_hip_lock_arrays = {nullptr, nullptr, 0};

void initialize_host_hip_lock_arrays() {
  if (g_host_hip_lock_arrays.atomic != nullptr) return;
  HIP_SAFE_CALL(hipMalloc(
      &g_host_hip_lock_arrays.atomic,
      sizeof(std::int32_t) * (KOKKOS_IMPL_HIP_SPACE_ATOMIC_MASK + 1)));
  HIP_SAFE_CALL(hipMalloc(
      &g_host_hip_lock_arrays.scratch,
      sizeof(std::int32_t) * (::Kokkos::Experimental::HIP::concurrency())));

  g_host_hip_lock_arrays.n = ::Kokkos::Experimental::HIP::concurrency();

  KOKKOS_COPY_HIP_LOCK_ARRAYS_TO_DEVICE();
  init_lock_array_kernel_atomic<<<
      (KOKKOS_IMPL_HIP_SPACE_ATOMIC_MASK + 1 + 255) / 256, 256, 0, nullptr>>>();
  init_lock_array_kernel_threadid<<<
      (::Kokkos::Experimental::HIP::concurrency() + 255) / 256, 256, 0,
      nullptr>>>(::Kokkos::Experimental::HIP::concurrency());
}

void finalize_host_hip_lock_arrays() {
  if (g_host_hip_lock_arrays.atomic == nullptr) return;
  HIP_SAFE_CALL(hipFree(g_host_hip_lock_arrays.atomic));
  g_host_hip_lock_arrays.atomic = nullptr;
  HIP_SAFE_CALL(hipFree(g_host_hip_lock_arrays.scratch));
  g_host_hip_lock_arrays.scratch = nullptr;
  g_host_hip_lock_arrays.n       = 0;
#ifdef KOKKOS_ENABLE_HIP_RELOCATABLE_DEVICE_CODE
  KOKKOS_COPY_HIP_LOCK_ARRAYS_TO_DEVICE();
#endif
}

}  // namespace Impl

}  // namespace Kokkos
