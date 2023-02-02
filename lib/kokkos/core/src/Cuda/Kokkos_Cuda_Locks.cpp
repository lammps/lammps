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

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#endif

#include <Kokkos_Core.hpp>
#ifdef KOKKOS_ENABLE_CUDA
#include <Cuda/Kokkos_Cuda_Locks.hpp>
#include <Cuda/Kokkos_Cuda_Error.hpp>

#ifdef KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE
namespace Kokkos {
namespace Impl {
__device__ __constant__ CudaLockArrays g_device_cuda_lock_arrays = {nullptr, 0};
}
}  // namespace Kokkos
#endif

namespace Kokkos {

namespace {

__global__ void init_lock_array_kernel_atomic() {
  unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < CUDA_SPACE_ATOMIC_MASK + 1) {
    Kokkos::Impl::g_device_cuda_lock_arrays.atomic[i] = 0;
  }
}

}  // namespace

namespace Impl {

CudaLockArrays g_host_cuda_lock_arrays = {nullptr, 0};

void initialize_host_cuda_lock_arrays() {
#ifdef KOKKOS_ENABLE_IMPL_DESUL_ATOMICS
  desul::Impl::init_lock_arrays();

  DESUL_ENSURE_CUDA_LOCK_ARRAYS_ON_DEVICE();
#endif
  if (g_host_cuda_lock_arrays.atomic != nullptr) return;
  KOKKOS_IMPL_CUDA_SAFE_CALL(
      cudaMalloc(&g_host_cuda_lock_arrays.atomic,
                 sizeof(int) * (CUDA_SPACE_ATOMIC_MASK + 1)));
  Impl::cuda_device_synchronize(
      "Kokkos::Impl::initialize_host_cuda_lock_arrays: Pre Init Lock Arrays");
  g_host_cuda_lock_arrays.n = Cuda::concurrency();
  KOKKOS_COPY_CUDA_LOCK_ARRAYS_TO_DEVICE();
  init_lock_array_kernel_atomic<<<(CUDA_SPACE_ATOMIC_MASK + 1 + 255) / 256,
                                  256>>>();
  Impl::cuda_device_synchronize(
      "Kokkos::Impl::initialize_host_cuda_lock_arrays: Post Init Lock Arrays");
}

void finalize_host_cuda_lock_arrays() {
#ifdef KOKKOS_ENABLE_IMPL_DESUL_ATOMICS
  desul::Impl::finalize_lock_arrays();
#endif

  if (g_host_cuda_lock_arrays.atomic == nullptr) return;
  cudaFree(g_host_cuda_lock_arrays.atomic);
  g_host_cuda_lock_arrays.atomic = nullptr;
  g_host_cuda_lock_arrays.n      = 0;
#ifdef KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE
  KOKKOS_COPY_CUDA_LOCK_ARRAYS_TO_DEVICE();
#endif
}

}  // namespace Impl

}  // namespace Kokkos

#else

void KOKKOS_CORE_SRC_CUDA_CUDA_LOCKS_PREVENT_LINK_ERROR() {}

#endif
