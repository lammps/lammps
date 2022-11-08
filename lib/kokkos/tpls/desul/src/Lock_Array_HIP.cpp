/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#include <cinttypes>
#include <desul/atomics/Lock_Array.hpp>
#include <sstream>
#include <string>

#ifdef DESUL_HAVE_HIP_ATOMICS
#ifdef DESUL_HIP_RDC
namespace desul {
namespace Impl {
__device__ __constant__ int32_t* HIP_SPACE_ATOMIC_LOCKS_DEVICE = nullptr;
__device__ __constant__ int32_t* HIP_SPACE_ATOMIC_LOCKS_NODE = nullptr;
}  // namespace Impl
}  // namespace desul
#endif

namespace desul {

namespace {

__global__ void init_lock_arrays_hip_kernel() {
  unsigned i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < HIP_SPACE_ATOMIC_MASK + 1) {
    Impl::HIP_SPACE_ATOMIC_LOCKS_DEVICE[i] = 0;
    Impl::HIP_SPACE_ATOMIC_LOCKS_NODE[i] = 0;
  }
}

}  // namespace

namespace Impl {

int32_t* HIP_SPACE_ATOMIC_LOCKS_DEVICE_h = nullptr;
int32_t* HIP_SPACE_ATOMIC_LOCKS_NODE_h = nullptr;

// Putting this into anonymous namespace so we don't have multiple defined symbols
// When linking in more than one copy of the object file
namespace {

void check_error_and_throw_hip(hipError_t e, const std::string msg) {
  if (e != hipSuccess) {
    std::ostringstream out;
    out << "Desul::Error: " << msg << " error(" << hipGetErrorName(e)
        << "): " << hipGetErrorString(e);
    throw std::runtime_error(out.str());
  }
}

}  // namespace

template <typename T>
void init_lock_arrays_hip() {
  if (HIP_SPACE_ATOMIC_LOCKS_DEVICE_h != nullptr) return;

  auto error_malloc1 = hipMalloc(&HIP_SPACE_ATOMIC_LOCKS_DEVICE_h,
                                 sizeof(int32_t) * (HIP_SPACE_ATOMIC_MASK + 1));
  check_error_and_throw_hip(error_malloc1,
                            "init_lock_arrays_hip: hipMalloc device locks");

  auto error_malloc2 = hipHostMalloc(&HIP_SPACE_ATOMIC_LOCKS_NODE_h,
                                     sizeof(int32_t) * (HIP_SPACE_ATOMIC_MASK + 1));
  check_error_and_throw_hip(error_malloc2,
                            "init_lock_arrays_hip: hipMallocHost host locks");

  auto error_sync1 = hipDeviceSynchronize();
  DESUL_IMPL_COPY_HIP_LOCK_ARRAYS_TO_DEVICE();
  check_error_and_throw_hip(error_sync1, "init_lock_arrays_hip: post malloc");

  init_lock_arrays_hip_kernel<<<(HIP_SPACE_ATOMIC_MASK + 1 + 255) / 256, 256>>>();

  auto error_sync2 = hipDeviceSynchronize();
  check_error_and_throw_hip(error_sync2, "init_lock_arrays_hip: post init");
}

template <typename T>
void finalize_lock_arrays_hip() {
  if (HIP_SPACE_ATOMIC_LOCKS_DEVICE_h == nullptr) return;
  auto error_free1 = hipFree(HIP_SPACE_ATOMIC_LOCKS_DEVICE_h);
  check_error_and_throw_hip(error_free1, "finalize_lock_arrays_hip: free device locks");
  auto error_free2 = hipHostFree(HIP_SPACE_ATOMIC_LOCKS_NODE_h);
  check_error_and_throw_hip(error_free2, "finalize_lock_arrays_hip: free host locks");
  HIP_SPACE_ATOMIC_LOCKS_DEVICE_h = nullptr;
  HIP_SPACE_ATOMIC_LOCKS_NODE_h = nullptr;
#ifdef DESUL_HIP_RDC
  DESUL_IMPL_COPY_HIP_LOCK_ARRAYS_TO_DEVICE();
#endif
}

template void init_lock_arrays_hip<int>();
template void finalize_lock_arrays_hip<int>();

}  // namespace Impl

}  // namespace desul
#endif
