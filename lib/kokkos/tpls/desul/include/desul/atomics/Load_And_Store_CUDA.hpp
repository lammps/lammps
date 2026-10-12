/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_LOAD_AND_STORE_CUDA_HPP_
#define DESUL_ATOMICS_LOAD_AND_STORE_CUDA_HPP_

#include <desul/atomics/Compare_Exchange_CUDA.hpp>
#include <desul/atomics/Lock_Free_Types_CUDA.hpp>
#include <desul/atomics/Macros.hpp>

// Including CUDA ptx based exchange atomics
// When building with clang we need to include the device functions always
// since clang must see a consistent overload set in both device and host compilation
// but that means we need to know on the host what to make visible, i.e. we need
// a host side compile knowledge of architecture.
// We simply can say DESUL proper doesn't support clang CUDA build pre Volta,
// Kokkos has that knowledge and so I use it here, allowing in Kokkos to use
// clang with pre Volta as CUDA compiler
#ifndef DESUL_CUDA_ARCH_IS_PRE_VOLTA

#include <desul/atomics/cuda/CUDA_asm_loadstore.hpp>

// use lock-tables for the sizes that do not have an intrinsic.
namespace desul {
namespace Impl {
template <class T,
          class MemoryOrder,
          class MemoryScope,
          std::enable_if_t<!device_atomic_always_lock_free<T>, int> = 0>
__device__ T device_atomic_load(T const* ptr,
                                MemoryOrder /*order*/,
                                MemoryScope scope) {
  // This is a way to avoid deadlock in a warp or wave front
  T ret;
  int done = 0;
  unsigned int mask = __activemask();
  unsigned int active = __ballot_sync(mask, 1);
  unsigned int done_active = 0;
  while (active != done_active) {
    if (!done) {
      if (lock_address_cuda((void*)ptr, scope)) {
        device_atomic_thread_fence(MemoryOrderAcquire(), scope);
        ret = *ptr;
        device_atomic_thread_fence(MemoryOrderRelease(), scope);
        unlock_address_cuda((void*)ptr, scope);
        done = 1;
      }
    }
    done_active = __ballot_sync(mask, done);
  }
  return ret;
}

template <class T,
          class MemoryOrder,
          class MemoryScope,
          std::enable_if_t<!device_atomic_always_lock_free<T>, int> = 0>
__device__ void device_atomic_store(T* ptr,
                                    T val,
                                    MemoryOrder /*order*/,
                                    MemoryScope scope) {
  // This is a way to avoid deadlock in a warp or wave front
  int done = 0;
  unsigned int mask = __activemask();
  unsigned int active = __ballot_sync(mask, 1);
  unsigned int done_active = 0;
  while (active != done_active) {
    if (!done) {
      if (lock_address_cuda((void*)ptr, scope)) {
        device_atomic_thread_fence(MemoryOrderAcquire(), scope);
        *ptr = val;
        device_atomic_thread_fence(MemoryOrderRelease(), scope);
        unlock_address_cuda((void*)ptr, scope);
        done = 1;
      }
    }
    done_active = __ballot_sync(mask, done);
  }
}

}  // namespace Impl
}  // namespace desul

#else
// pre-Volta there was no load and store instructions, so we have to use CAS for legacy
// support. This is violating C++ semantics, as we potentially do a RMW op on a const.
#include <desul/atomics/Compare_Exchange_CUDA.hpp>
namespace desul {
namespace Impl {

DESUL_IMPL_ATOMIC_LOAD_AND_STORE_WITH_CAS(DESUL_IMPL_DEVICE_FUNCTION, device)

}
}  // namespace desul
#endif
#endif
