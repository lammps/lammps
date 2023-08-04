/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

// FIXME_SYCL Use SYCL_EXT_ONEAPI_DEVICE_GLOBAL when available instead
#ifdef DESUL_SYCL_DEVICE_GLOBAL_SUPPORTED

#include <cinttypes>
#include <desul/atomics/Lock_Array_SYCL.hpp>

namespace desul::Impl {

SYCL_EXTERNAL
sycl_device_global<int32_t*> SYCL_SPACE_ATOMIC_LOCKS_DEVICE;
SYCL_EXTERNAL
sycl_device_global<int32_t*> SYCL_SPACE_ATOMIC_LOCKS_NODE;

int32_t* SYCL_SPACE_ATOMIC_LOCKS_DEVICE_h = nullptr;
int32_t* SYCL_SPACE_ATOMIC_LOCKS_NODE_h = nullptr;

template <>
void init_lock_arrays_sycl<int>(sycl::queue q) {
  if (SYCL_SPACE_ATOMIC_LOCKS_DEVICE_h != nullptr) return;

  SYCL_SPACE_ATOMIC_LOCKS_DEVICE_h =
      sycl::malloc_device<int32_t>(SYCL_SPACE_ATOMIC_MASK + 1, q);
  SYCL_SPACE_ATOMIC_LOCKS_NODE_h =
      sycl::malloc_host<int32_t>(SYCL_SPACE_ATOMIC_MASK + 1, q);

  // FIXME_SYCL Once supported, the following should be replaced by
  // q.memcpy(SYCL_SPACE_ATOMIC_LOCKS_DEVICE,
  //          &SYCL_SPACE_ATOMIC_LOCKS_DEVICE_h,
  //          sizeof(int32_t*));
  // q.memcpy(SYCL_SPACE_ATOMIC_LOCKS_NODE,
  //          &SYCL_SPACE_ATOMIC_LOCKS_NODE_h,
  //          sizeof(int32_t*));
  auto device_ptr = SYCL_SPACE_ATOMIC_LOCKS_DEVICE_h;
  auto node_ptr = SYCL_SPACE_ATOMIC_LOCKS_NODE_h;
  q.single_task([=] {
    SYCL_SPACE_ATOMIC_LOCKS_DEVICE.get() = device_ptr;
    SYCL_SPACE_ATOMIC_LOCKS_NODE.get() = node_ptr;
  });

  q.memset(SYCL_SPACE_ATOMIC_LOCKS_DEVICE_h,
           0,
           sizeof(int32_t) * (SYCL_SPACE_ATOMIC_MASK + 1));
  q.memset(SYCL_SPACE_ATOMIC_LOCKS_NODE_h,
           0,
           sizeof(int32_t) * (SYCL_SPACE_ATOMIC_MASK + 1));

  q.wait_and_throw();
}

template <>
void finalize_lock_arrays_sycl<int>(sycl::queue q) {
  if (SYCL_SPACE_ATOMIC_LOCKS_DEVICE_h == nullptr) return;

  sycl::free(SYCL_SPACE_ATOMIC_LOCKS_DEVICE_h, q);
  sycl::free(SYCL_SPACE_ATOMIC_LOCKS_NODE_h, q);
  SYCL_SPACE_ATOMIC_LOCKS_DEVICE_h = nullptr;
  SYCL_SPACE_ATOMIC_LOCKS_NODE_h = nullptr;
}

} // namespace desul::Impl
#endif
