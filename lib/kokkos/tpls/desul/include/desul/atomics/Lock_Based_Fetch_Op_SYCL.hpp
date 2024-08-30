/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_LOCK_BASED_FETCH_OP_SYCL_HPP_
#define DESUL_ATOMICS_LOCK_BASED_FETCH_OP_SYCL_HPP_

#include <desul/atomics/Common.hpp>
#include <desul/atomics/Lock_Array_SYCL.hpp>
#include <desul/atomics/Thread_Fence_SYCL.hpp>
#include <type_traits>

namespace desul {
namespace Impl {

template <class Oper,
          class T,
          class MemoryOrder,
          class MemoryScope,
          // equivalent to:
          //   requires !atomic_always_lock_free(sizeof(T))
          std::enable_if_t<!atomic_always_lock_free(sizeof(T)), int> = 0>
T device_atomic_fetch_oper(const Oper& op,
                           T* const dest,
                           dont_deduce_this_parameter_t<const T> val,
                           MemoryOrder /*order*/,
                           MemoryScope scope) {
  // This is a way to avoid deadlock in a subgroup
  T return_val;
  int done = 0;
  auto sg = sycl::ext::oneapi::experimental::this_sub_group();
  using sycl::ext::oneapi::group_ballot;
  using sycl::ext::oneapi::sub_group_mask;
  sub_group_mask active = group_ballot(sg, 1);
  sub_group_mask done_active = group_ballot(sg, 0);
  while (active != done_active) {
    if (!done) {
      if (lock_address_sycl((void*)dest, scope)) {
        device_atomic_thread_fence(MemoryOrderAcquire(), scope);
        return_val = *dest;
        *dest = op.apply(return_val, val);
        device_atomic_thread_fence(MemoryOrderRelease(), scope);
        unlock_address_sycl((void*)dest, scope);
        done = 1;
      }
    }
    done_active = group_ballot(sg, done);
  }
  return return_val;
}

template <class Oper,
          class T,
          class MemoryOrder,
          class MemoryScope,
          // equivalent to:
          //   requires !atomic_always_lock_free(sizeof(T))
          std::enable_if_t<!atomic_always_lock_free(sizeof(T)), int> = 0>
T device_atomic_oper_fetch(const Oper& op,
                           T* const dest,
                           dont_deduce_this_parameter_t<const T> val,
                           MemoryOrder /*order*/,
                           MemoryScope scope) {
  // This is a way to avoid deadlock in a subgroup
  T return_val;
  int done = 0;
  auto sg = sycl::ext::oneapi::experimental::this_sub_group();
  using sycl::ext::oneapi::group_ballot;
  using sycl::ext::oneapi::sub_group_mask;
  sub_group_mask active = group_ballot(sg, 1);
  sub_group_mask done_active = group_ballot(sg, 0);
  while (active != done_active) {
    if (!done) {
      if (lock_address_sycl((void*)dest, scope)) {
        device_atomic_thread_fence(MemoryOrderAcquire(), scope);
        return_val = op.apply(*dest, val);
        *dest = return_val;
        device_atomic_thread_fence(MemoryOrderRelease(), scope);
        unlock_address_sycl((void*)dest, scope);
        done = 1;
      }
    }
    done_active = group_ballot(sg, done);
  }
  return return_val;
}
}  // namespace Impl
}  // namespace desul

#endif
