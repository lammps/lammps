/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_LOCK_BASED_FETCH_OP_HOST_HPP_
#define DESUL_ATOMICS_LOCK_BASED_FETCH_OP_HOST_HPP_

#include <desul/atomics/Common.hpp>
#ifdef DESUL_HAVE_GCC_ATOMICS
#include <desul/atomics/Compare_Exchange_GCC.hpp>
#endif
#ifdef DESUL_HAVE_MSVC_ATOMICS
#include <desul/atomics/Compare_Exchange_MSVC.hpp>
#endif
#ifdef DESUL_HAVE_OPENMP_ATOMICS
#include <desul/atomics/Compare_Exchange_OpenMP.hpp>
#endif
#include <desul/atomics/Lock_Array.hpp>
#include <desul/atomics/Thread_Fence.hpp>
#include <type_traits>

namespace desul {
namespace Impl {

template <class Oper,
          class T,
          class MemoryOrder,
          class MemoryScope,
          std::enable_if_t<!host_atomic_always_lock_free<T>, int> = 0>
inline T host_atomic_fetch_oper(const Oper& op,
                                T* const dest,
                                dont_deduce_this_parameter_t<const T> val,
                                MemoryOrder /*order*/,
                                MemoryScope scope) {
  // Acquire a lock for the address
  while (!lock_address((void*)dest, scope)) {
  }

  host_atomic_thread_fence(MemoryOrderAcquire(), scope);
  T return_val = *dest;
  *dest = op.apply(return_val, val);
  host_atomic_thread_fence(MemoryOrderRelease(), scope);
  unlock_address((void*)dest, scope);
  return return_val;
}

}  // namespace Impl
}  // namespace desul

#endif
