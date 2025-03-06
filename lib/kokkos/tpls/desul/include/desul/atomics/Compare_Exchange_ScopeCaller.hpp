/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_COMPARE_EXCHANGE_SCOPECALLER_HPP_
#define DESUL_ATOMICS_COMPARE_EXCHANGE_SCOPECALLER_HPP_

#include <desul/atomics/Common.hpp>

namespace desul {

#define DESUL_ATOMIC_EXCHANGE_SCOPECALLER(MEMORY_ORDER)               \
  template <class T>                                                  \
  DESUL_INLINE_FUNCTION T atomic_exchange(                            \
      T* dest, T value, MEMORY_ORDER, MemoryScopeCaller) {            \
    T return_val = *dest;                                             \
    *dest = value;                                                    \
    return return_val;                                                \
  }                                                                   \
                                                                      \
  template <class T>                                                  \
  DESUL_INLINE_FUNCTION T atomic_compare_exchange(                    \
      T* dest, T compare, T value, MEMORY_ORDER, MemoryScopeCaller) { \
    T current_val = *dest;                                            \
    if (current_val == compare) *dest = value;                        \
    return current_val;                                               \
  }

DESUL_ATOMIC_EXCHANGE_SCOPECALLER(MemoryOrderSeqCst)
DESUL_ATOMIC_EXCHANGE_SCOPECALLER(MemoryOrderAcqRel)
DESUL_ATOMIC_EXCHANGE_SCOPECALLER(MemoryOrderRelease)
DESUL_ATOMIC_EXCHANGE_SCOPECALLER(MemoryOrderAcquire)
DESUL_ATOMIC_EXCHANGE_SCOPECALLER(MemoryOrderRelaxed)

#undef DESUL_ATOMIC_EXCHANGE_SCOPECALLER

}  // namespace desul

#endif
