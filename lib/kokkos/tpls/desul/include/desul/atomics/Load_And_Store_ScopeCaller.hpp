/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_LOAD_AND_STORE_SCOPECALLER_HPP_
#define DESUL_ATOMICS_LOAD_AND_STORE_SCOPECALLER_HPP_

#include <desul/atomics/Common.hpp>
#include <desul/atomics/Macros.hpp>

namespace desul {

#define DESUL_IMPL_LOAD_AND_STORE_SCOPECALLER(MEMORY_ORDER)                            \
  template <class T>                                                                   \
  DESUL_INLINE_FUNCTION T atomic_load(T const* ptr, MEMORY_ORDER, MemoryScopeCaller) { \
    return *ptr;                                                                       \
  }                                                                                    \
                                                                                       \
  template <class T>                                                                   \
  DESUL_INLINE_FUNCTION void atomic_store(                                             \
      T* ptr, T val, MEMORY_ORDER, MemoryScopeCaller) {                                \
    *ptr = val;                                                                        \
  }

DESUL_IMPL_LOAD_AND_STORE_SCOPECALLER(MemoryOrderSeqCst)
DESUL_IMPL_LOAD_AND_STORE_SCOPECALLER(MemoryOrderAcqRel)
DESUL_IMPL_LOAD_AND_STORE_SCOPECALLER(MemoryOrderRelease)
DESUL_IMPL_LOAD_AND_STORE_SCOPECALLER(MemoryOrderAcquire)
DESUL_IMPL_LOAD_AND_STORE_SCOPECALLER(MemoryOrderRelaxed)

#undef DESUL_IMPL_LOAD_AND_STORE_SCOPECALLER

}  // namespace desul

#endif
