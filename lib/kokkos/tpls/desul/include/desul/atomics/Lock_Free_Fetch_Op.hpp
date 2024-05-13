/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_LOCK_FREE_FETCH_OP_HPP_
#define DESUL_ATOMICS_LOCK_FREE_FETCH_OP_HPP_

#include <desul/atomics/Common.hpp>
#include <desul/atomics/Compare_Exchange.hpp>
#include <type_traits>

#if defined(__GNUC__) && (!defined(__clang__))
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

namespace desul {
namespace Impl {

#define DESUL_IMPL_ATOMIC_FETCH_OPER(ANNOTATION, HOST_OR_DEVICE)                     \
  template <class Oper,                                                              \
            class T,                                                                 \
            class MemoryOrder,                                                       \
            class MemoryScope,                                                       \
            std::enable_if_t<atomic_always_lock_free(sizeof(T)), int> = 0>           \
  ANNOTATION T HOST_OR_DEVICE##_atomic_fetch_oper(                                   \
      const Oper& op,                                                                \
      T* const dest,                                                                 \
      dont_deduce_this_parameter_t<const T> val,                                     \
      MemoryOrder order,                                                             \
      MemoryScope scope) {                                                           \
    using cas_t = atomic_compare_exchange_t<T>;                                      \
    cas_t oldval = reinterpret_cast<cas_t&>(*dest);                                  \
    cas_t assume = oldval;                                                           \
                                                                                     \
    do {                                                                             \
      if (check_early_exit(op, reinterpret_cast<T&>(oldval), val))                   \
        return reinterpret_cast<T&>(oldval);                                         \
      assume = oldval;                                                               \
      T newval = op.apply(reinterpret_cast<T&>(assume), val);                        \
      oldval =                                                                       \
          HOST_OR_DEVICE##_atomic_compare_exchange(reinterpret_cast<cas_t*>(dest),   \
                                                   assume,                           \
                                                   reinterpret_cast<cas_t&>(newval), \
                                                   order,                            \
                                                   scope);                           \
    } while (assume != oldval);                                                      \
                                                                                     \
    return reinterpret_cast<T&>(oldval);                                             \
  }                                                                                  \
                                                                                     \
  template <class Oper,                                                              \
            class T,                                                                 \
            class MemoryOrder,                                                       \
            class MemoryScope,                                                       \
            std::enable_if_t<atomic_always_lock_free(sizeof(T)), int> = 0>           \
  ANNOTATION T HOST_OR_DEVICE##_atomic_oper_fetch(                                   \
      const Oper& op,                                                                \
      T* const dest,                                                                 \
      dont_deduce_this_parameter_t<const T> val,                                     \
      MemoryOrder order,                                                             \
      MemoryScope scope) {                                                           \
    using cas_t = atomic_compare_exchange_t<T>;                                      \
    cas_t oldval = reinterpret_cast<cas_t&>(*dest);                                  \
    T newval = val;                                                                  \
    cas_t assume = oldval;                                                           \
    do {                                                                             \
      if (check_early_exit(op, reinterpret_cast<T&>(oldval), val))                   \
        return reinterpret_cast<T&>(oldval);                                         \
      assume = oldval;                                                               \
      newval = op.apply(reinterpret_cast<T&>(assume), val);                          \
      oldval =                                                                       \
          HOST_OR_DEVICE##_atomic_compare_exchange(reinterpret_cast<cas_t*>(dest),   \
                                                   assume,                           \
                                                   reinterpret_cast<cas_t&>(newval), \
                                                   order,                            \
                                                   scope);                           \
    } while (assume != oldval);                                                      \
                                                                                     \
    return newval;                                                                   \
  }

DESUL_IMPL_ATOMIC_FETCH_OPER(DESUL_IMPL_HOST_FUNCTION, host)
DESUL_IMPL_ATOMIC_FETCH_OPER(DESUL_IMPL_DEVICE_FUNCTION, device)

#undef DESUL_IMPL_ATOMIC_FETCH_OPER

}  // namespace Impl
}  // namespace desul

#if defined(__GNUC__) && (!defined(__clang__))
#pragma GCC diagnostic pop
#endif

#endif
