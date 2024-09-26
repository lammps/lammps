/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_FETCH_OP_GENERIC_HPP_
#define DESUL_ATOMICS_FETCH_OP_GENERIC_HPP_

#include <desul/atomics/Common.hpp>
#include <desul/atomics/Lock_Based_Fetch_Op.hpp>
#include <desul/atomics/Lock_Free_Fetch_Op.hpp>
#include <desul/atomics/Operator_Function_Objects.hpp>
#include <type_traits>

namespace desul {
namespace Impl {

#define DESUL_IMPL_ATOMIC_FETCH_OP(ANNOTATION, HOST_OR_DEVICE, FETCH_OP, OP_FETCH) \
  template <class T, class MemoryOrder, class MemoryScope>                         \
  ANNOTATION T HOST_OR_DEVICE##_atomic_##FETCH_OP(                                 \
      T* const dest, const T val, MemoryOrder order, MemoryScope scope) {          \
    return HOST_OR_DEVICE##_atomic_fetch_oper(                                     \
        OP_FETCH##_operator<T, const T>(), dest, val, order, scope);               \
  }                                                                                \
  template <class T, class MemoryOrder, class MemoryScope>                         \
  ANNOTATION T HOST_OR_DEVICE##_atomic_##OP_FETCH(                                 \
      T* const dest, const T val, MemoryOrder order, MemoryScope scope) {          \
    return HOST_OR_DEVICE##_atomic_oper_fetch(                                     \
        OP_FETCH##_operator<T, const T>(), dest, val, order, scope);               \
  }

#define DESUL_IMPL_ATOMIC_FETCH_OP_HOST_AND_DEVICE(FETCH_OP, OP_FETCH)           \
  DESUL_IMPL_ATOMIC_FETCH_OP(DESUL_IMPL_HOST_FUNCTION, host, FETCH_OP, OP_FETCH) \
  DESUL_IMPL_ATOMIC_FETCH_OP(DESUL_IMPL_DEVICE_FUNCTION, device, FETCH_OP, OP_FETCH)

DESUL_IMPL_ATOMIC_FETCH_OP_HOST_AND_DEVICE(fetch_add, add_fetch)
DESUL_IMPL_ATOMIC_FETCH_OP_HOST_AND_DEVICE(fetch_sub, sub_fetch)
DESUL_IMPL_ATOMIC_FETCH_OP_HOST_AND_DEVICE(fetch_max, max_fetch)
DESUL_IMPL_ATOMIC_FETCH_OP_HOST_AND_DEVICE(fetch_min, min_fetch)
DESUL_IMPL_ATOMIC_FETCH_OP_HOST_AND_DEVICE(fetch_mul, mul_fetch)
DESUL_IMPL_ATOMIC_FETCH_OP_HOST_AND_DEVICE(fetch_div, div_fetch)
DESUL_IMPL_ATOMIC_FETCH_OP_HOST_AND_DEVICE(fetch_mod, mod_fetch)
DESUL_IMPL_ATOMIC_FETCH_OP_HOST_AND_DEVICE(fetch_and, and_fetch)
DESUL_IMPL_ATOMIC_FETCH_OP_HOST_AND_DEVICE(fetch_or, or_fetch)
DESUL_IMPL_ATOMIC_FETCH_OP_HOST_AND_DEVICE(fetch_xor, xor_fetch)
DESUL_IMPL_ATOMIC_FETCH_OP_HOST_AND_DEVICE(fetch_nand, nand_fetch)

DESUL_IMPL_ATOMIC_FETCH_OP_HOST_AND_DEVICE(fetch_inc_mod, inc_mod_fetch)
DESUL_IMPL_ATOMIC_FETCH_OP_HOST_AND_DEVICE(fetch_dec_mod, dec_mod_fetch)

#undef DESUL_IMPL_ATOMIC_FETCH_OP_HOST_AND_DEVICE
#undef DESUL_IMPL_ATOMIC_FETCH_OP

#define DESUL_IMPL_ATOMIC_FETCH_OP_SHIFT(ANNOTATION, HOST_OR_DEVICE, OP)             \
  template <class T, class MemoryOrder, class MemoryScope>                           \
  ANNOTATION T HOST_OR_DEVICE##_atomic_fetch_##OP(                                   \
      T* const dest, const unsigned int val, MemoryOrder order, MemoryScope scope) { \
    return HOST_OR_DEVICE##_atomic_fetch_oper(                                       \
        OP##_fetch_operator<T, const unsigned int>(), dest, val, order, scope);      \
  }                                                                                  \
  template <class T, class MemoryOrder, class MemoryScope>                           \
  ANNOTATION T HOST_OR_DEVICE##_atomic_##OP##_fetch(                                 \
      T* const dest, const unsigned int val, MemoryOrder order, MemoryScope scope) { \
    return HOST_OR_DEVICE##_atomic_oper_fetch(                                       \
        OP##_fetch_operator<T, const unsigned int>(), dest, val, order, scope);      \
  }

#define DESUL_IMPL_ATOMIC_FETCH_OP_SHIFT_HOST_AND_DEVICE(OP)           \
  DESUL_IMPL_ATOMIC_FETCH_OP_SHIFT(DESUL_IMPL_HOST_FUNCTION, host, OP) \
  DESUL_IMPL_ATOMIC_FETCH_OP_SHIFT(DESUL_IMPL_DEVICE_FUNCTION, device, OP)

DESUL_IMPL_ATOMIC_FETCH_OP_SHIFT_HOST_AND_DEVICE(lshift)
DESUL_IMPL_ATOMIC_FETCH_OP_SHIFT_HOST_AND_DEVICE(rshift)

#undef DESUL_IMPL_ATOMIC_FETCH_OP_SHIFT_HOST_AND_DEVICE
#undef DESUL_IMPL_ATOMIC_FETCH_OP_SHIFT

// NOTE: using atomic_oper_fetch in the fallback implementation of atomic_store to avoid
// reading potentially uninitialized values which would yield undefined behavior.
#define DESUL_IMPL_ATOMIC_LOAD_AND_STORE(ANNOTATION, HOST_OR_DEVICE)                 \
  template <class T, class MemoryOrder, class MemoryScope>                           \
  ANNOTATION T HOST_OR_DEVICE##_atomic_load(                                         \
      const T* const dest, MemoryOrder order, MemoryScope scope) {                   \
    return HOST_OR_DEVICE##_atomic_fetch_oper(                                       \
        load_fetch_operator<T, const T>(), const_cast<T*>(dest), T(), order, scope); \
  }                                                                                  \
                                                                                     \
  template <class T, class MemoryOrder, class MemoryScope>                           \
  ANNOTATION void HOST_OR_DEVICE##_atomic_store(                                     \
      T* const dest, const T val, MemoryOrder order, MemoryScope scope) {            \
    (void)HOST_OR_DEVICE##_atomic_oper_fetch(                                        \
        store_fetch_operator<T, const T>(), dest, val, order, scope);                \
  }

DESUL_IMPL_ATOMIC_LOAD_AND_STORE(DESUL_IMPL_HOST_FUNCTION, host)
DESUL_IMPL_ATOMIC_LOAD_AND_STORE(DESUL_IMPL_DEVICE_FUNCTION, device)

#undef DESUL_IMPL_ATOMIC_LOAD_AND_STORE

#define DESUL_IMPL_ATOMIC_OP(ANNOTATION, HOST_OR_DEVICE, OP)              \
  template <class T, class MemoryOrder, class MemoryScope>                \
  ANNOTATION void HOST_OR_DEVICE##_atomic_##OP(                           \
      T* const dest, const T val, MemoryOrder order, MemoryScope scope) { \
    (void)HOST_OR_DEVICE##_atomic_fetch_##OP(dest, val, order, scope);    \
  }

#define DESUL_IMPL_ATOMIC_OP_HOST_AND_DEVICE(OP)           \
  DESUL_IMPL_ATOMIC_OP(DESUL_IMPL_HOST_FUNCTION, host, OP) \
  DESUL_IMPL_ATOMIC_OP(DESUL_IMPL_DEVICE_FUNCTION, device, OP)

DESUL_IMPL_ATOMIC_OP_HOST_AND_DEVICE(add)
DESUL_IMPL_ATOMIC_OP_HOST_AND_DEVICE(sub)
DESUL_IMPL_ATOMIC_OP_HOST_AND_DEVICE(mul)
DESUL_IMPL_ATOMIC_OP_HOST_AND_DEVICE(div)
DESUL_IMPL_ATOMIC_OP_HOST_AND_DEVICE(min)
DESUL_IMPL_ATOMIC_OP_HOST_AND_DEVICE(max)

#undef DESUL_IMPL_ATOMIC_OP_HOST_AND_DEVICE
#undef DESUL_IMPL_ATOMIC_OP

#define DESUL_IMPL_ATOMIC_INCREMENT_DECREMENT(ANNOTATION, HOST_OR_DEVICE) \
  template <class T, class MemoryOrder, class MemoryScope>                \
  ANNOTATION T HOST_OR_DEVICE##_atomic_inc_fetch(                         \
      T* const dest, MemoryOrder order, MemoryScope scope) {              \
    return HOST_OR_DEVICE##_atomic_add_fetch(dest, T(1), order, scope);   \
  }                                                                       \
                                                                          \
  template <class T, class MemoryOrder, class MemoryScope>                \
  ANNOTATION T HOST_OR_DEVICE##_atomic_dec_fetch(                         \
      T* const dest, MemoryOrder order, MemoryScope scope) {              \
    return HOST_OR_DEVICE##_atomic_sub_fetch(dest, T(1), order, scope);   \
  }                                                                       \
                                                                          \
  template <class T, class MemoryOrder, class MemoryScope>                \
  ANNOTATION T HOST_OR_DEVICE##_atomic_fetch_inc(                         \
      T* const dest, MemoryOrder order, MemoryScope scope) {              \
    return HOST_OR_DEVICE##_atomic_fetch_add(dest, T(1), order, scope);   \
  }                                                                       \
                                                                          \
  template <class T, class MemoryOrder, class MemoryScope>                \
  ANNOTATION T HOST_OR_DEVICE##_atomic_fetch_dec(                         \
      T* const dest, MemoryOrder order, MemoryScope scope) {              \
    return HOST_OR_DEVICE##_atomic_fetch_sub(dest, T(1), order, scope);   \
  }                                                                       \
                                                                          \
  template <class T, class MemoryOrder, class MemoryScope>                \
  ANNOTATION void HOST_OR_DEVICE##_atomic_inc(                            \
      T* const dest, MemoryOrder order, MemoryScope scope) {              \
    return HOST_OR_DEVICE##_atomic_add(dest, T(1), order, scope);         \
  }                                                                       \
                                                                          \
  template <class T, class MemoryOrder, class MemoryScope>                \
  ANNOTATION void HOST_OR_DEVICE##_atomic_dec(                            \
      T* const dest, MemoryOrder order, MemoryScope scope) {              \
    return HOST_OR_DEVICE##_atomic_sub(dest, T(1), order, scope);         \
  }

DESUL_IMPL_ATOMIC_INCREMENT_DECREMENT(DESUL_IMPL_HOST_FUNCTION, host)
DESUL_IMPL_ATOMIC_INCREMENT_DECREMENT(DESUL_IMPL_DEVICE_FUNCTION, device)

#undef DESUL_IMPL_ATOMIC_INCREMENT_DECREMENT

}  // namespace Impl
}  // namespace desul

#endif
