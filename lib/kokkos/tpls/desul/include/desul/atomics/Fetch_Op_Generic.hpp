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

#define DESUL_IMPL_ATOMIC_FETCH_OP(ANNOTATION, HOST_OR_DEVICE, _OP)        \
  template <class T, class MemoryOrder, class MemoryScope>                 \
  ANNOTATION T HOST_OR_DEVICE##_atomic_fetch##_OP(                         \
      T* const dest, const T val, MemoryOrder order, MemoryScope scope) {  \
    return HOST_OR_DEVICE##_atomic_fetch_oper(                             \
        _OP##_fetch_operator<T, const T>(), dest, val, order, scope);      \
  }                                                                        \
  template <class T, class MemoryOrder, class MemoryScope>                 \
  ANNOTATION T HOST_OR_DEVICE##_atomic##_OP##_fetch(                       \
      T* const dest, const T val, MemoryOrder order, MemoryScope scope) {  \
    return _OP##_fetch_operator<T, const T>::apply(                        \
        HOST_OR_DEVICE##_atomic_fetch##_OP(dest, val, order, scope), val); \
  }                                                                        \
  template <class T, class MemoryOrder, class MemoryScope>                 \
  ANNOTATION void HOST_OR_DEVICE##_atomic##_OP(                            \
      T* const dest, const T val, MemoryOrder order, MemoryScope scope) {  \
    (void)HOST_OR_DEVICE##_atomic_fetch##_OP(dest, val, order, scope);     \
  }

#define DESUL_IMPL_ATOMIC_FETCH_OP_HOST_AND_DEVICE(_OP)           \
  DESUL_IMPL_ATOMIC_FETCH_OP(DESUL_IMPL_HOST_FUNCTION, host, _OP) \
  DESUL_IMPL_ATOMIC_FETCH_OP(DESUL_IMPL_DEVICE_FUNCTION, device, _OP)

DESUL_IMPL_ATOMIC_FETCH_OP_HOST_AND_DEVICE(_add)
DESUL_IMPL_ATOMIC_FETCH_OP_HOST_AND_DEVICE(_sub)
DESUL_IMPL_ATOMIC_FETCH_OP_HOST_AND_DEVICE(_max)
DESUL_IMPL_ATOMIC_FETCH_OP_HOST_AND_DEVICE(_min)
DESUL_IMPL_ATOMIC_FETCH_OP_HOST_AND_DEVICE(_mul)
DESUL_IMPL_ATOMIC_FETCH_OP_HOST_AND_DEVICE(_div)
DESUL_IMPL_ATOMIC_FETCH_OP_HOST_AND_DEVICE(_mod)
DESUL_IMPL_ATOMIC_FETCH_OP_HOST_AND_DEVICE(_and)
DESUL_IMPL_ATOMIC_FETCH_OP_HOST_AND_DEVICE(_or)
DESUL_IMPL_ATOMIC_FETCH_OP_HOST_AND_DEVICE(_xor)
DESUL_IMPL_ATOMIC_FETCH_OP_HOST_AND_DEVICE(_nand)

DESUL_IMPL_ATOMIC_FETCH_OP_HOST_AND_DEVICE(_inc_mod)
DESUL_IMPL_ATOMIC_FETCH_OP_HOST_AND_DEVICE(_dec_mod)

#undef DESUL_IMPL_ATOMIC_FETCH_OP_HOST_AND_DEVICE
#undef DESUL_IMPL_ATOMIC_FETCH_OP

#define DESUL_IMPL_ATOMIC_FETCH_OP_SHIFT(ANNOTATION, HOST_OR_DEVICE, _OP)            \
  template <class T, class MemoryOrder, class MemoryScope>                           \
  ANNOTATION T HOST_OR_DEVICE##_atomic_fetch##_OP(                                   \
      T* const dest, const unsigned int val, MemoryOrder order, MemoryScope scope) { \
    return HOST_OR_DEVICE##_atomic_fetch_oper(                                       \
        _OP##_fetch_operator<T, const unsigned int>(), dest, val, order, scope);     \
  }                                                                                  \
  template <class T, class MemoryOrder, class MemoryScope>                           \
  ANNOTATION T HOST_OR_DEVICE##_atomic##_OP##_fetch(                                 \
      T* const dest, const unsigned int val, MemoryOrder order, MemoryScope scope) { \
    return _OP##_fetch_operator<T, const unsigned int>::apply(                       \
        HOST_OR_DEVICE##_atomic_fetch##_OP(dest, val, order, scope), val);           \
  }                                                                                  \
  template <class T, class MemoryOrder, class MemoryScope>                           \
  ANNOTATION void HOST_OR_DEVICE##_atomic##_OP(                                      \
      T* const dest, const unsigned int val, MemoryOrder order, MemoryScope scope) { \
    (void)HOST_OR_DEVICE##_atomic_fetch##_OP(dest, val, order, scope);               \
  }

#define DESUL_IMPL_ATOMIC_FETCH_OP_SHIFT_HOST_AND_DEVICE(_OP)           \
  DESUL_IMPL_ATOMIC_FETCH_OP_SHIFT(DESUL_IMPL_HOST_FUNCTION, host, _OP) \
  DESUL_IMPL_ATOMIC_FETCH_OP_SHIFT(DESUL_IMPL_DEVICE_FUNCTION, device, _OP)

DESUL_IMPL_ATOMIC_FETCH_OP_SHIFT_HOST_AND_DEVICE(_lshift)
DESUL_IMPL_ATOMIC_FETCH_OP_SHIFT_HOST_AND_DEVICE(_rshift)

#undef DESUL_IMPL_ATOMIC_FETCH_OP_SHIFT_HOST_AND_DEVICE
#undef DESUL_IMPL_ATOMIC_FETCH_OP_SHIFT

// NOTE: We would want to use atomic_oper_fetch in the fallback implementation of
// atomic_store to avoid reading potentially uninitialized values which would yield
// undefined behavior. As atomic_oper_fetch is not implemented, we have specializations
// of the lock based fetch_oper for _store_fetch_operator that uses a default
// constructed value instead of reading from a potentially uninitialized address.
#define DESUL_IMPL_ATOMIC_LOAD_AND_STORE(ANNOTATION, HOST_OR_DEVICE)                  \
  template <class T, class MemoryOrder, class MemoryScope>                            \
  ANNOTATION T HOST_OR_DEVICE##_atomic_load(                                          \
      const T* const dest, MemoryOrder order, MemoryScope scope) {                    \
    return HOST_OR_DEVICE##_atomic_fetch_oper(                                        \
        _load_fetch_operator<T, const T>(), const_cast<T*>(dest), T(), order, scope); \
  }                                                                                   \
                                                                                      \
  template <class T, class MemoryOrder, class MemoryScope>                            \
  ANNOTATION void HOST_OR_DEVICE##_atomic_store(                                      \
      T* const dest, const T val, MemoryOrder order, MemoryScope scope) {             \
    (void)HOST_OR_DEVICE##_atomic_fetch_oper(                                         \
        _store_fetch_operator<T, const T>(), dest, val, order, scope);                \
  }

DESUL_IMPL_ATOMIC_LOAD_AND_STORE(DESUL_IMPL_HOST_FUNCTION, host)
DESUL_IMPL_ATOMIC_LOAD_AND_STORE(DESUL_IMPL_DEVICE_FUNCTION, device)

#undef DESUL_IMPL_ATOMIC_LOAD_AND_STORE

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
