/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_COMPARE_EXCHANGE_MSVC_HPP_
#define DESUL_ATOMICS_COMPARE_EXCHANGE_MSVC_HPP_
#include <type_traits>

#include "desul/atomics/Common.hpp"
#ifdef DESUL_HAVE_MSVC_ATOMICS

#ifndef DESUL_HAVE_16BYTE_COMPARE_AND_SWAP
#define DESUL_HAVE_16BYTE_COMPARE_AND_SWAP
#endif

namespace desul {

// Forward declare these functions. They use compare_exchange themselves
// so the actual header file with them comes after this file is included.
namespace Impl {
template <typename MemoryScope>
inline bool lock_address(void* ptr, MemoryScope ms);

template <typename MemoryScope>
void unlock_address(void* ptr, MemoryScope ms);
}  // namespace Impl

template <class MemoryOrder, class MemoryScope>
void atomic_thread_fence(MemoryOrder, MemoryScope) {
  std::atomic_thread_fence(CXXMemoryOrder<MemoryOrder>::value);
}

template <typename T, class MemoryScope>
typename std::enable_if<sizeof(T) == 1, T>::type atomic_exchange(T* const dest,
                                                                 T val,
                                                                 MemoryOrderRelaxed,
                                                                 MemoryScope) {
  char return_val = _InterlockedExchange8((char*)dest, *((char*)&val));
  return *(reinterpret_cast<T*>(&return_val));
}

template <typename T, class MemoryScope>
typename std::enable_if<sizeof(T) == 2, T>::type atomic_exchange(T* const dest,
                                                                 T val,
                                                                 MemoryOrderRelaxed,
                                                                 MemoryScope) {
  short return_val = _InterlockedExchange16((short*)dest, *((short*)&val));
  return *(reinterpret_cast<T*>(&return_val));
}

template <typename T, class MemoryScope>
typename std::enable_if<sizeof(T) == 4, T>::type atomic_exchange(T* const dest,
                                                                 T val,
                                                                 MemoryOrderRelaxed,
                                                                 MemoryScope) {
  long return_val = _InterlockedExchange((long*)dest, *((long*)&val));
  return *(reinterpret_cast<T*>(&return_val));
}

template <typename T, class MemoryScope>
typename std::enable_if<sizeof(T) == 8, T>::type atomic_exchange(T* const dest,
                                                                 T val,
                                                                 MemoryOrderRelaxed,
                                                                 MemoryScope) {
  __int64 return_val = _InterlockedExchange64((__int64*)dest, *((__int64*)&val));
  return *(reinterpret_cast<T*>(&return_val));
}

template <typename T, class MemoryScope>
typename std::enable_if<sizeof(T) == 1, T>::type atomic_exchange(T* const dest,
                                                                 T val,
                                                                 MemoryOrderSeqCst,
                                                                 MemoryScope) {
  char return_val = _InterlockedExchange8((char*)dest, *((char*)&val));
  return *(reinterpret_cast<T*>(&return_val));
}

template <typename T, class MemoryScope>
typename std::enable_if<sizeof(T) == 2, T>::type atomic_exchange(T* const dest,
                                                                 T val,
                                                                 MemoryOrderSeqCst,
                                                                 MemoryScope) {
  short return_val = _InterlockedExchange16((short*)dest, *((short*)&val));
  return *(reinterpret_cast<T*>(&return_val));
}

template <typename T, class MemoryScope>
typename std::enable_if<sizeof(T) == 4, T>::type atomic_exchange(T* const dest,
                                                                 T val,
                                                                 MemoryOrderSeqCst,
                                                                 MemoryScope) {
  long return_val = _InterlockedExchange((long*)dest, *((long*)&val));
  return *(reinterpret_cast<T*>(&return_val));
}

template <typename T, class MemoryScope>
typename std::enable_if<sizeof(T) == 8, T>::type atomic_exchange(T* const dest,
                                                                 T val,
                                                                 MemoryOrderSeqCst,
                                                                 MemoryScope) {
  __int64 return_val = _InterlockedExchange64((__int64*)dest, *((__int64*)&val));
  return *(reinterpret_cast<T*>(&return_val));
}

template <typename T, class MemoryOrder, class MemoryScope>
typename std::enable_if<(sizeof(T) != 1 && sizeof(T) != 2 && sizeof(T) != 4 &&
                         sizeof(T) != 8),
                        T>::type
atomic_exchange(T* const dest, T val, MemoryOrder, MemoryScope scope) {
  while (!Impl::lock_address((void*)dest, scope)) {
  }
  if (std::is_same<MemoryOrder, MemoryOrderSeqCst>::value)
    atomic_thread_fence(MemoryOrderRelease(), scope);
  atomic_thread_fence(MemoryOrderAcquire(), scope);
  T return_val = *dest;
  *dest = val;
  atomic_thread_fence(MemoryOrderRelease(), scope);

  Impl::unlock_address((void*)dest, scope);
  return return_val;
}

template <typename T, class MemoryScope>
typename std::enable_if<sizeof(T) == 1, T>::type atomic_compare_exchange(
    T* const dest, T compare, T val, MemoryOrderRelaxed, MemoryScope) {
  char return_val =
      _InterlockedCompareExchange8((char*)dest, *((char*)&val), *((char*)&compare));
  return *(reinterpret_cast<T*>(&return_val));
}

template <typename T, class MemoryScope>
typename std::enable_if<sizeof(T) == 2, T>::type atomic_compare_exchange(
    T* const dest, T compare, T val, MemoryOrderRelaxed, MemoryScope) {
  short return_val =
      _InterlockedCompareExchange16((short*)dest, *((short*)&val), *((short*)&compare));
  return *(reinterpret_cast<T*>(&return_val));
}

template <typename T, class MemoryScope>
typename std::enable_if<sizeof(T) == 4, T>::type atomic_compare_exchange(
    T* const dest, T compare, T val, MemoryOrderRelaxed, MemoryScope) {
  long return_val =
      _InterlockedCompareExchange((long*)dest, *((long*)&val), *((long*)&compare));
  return *(reinterpret_cast<T*>(&return_val));
}

template <typename T, class MemoryScope>
typename std::enable_if<sizeof(T) == 8, T>::type atomic_compare_exchange(
    T* const dest, T compare, T val, MemoryOrderRelaxed, MemoryScope) {
  __int64 return_val = _InterlockedCompareExchange64(
      (__int64*)dest, *((__int64*)&val), *((__int64*)&compare));
  return *(reinterpret_cast<T*>(&return_val));
}

template <typename T, class MemoryScope>
typename std::enable_if<sizeof(T) == 16, T>::type atomic_compare_exchange(
    T* const dest, T compare, T val, MemoryOrderRelaxed, MemoryScope) {
  Dummy16ByteValue* val16 = reinterpret_cast<Dummy16ByteValue*>(&val);
  (void)_InterlockedCompareExchange128(reinterpret_cast<__int64*>(dest),
                                       val16->value2,
                                       val16->value1,
                                       (reinterpret_cast<__int64*>(&compare)));
  return compare;
}

template <typename T, class MemoryScope>
typename std::enable_if<sizeof(T) == 1, T>::type atomic_compare_exchange(
    T* const dest, T compare, T val, MemoryOrderSeqCst, MemoryScope) {
  char return_val =
      _InterlockedCompareExchange8((char*)dest, *((char*)&val), *((char*)&compare));
  return *(reinterpret_cast<T*>(&return_val));
}

template <typename T, class MemoryScope>
typename std::enable_if<sizeof(T) == 2, T>::type atomic_compare_exchange(
    T* const dest, T compare, T val, MemoryOrderSeqCst, MemoryScope) {
  short return_val =
      _InterlockedCompareExchange16((short*)dest, *((short*)&val), *((short*)&compare));
  return *(reinterpret_cast<T*>(&return_val));
}

template <typename T, class MemoryScope>
typename std::enable_if<sizeof(T) == 4, T>::type atomic_compare_exchange(
    T* const dest, T compare, T val, MemoryOrderSeqCst, MemoryScope) {
  long return_val =
      _InterlockedCompareExchange((long*)dest, *((long*)&val), *((long*)&compare));
  return *(reinterpret_cast<T*>(&return_val));
}

template <typename T, class MemoryScope>
typename std::enable_if<sizeof(T) == 8, T>::type atomic_compare_exchange(
    T* const dest, T compare, T val, MemoryOrderSeqCst, MemoryScope) {
  __int64 return_val = _InterlockedCompareExchange64(
      (__int64*)dest, *((__int64*)&val), *((__int64*)&compare));
  return *(reinterpret_cast<T*>(&return_val));
}

template <typename T, class MemoryScope>
typename std::enable_if<sizeof(T) == 16, T>::type atomic_compare_exchange(
    T* const dest, T compare, T val, MemoryOrderSeqCst, MemoryScope) {
  Dummy16ByteValue* val16 = reinterpret_cast<Dummy16ByteValue*>(&val);
  (void)_InterlockedCompareExchange128(reinterpret_cast<__int64*>(dest),
                                       val16->value2,
                                       val16->value1,
                                       (reinterpret_cast<__int64*>(&compare)));
  return compare;
}

template <typename T, class MemoryOrder, class MemoryScope>
typename std::enable_if<(sizeof(T) != 1 && sizeof(T) != 2 && sizeof(T) != 4 &&
                         sizeof(T) != 8 && sizeof(T) != 16),
                        T>::type
atomic_compare_exchange(
    T* const dest, T compare, T val, MemoryOrder, MemoryScope scope) {
  while (!Impl::lock_address((void*)dest, scope)) {
  }
  if (std::is_same<MemoryOrder, MemoryOrderSeqCst>::value)
    atomic_thread_fence(MemoryOrderRelease(), scope);
  atomic_thread_fence(MemoryOrderAcquire(), scope);
  T return_val = *dest;
  if (return_val == compare) {
    *dest = val;
    atomic_thread_fence(MemoryOrderRelease(), scope);
  }

  Impl::unlock_address((void*)dest, scope);
  return return_val;
}

}  // namespace desul

#endif
#endif
