/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/
#ifndef DESUL_ATOMICS_CUDA_HPP_
#define DESUL_ATOMICS_CUDA_HPP_

#ifdef DESUL_HAVE_CUDA_ATOMICS
// When building with Clang we need to include the device functions always since Clang
// must see a consistent overload set in both device and host compilation, but that
// means we need to know on the host what to make visible, i.e. we need a host side
// compile knowledge of architecture.
#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 700)) || \
    (!defined(__NVCC__) && !defined(DESUL_CUDA_ARCH_IS_PRE_VOLTA))
#define DESUL_HAVE_CUDA_ATOMICS_ASM
#include <desul/atomics/cuda/CUDA_asm.hpp>
#endif

#if (defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 700)) || \
    (!defined(__NVCC__) && !defined(DESUL_HAVE_CUDA_ATOMICS_ASM))
namespace desul {
namespace Impl {
template <class T>
struct is_cuda_atomic_integer_type {
  static constexpr bool value = std::is_same<T, int>::value ||
                                std::is_same<T, unsigned int>::value ||
                                std::is_same<T, unsigned long long int>::value;
};

template <class T>
struct is_cuda_atomic_add_type {
  static constexpr bool value = is_cuda_atomic_integer_type<T>::value ||
#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 600)
                                std::is_same<T, double>::value ||
#endif
                                std::is_same<T, float>::value;
};

template <class T>
struct is_cuda_atomic_sub_type {
  static constexpr bool value =
      std::is_same<T, int>::value || std::is_same<T, unsigned int>::value;
};
}  // namespace Impl

// Atomic Add
template <class T>
__device__ inline
    std::enable_if_t<Impl::is_cuda_atomic_add_type<T>::value, T>
    atomic_fetch_add(T* dest, T val, MemoryOrderRelaxed, MemoryScopeDevice) {
  return atomicAdd(dest, val);
}

template <class T, class MemoryOrder>
__device__ inline
    std::enable_if_t<Impl::is_cuda_atomic_add_type<T>::value, T>
    atomic_fetch_add(T* dest, T val, MemoryOrder, MemoryScopeDevice) {
  __threadfence();
  T return_val = atomicAdd(dest, val);
  __threadfence();
  return return_val;
}

template <class T, class MemoryOrder>
__device__ inline
    std::enable_if_t<Impl::is_cuda_atomic_add_type<T>::value, T>
    atomic_fetch_add(T* dest, T val, MemoryOrder, MemoryScopeCore) {
  return atomic_fetch_add(dest, val, MemoryOrder(), MemoryScopeDevice());
}

// Atomic Sub
template <class T>
__device__ inline
    std::enable_if_t<Impl::is_cuda_atomic_sub_type<T>::value, T>
    atomic_fetch_sub(T* dest, T val, MemoryOrderRelaxed, MemoryScopeDevice) {
  return atomicSub(dest, val);
}

template <class T, class MemoryOrder>
__device__ inline
    std::enable_if_t<Impl::is_cuda_atomic_sub_type<T>::value, T>
    atomic_fetch_sub(T* dest, T val, MemoryOrder, MemoryScopeDevice) {
  __threadfence();
  T return_val = atomicSub(dest, val);
  __threadfence();
  return return_val;
}

template <class T, class MemoryOrder>
__device__ inline
    std::enable_if_t<Impl::is_cuda_atomic_sub_type<T>::value, T>
    atomic_fetch_sub(T* dest, T val, MemoryOrder, MemoryScopeCore) {
  return atomic_fetch_sub(dest, val, MemoryOrder(), MemoryScopeDevice());
}

// Wrap around atomic add
__device__ inline unsigned int atomic_fetch_inc_mod(unsigned int* dest,
                                                    unsigned int val,
                                                    MemoryOrderRelaxed,
                                                    MemoryScopeDevice) {
  return atomicInc(dest, val);
}

template <typename MemoryOrder>
__device__ inline unsigned int atomic_fetch_inc_mod(unsigned int* dest,
                                                    unsigned int val,
                                                    MemoryOrder,
                                                    MemoryScopeDevice) {
  __threadfence();
  unsigned int return_val = atomicInc(dest, val);
  __threadfence();
  return return_val;
}

template <typename MemoryOrder>
__device__ inline unsigned int atomic_fetch_inc_mod(unsigned int* dest,
                                                    unsigned int val,
                                                    MemoryOrder,
                                                    MemoryScopeCore) {
  return atomic_fetch_inc_mod(dest, val, MemoryOrder(), MemoryScopeDevice());
}

// Wrap around atomic sub
__device__ inline unsigned int atomic_fetch_dec_mod(unsigned int* dest,
                                                    unsigned int val,
                                                    MemoryOrderRelaxed,
                                                    MemoryScopeDevice) {
  return atomicDec(dest, val);
}

template <typename MemoryOrder>
__device__ inline unsigned int atomic_fetch_dec_mod(unsigned int* dest,
                                                    unsigned int val,
                                                    MemoryOrder,
                                                    MemoryScopeDevice) {
  __threadfence();
  unsigned int return_val = atomicDec(dest, val);
  __threadfence();
  return return_val;
}

template <typename MemoryOrder>
__device__ inline unsigned int atomic_fetch_dec_mod(unsigned int* dest,
                                                    unsigned int val,
                                                    MemoryOrder,
                                                    MemoryScopeCore) {
  return atomic_fetch_dec_mod(dest, val, MemoryOrder(), MemoryScopeDevice());
}

// Atomic Inc
template <typename T>
__device__ inline
    std::enable_if_t<Impl::is_cuda_atomic_add_type<T>::value, T>
    atomic_fetch_inc(T* dest, MemoryOrderRelaxed, MemoryScopeDevice) {
  return atomicAdd(dest, T(1));
}

template <typename T, typename MemoryOrder>
__device__ inline
    std::enable_if_t<Impl::is_cuda_atomic_add_type<T>::value, T>
    atomic_fetch_inc(T* dest, MemoryOrder, MemoryScopeDevice) {
  __threadfence();
  T return_val = atomicAdd(dest, T(1));
  __threadfence();

  return return_val;
}

template <typename T, typename MemoryOrder>
__device__ inline
    std::enable_if_t<Impl::is_cuda_atomic_add_type<T>::value, T>
    atomic_fetch_inc(T* dest, MemoryOrder, MemoryScopeCore) {
  return atomic_fetch_add(dest, T(1), MemoryOrder(), MemoryScopeDevice());
}

// Atomic Dec
template <typename T>
__device__ inline
    std::enable_if_t<Impl::is_cuda_atomic_sub_type<T>::value, T>
    atomic_fetch_dec(T* dest, MemoryOrderRelaxed, MemoryScopeDevice) {
  return atomicSub(dest, T(1));
}

template <typename T, typename MemoryOrder>
__device__ inline
    std::enable_if_t<Impl::is_cuda_atomic_sub_type<T>::value, T>
    atomic_fetch_dec(T* dest, MemoryOrder, MemoryScopeDevice) {
  __threadfence();
  T return_val = atomicSub(dest, T(1));
  __threadfence();
  return return_val;
}

template <typename T, typename MemoryOrder>
__device__ inline
    std::enable_if_t<Impl::is_cuda_atomic_sub_type<T>::value, T>
    atomic_fetch_dec(T* dest, MemoryOrder, MemoryScopeCore) {
  return atomic_fetch_sub(dest, T(1), MemoryOrder(), MemoryScopeDevice());
}

// Atomic Max
template <class T>
__device__ inline
    std::enable_if_t<Impl::is_cuda_atomic_integer_type<T>::value, T>
    atomic_fetch_max(T* dest, T val, MemoryOrderRelaxed, MemoryScopeDevice) {
  return atomicMax(dest, val);
}

template <class T, class MemoryOrder>
__device__ inline
    std::enable_if_t<Impl::is_cuda_atomic_integer_type<T>::value, T>
    atomic_fetch_max(T* dest, T val, MemoryOrder, MemoryScopeDevice) {
  __threadfence();
  T return_val = atomicMax(dest, val);
  __threadfence();
  return return_val;
}

template <class T, class MemoryOrder>
__device__ inline
    std::enable_if_t<Impl::is_cuda_atomic_integer_type<T>::value, T>
    atomic_fetch_max(T* dest, T val, MemoryOrder, MemoryScopeCore) {
  return atomic_fetch_max(dest, val, MemoryOrder(), MemoryScopeDevice());
}

// Atomic Min
template <class T>
__device__ inline
    std::enable_if_t<Impl::is_cuda_atomic_integer_type<T>::value, T>
    atomic_fetch_min(T* dest, T val, MemoryOrderRelaxed, MemoryScopeDevice) {
  return atomicMin(dest, val);
}

template <class T, class MemoryOrder>
__device__ inline
    std::enable_if_t<Impl::is_cuda_atomic_integer_type<T>::value, T>
    atomic_fetch_min(T* dest, T val, MemoryOrder, MemoryScopeDevice) {
  __threadfence();
  T return_val = atomicMin(dest, val);
  __threadfence();
  return return_val;
}

template <class T, class MemoryOrder>
__device__ inline
    std::enable_if_t<Impl::is_cuda_atomic_integer_type<T>::value, T>
    atomic_fetch_min(T* dest, T val, MemoryOrder, MemoryScopeCore) {
  return atomic_fetch_min(dest, val, MemoryOrder(), MemoryScopeDevice());
}

// Atomic And
template <class T>
__device__ inline
    std::enable_if_t<Impl::is_cuda_atomic_integer_type<T>::value, T>
    atomic_fetch_and(T* dest, T val, MemoryOrderRelaxed, MemoryScopeDevice) {
  return atomicAnd(dest, val);
}

template <class T, class MemoryOrder>
__device__ inline
    std::enable_if_t<Impl::is_cuda_atomic_integer_type<T>::value, T>
    atomic_fetch_and(T* dest, T val, MemoryOrder, MemoryScopeDevice) {
  __threadfence();
  T return_val = atomicAnd(dest, val);
  __threadfence();
  return return_val;
}

template <class T, class MemoryOrder>
__device__ inline
    std::enable_if_t<Impl::is_cuda_atomic_integer_type<T>::value, T>
    atomic_fetch_and(T* dest, T val, MemoryOrder, MemoryScopeCore) {
  return atomic_fetch_and(dest, val, MemoryOrder(), MemoryScopeDevice());
}

// Atomic XOR
template <class T>
__device__ inline
    std::enable_if_t<Impl::is_cuda_atomic_integer_type<T>::value, T>
    atomic_fetch_xor(T* dest, T val, MemoryOrderRelaxed, MemoryScopeDevice) {
  return atomicXor(dest, val);
}

template <class T, class MemoryOrder>
__device__ inline
    std::enable_if_t<Impl::is_cuda_atomic_integer_type<T>::value, T>
    atomic_fetch_xor(T* dest, T val, MemoryOrder, MemoryScopeDevice) {
  __threadfence();
  T return_val = atomicXor(dest, val);
  __threadfence();
  return return_val;
}

template <class T, class MemoryOrder>
__device__ inline
    std::enable_if_t<Impl::is_cuda_atomic_integer_type<T>::value, T>
    atomic_fetch_xor(T* dest, T val, MemoryOrder, MemoryScopeCore) {
  return atomic_fetch_xor(dest, val, MemoryOrder(), MemoryScopeDevice());
}

// Atomic OR
template <class T>
__device__ inline
    std::enable_if_t<Impl::is_cuda_atomic_integer_type<T>::value, T>
    atomic_fetch_or(T* dest, T val, MemoryOrderRelaxed, MemoryScopeDevice) {
  return atomicOr(dest, val);
}

template <class T, class MemoryOrder>
__device__ inline
    std::enable_if_t<Impl::is_cuda_atomic_integer_type<T>::value, T>
    atomic_fetch_or(T* dest, T val, MemoryOrder, MemoryScopeDevice) {
  __threadfence();
  T return_val = atomicOr(dest, val);
  __threadfence();
  return return_val;
}

template <class T, class MemoryOrder>
__device__ inline
    std::enable_if_t<Impl::is_cuda_atomic_integer_type<T>::value, T>
    atomic_fetch_or(T* dest, T val, MemoryOrder, MemoryScopeCore) {
  return atomic_fetch_or(dest, val, MemoryOrder(), MemoryScopeDevice());
}
}  // namespace desul
#endif

#if !defined(__NVCC__)
// Functions defined as device functions in CUDA which don't exist in the GCC overload
// set
namespace desul {

#if defined(DESUL_HAVE_CUDA_ATOMICS_ASM)
#define DESUL_IMPL_CUDA_HOST_ATOMIC_ADD(TYPE, ORDER, SCOPE)                      \
  inline void atomic_add(TYPE* const dest, TYPE val, ORDER order, SCOPE scope) { \
    (void)atomic_fetch_add(dest, val, order, scope);                             \
  }
DESUL_IMPL_CUDA_HOST_ATOMIC_ADD(int32_t, MemoryOrderRelaxed, MemoryScopeDevice);
DESUL_IMPL_CUDA_HOST_ATOMIC_ADD(long,
                                MemoryOrderRelaxed,
                                MemoryScopeDevice);  // only for ASM?
DESUL_IMPL_CUDA_HOST_ATOMIC_ADD(unsigned int, MemoryOrderRelaxed, MemoryScopeDevice);
DESUL_IMPL_CUDA_HOST_ATOMIC_ADD(unsigned long long,
                                MemoryOrderRelaxed,
                                MemoryScopeDevice);
DESUL_IMPL_CUDA_HOST_ATOMIC_ADD(float, MemoryOrderRelaxed, MemoryScopeDevice);
DESUL_IMPL_CUDA_HOST_ATOMIC_ADD(double, MemoryOrderRelaxed, MemoryScopeDevice);

#define DESUL_IMPL_CUDA_HOST_ATOMIC_SUB(TYPE, ORDER, SCOPE)                      \
  inline void atomic_sub(TYPE* const dest, TYPE val, ORDER order, SCOPE scope) { \
    (void)atomic_fetch_sub(dest, val, order, scope);                             \
  }
DESUL_IMPL_CUDA_HOST_ATOMIC_SUB(int32_t, MemoryOrderRelaxed, MemoryScopeDevice);
DESUL_IMPL_CUDA_HOST_ATOMIC_SUB(long,
                                MemoryOrderRelaxed,
                                MemoryScopeDevice);  // only for ASM?
DESUL_IMPL_CUDA_HOST_ATOMIC_SUB(unsigned int, MemoryOrderRelaxed, MemoryScopeDevice);
DESUL_IMPL_CUDA_HOST_ATOMIC_SUB(float, MemoryOrderRelaxed, MemoryScopeDevice);
DESUL_IMPL_CUDA_HOST_ATOMIC_SUB(double, MemoryOrderRelaxed, MemoryScopeDevice);

#define DESUL_IMPL_CUDA_HOST_ATOMIC_INC(TYPE, ORDER, SCOPE)            \
  inline void atomic_inc(TYPE* const dest, ORDER order, SCOPE scope) { \
    (void)atomic_fetch_inc(dest, order, scope);                        \
  }
DESUL_IMPL_CUDA_HOST_ATOMIC_INC(unsigned int,
                                MemoryOrderRelaxed,
                                MemoryScopeDevice);  // only for ASM?

#define DESUL_IMPL_CUDA_HOST_ATOMIC_DEC(TYPE, ORDER, SCOPE)            \
  inline void atomic_dec(TYPE* const dest, ORDER order, SCOPE scope) { \
    (void)atomic_fetch_dec(dest, order, scope);                        \
  }
DESUL_IMPL_CUDA_HOST_ATOMIC_DEC(unsigned,
                                MemoryOrderRelaxed,
                                MemoryScopeDevice);  // only for ASM?

#endif  // DESUL_HAVE_CUDA_ATOMICS_ASM

#define DESUL_IMPL_CUDA_HOST_ATOMIC_INC_MOD(TYPE, ORDER, SCOPE)                      \
  inline TYPE atomic_fetch_inc_mod(TYPE* dest, TYPE val, ORDER order, SCOPE scope) { \
    using cas_t = typename Impl::atomic_compare_exchange_type<sizeof(TYPE)>::type;   \
    cas_t oldval = reinterpret_cast<cas_t&>(*dest);                                  \
    cas_t assume = oldval;                                                           \
    do {                                                                             \
      assume = oldval;                                                               \
      TYPE newval = (reinterpret_cast<TYPE&>(assume) >= val)                         \
                        ? static_cast<TYPE>(0)                                       \
                        : reinterpret_cast<TYPE&>(assume) + static_cast<TYPE>(1);    \
      oldval = desul::atomic_compare_exchange(reinterpret_cast<cas_t*>(dest),        \
                                              assume,                                \
                                              reinterpret_cast<cas_t&>(newval),      \
                                              order,                                 \
                                              scope);                                \
    } while (assume != oldval);                                                      \
    return reinterpret_cast<TYPE&>(oldval);                                          \
  }
DESUL_IMPL_CUDA_HOST_ATOMIC_INC_MOD(unsigned int,
                                    MemoryOrderRelaxed,
                                    MemoryScopeDevice);
#define DESUL_IMPL_CUDA_HOST_ATOMIC_DEC_MOD(TYPE, ORDER, SCOPE)                      \
  inline TYPE atomic_fetch_dec_mod(TYPE* dest, TYPE val, ORDER order, SCOPE scope) { \
    using cas_t = typename Impl::atomic_compare_exchange_type<sizeof(TYPE)>::type;   \
    cas_t oldval = reinterpret_cast<cas_t&>(*dest);                                  \
    cas_t assume = oldval;                                                           \
    do {                                                                             \
      assume = oldval;                                                               \
      TYPE newval = ((reinterpret_cast<TYPE&>(assume) == static_cast<TYPE>(0)) |     \
                     (reinterpret_cast<TYPE&>(assume) > val))                        \
                        ? val                                                        \
                        : reinterpret_cast<TYPE&>(assume) - static_cast<TYPE>(1);    \
      oldval = desul::atomic_compare_exchange(reinterpret_cast<cas_t*>(dest),        \
                                              assume,                                \
                                              reinterpret_cast<cas_t&>(newval),      \
                                              order,                                 \
                                              scope);                                \
    } while (assume != oldval);                                                      \
    return reinterpret_cast<TYPE&>(oldval);                                          \
  }
DESUL_IMPL_CUDA_HOST_ATOMIC_DEC_MOD(unsigned int,
                                    MemoryOrderRelaxed,
                                    MemoryScopeDevice);

#define DESUL_IMPL_CUDA_HOST_ATOMIC_FETCH_ADD(TYPE, ORDER, SCOPE)                      \
  inline TYPE atomic_fetch_add(TYPE* const dest, TYPE val, ORDER order, SCOPE scope) { \
    return Impl::atomic_fetch_oper(                                                    \
        Impl::AddOper<TYPE, const TYPE>(), dest, val, order, scope);                   \
  }
DESUL_IMPL_CUDA_HOST_ATOMIC_FETCH_ADD(float, MemoryOrderRelaxed, MemoryScopeDevice);
DESUL_IMPL_CUDA_HOST_ATOMIC_FETCH_ADD(double, MemoryOrderRelaxed, MemoryScopeDevice);

#define DESUL_IMPL_CUDA_HOST_ATOMIC_FETCH_SUB(TYPE, ORDER, SCOPE)                      \
  inline TYPE atomic_fetch_sub(TYPE* const dest, TYPE val, ORDER order, SCOPE scope) { \
    return Impl::atomic_fetch_oper(                                                    \
        Impl::SubOper<TYPE, const TYPE>(), dest, val, order, scope);                   \
  }
DESUL_IMPL_CUDA_HOST_ATOMIC_FETCH_SUB(float, MemoryOrderRelaxed, MemoryScopeDevice);
DESUL_IMPL_CUDA_HOST_ATOMIC_FETCH_SUB(double, MemoryOrderRelaxed, MemoryScopeDevice);

#define DESUL_IMPL_CUDA_HOST_ATOMIC_FETCH_MAX(TYPE, ORDER, SCOPE)                      \
  inline TYPE atomic_fetch_max(TYPE* const dest, TYPE val, ORDER order, SCOPE scope) { \
    return Impl::atomic_fetch_oper(                                                    \
        Impl::MaxOper<TYPE, const TYPE>(), dest, val, order, scope);                   \
  }
DESUL_IMPL_CUDA_HOST_ATOMIC_FETCH_MAX(int, MemoryOrderRelaxed, MemoryScopeDevice);
DESUL_IMPL_CUDA_HOST_ATOMIC_FETCH_MAX(long,
                                      MemoryOrderRelaxed,
                                      MemoryScopeDevice);  // only for ASM?
DESUL_IMPL_CUDA_HOST_ATOMIC_FETCH_MAX(unsigned int,
                                      MemoryOrderRelaxed,
                                      MemoryScopeDevice);
DESUL_IMPL_CUDA_HOST_ATOMIC_FETCH_MAX(unsigned long,
                                      MemoryOrderRelaxed,
                                      MemoryScopeDevice);
//  DESUL_IMPL_CUDA_HOST_ATOMIC_FETCH_MAX(unsigned long
//  long,MemoryOrderRelaxed,MemoryScopeDevice);

#define DESUL_IMPL_CUDA_HOST_ATOMIC_FETCH_MIN(TYPE, ORDER, SCOPE)                      \
  inline TYPE atomic_fetch_min(TYPE* const dest, TYPE val, ORDER order, SCOPE scope) { \
    return Impl::atomic_fetch_oper(                                                    \
        Impl::MinOper<TYPE, const TYPE>(), dest, val, order, scope);                   \
  }
DESUL_IMPL_CUDA_HOST_ATOMIC_FETCH_MIN(int, MemoryOrderRelaxed, MemoryScopeDevice);
DESUL_IMPL_CUDA_HOST_ATOMIC_FETCH_MIN(long,
                                      MemoryOrderRelaxed,
                                      MemoryScopeDevice);  // only for ASM?
DESUL_IMPL_CUDA_HOST_ATOMIC_FETCH_MIN(unsigned int,
                                      MemoryOrderRelaxed,
                                      MemoryScopeDevice);
DESUL_IMPL_CUDA_HOST_ATOMIC_FETCH_MIN(unsigned long,
                                      MemoryOrderRelaxed,
                                      MemoryScopeDevice);
//  DESUL_IMPL_CUDA_HOST_ATOMIC_FETCH_MIN(unsigned long
//  long,MemoryOrderRelaxed,MemoryScopeDevice); inline void atomic_fetch_max(int32_t*
//  const dest, int32_t val, MemoryOrderRelaxed order, MemoryScopeDevice scope) {

}  // namespace desul

// Functions defined int the GCC overload set but not in the device overload set
namespace desul {
__device__ inline unsigned long long atomic_fetch_add(unsigned long long* const dest,
                                                      unsigned long long val,
                                                      MemoryOrderRelaxed order,
                                                      MemoryScopeDevice scope) {
  return Impl::atomic_fetch_oper(
      Impl::AddOper<unsigned long long, const unsigned long long>(),
      dest,
      val,
      order,
      scope);
}
__device__ inline long long atomic_fetch_add(long long* const dest,
                                             long long val,
                                             MemoryOrderRelaxed order,
                                             MemoryScopeDevice scope) {
  return Impl::atomic_fetch_oper(
      Impl::AddOper<long long, const long long>(), dest, val, order, scope);
}
__device__ inline long atomic_fetch_add(long* const dest,
                                        long val,
                                        MemoryOrderRelaxed order,
                                        MemoryScopeDevice scope) {
  return Impl::atomic_fetch_oper(
      Impl::AddOper<long, const long>(), dest, val, order, scope);
}
__device__ inline long long atomic_fetch_sub(long long* const dest,
                                             long long val,
                                             MemoryOrderRelaxed order,
                                             MemoryScopeDevice scope) {
  return Impl::atomic_fetch_oper(
      Impl::SubOper<long long, const long long>(), dest, val, order, scope);
}
__device__ inline long atomic_fetch_sub(long* const dest,
                                        long val,
                                        MemoryOrderRelaxed order,
                                        MemoryScopeDevice scope) {
  return Impl::atomic_fetch_oper(
      Impl::SubOper<long, const long>(), dest, val, order, scope);
}
__device__ inline long atomic_fetch_max(long* const dest,
                                        long val,
                                        MemoryOrderRelaxed order,
                                        MemoryScopeDevice scope) {
  return Impl::atomic_fetch_oper(
      Impl::MaxOper<long, const long>(), dest, val, order, scope);
}
__device__ inline long atomic_fetch_min(long* const dest,
                                        long val,
                                        MemoryOrderRelaxed order,
                                        MemoryScopeDevice scope) {
  return Impl::atomic_fetch_oper(
      Impl::MinOper<long, const long>(), dest, val, order, scope);
}
__device__ inline long atomic_fetch_or(long* const dest,
                                       long val,
                                       MemoryOrderRelaxed order,
                                       MemoryScopeDevice scope) {
  return Impl::atomic_fetch_oper(
      Impl::OrOper<long, const long>(), dest, val, order, scope);
}
__device__ inline long long atomic_fetch_or(long long* const dest,
                                            long long val,
                                            MemoryOrderRelaxed order,
                                            MemoryScopeDevice scope) {
  return Impl::atomic_fetch_oper(
      Impl::OrOper<long long, const long long>(), dest, val, order, scope);
}
__device__ inline long atomic_fetch_xor(long* const dest,
                                        long val,
                                        MemoryOrderRelaxed order,
                                        MemoryScopeDevice scope) {
  return Impl::atomic_fetch_oper(
      Impl::XorOper<long, const long>(), dest, val, order, scope);
}
__device__ inline long long atomic_fetch_xor(long long* const dest,
                                             long long val,
                                             MemoryOrderRelaxed order,
                                             MemoryScopeDevice scope) {
  return Impl::atomic_fetch_oper(
      Impl::XorOper<long long, const long long>(), dest, val, order, scope);
}
__device__ inline long atomic_fetch_and(long* const dest,
                                        long val,
                                        MemoryOrderRelaxed order,
                                        MemoryScopeDevice scope) {
  return Impl::atomic_fetch_oper(
      Impl::AndOper<long, const long>(), dest, val, order, scope);
}
__device__ inline long long atomic_fetch_and(long long* const dest,
                                             long long val,
                                             MemoryOrderRelaxed order,
                                             MemoryScopeDevice scope) {
  return Impl::atomic_fetch_oper(
      Impl::AndOper<long long, const long long>(), dest, val, order, scope);
}

__device__ inline unsigned long long atomic_add_fetch(unsigned long long* const dest,
                                                      unsigned long long val,
                                                      MemoryOrderRelaxed order,
                                                      MemoryScopeDevice scope) {
  return Impl::atomic_oper_fetch(
      Impl::AddOper<unsigned long long, const unsigned long long>(),
      dest,
      val,
      order,
      scope);
}
__device__ inline long long atomic_add_fetch(long long* const dest,
                                             long long val,
                                             MemoryOrderRelaxed order,
                                             MemoryScopeDevice scope) {
  return Impl::atomic_oper_fetch(
      Impl::AddOper<long long, const long long>(), dest, val, order, scope);
}
__device__ inline long atomic_add_fetch(long* const dest,
                                        long val,
                                        MemoryOrderRelaxed order,
                                        MemoryScopeDevice scope) {
  return Impl::atomic_oper_fetch(
      Impl::AddOper<long, const long>(), dest, val, order, scope);
}
__device__ inline long long atomic_sub_fetch(long long* const dest,
                                             long long val,
                                             MemoryOrderRelaxed order,
                                             MemoryScopeDevice scope) {
  return Impl::atomic_oper_fetch(
      Impl::SubOper<long long, const long long>(), dest, val, order, scope);
}
__device__ inline long atomic_sub_fetch(long* const dest,
                                        long val,
                                        MemoryOrderRelaxed order,
                                        MemoryScopeDevice scope) {
  return Impl::atomic_oper_fetch(
      Impl::SubOper<long, const long>(), dest, val, order, scope);
}
__device__ inline long long atomic_or_fetch(long long* const dest,
                                            long long val,
                                            MemoryOrderRelaxed order,
                                            MemoryScopeDevice scope) {
  return Impl::atomic_oper_fetch(
      Impl::OrOper<long long, const long long>(), dest, val, order, scope);
}
__device__ inline long atomic_or_fetch(long* const dest,
                                       long val,
                                       MemoryOrderRelaxed order,
                                       MemoryScopeDevice scope) {
  return Impl::atomic_oper_fetch(
      Impl::OrOper<long, const long>(), dest, val, order, scope);
}
__device__ inline long long atomic_xor_fetch(long long* const dest,
                                             long long val,
                                             MemoryOrderRelaxed order,
                                             MemoryScopeDevice scope) {
  return Impl::atomic_oper_fetch(
      Impl::XorOper<long long, const long long>(), dest, val, order, scope);
}
__device__ inline long atomic_xor_fetch(long* const dest,
                                        long val,
                                        MemoryOrderRelaxed order,
                                        MemoryScopeDevice scope) {
  return Impl::atomic_oper_fetch(
      Impl::XorOper<long, const long>(), dest, val, order, scope);
}
__device__ inline long long atomic_and_fetch(long long* const dest,
                                             long val,
                                             MemoryOrderRelaxed order,
                                             MemoryScopeDevice scope) {
  return Impl::atomic_oper_fetch(
      Impl::AndOper<long long, const long long>(), dest, val, order, scope);
}
__device__ inline long atomic_and_fetch(long* const dest,
                                        long val,
                                        MemoryOrderRelaxed order,
                                        MemoryScopeDevice scope) {
  return Impl::atomic_oper_fetch(
      Impl::AndOper<long, const long>(), dest, val, order, scope);
}
}  // namespace desul
#endif

#endif  // DESUL_HAVE_CUDA_ATOMICS
#endif
