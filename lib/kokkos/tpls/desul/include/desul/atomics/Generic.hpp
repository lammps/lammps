/*
Copyright (c) 2019, Lawrence Livermore National Security, LLC
and DESUL project contributors. See the COPYRIGHT file for details.
Source: https://github.com/desul/desul

SPDX-License-Identifier: (BSD-3-Clause)
*/

#ifndef DESUL_ATOMICS_GENERIC_HPP_
#define DESUL_ATOMICS_GENERIC_HPP_

#include <type_traits>
#if defined(__GNUC__) && (!defined(__clang__))
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif
#include "desul/atomics/Common.hpp"
#include "desul/atomics/Compare_Exchange.hpp"
#include "desul/atomics/Lock_Array.hpp"
#include "desul/atomics/Macros.hpp"
// Combination operands to be used in an Compare and Exchange based atomic
// operation
namespace desul {
namespace Impl {

template <class Scalar1, class Scalar2>
struct MaxOper {
  DESUL_FORCEINLINE_FUNCTION
  static Scalar1 apply(const Scalar1& val1, const Scalar2& val2) {
    return (val1 > val2 ? val1 : val2);
  }
  DESUL_FORCEINLINE_FUNCTION
  static constexpr bool check_early_exit(Scalar1 const& val1, Scalar2 const& val2) {
    return val1 > val2;
  }
};

template <class Scalar1, class Scalar2>
struct MinOper {
  DESUL_FORCEINLINE_FUNCTION
  static Scalar1 apply(const Scalar1& val1, const Scalar2& val2) {
    return (val1 < val2 ? val1 : val2);
  }
  DESUL_FORCEINLINE_FUNCTION
  static constexpr bool check_early_exit(Scalar1 const& val1, Scalar2 const& val2) {
    return val1 < val2;
  }
};

template <typename Op, typename Scalar1, typename Scalar2, typename = bool>
struct may_exit_early : std::false_type {};

// This exit early optimization causes weird compiler errors with MSVC 2019
#ifndef DESUL_HAVE_MSVC_ATOMICS
template <typename Op, typename Scalar1, typename Scalar2>
struct may_exit_early<Op,
                      Scalar1,
                      Scalar2,
                      decltype(Op::check_early_exit(std::declval<Scalar1 const&>(),
                                                    std::declval<Scalar2 const&>()))>
    : std::true_type {};
#endif

template <typename Op, typename Scalar1, typename Scalar2>
constexpr DESUL_FUNCTION
    typename std::enable_if<may_exit_early<Op, Scalar1, Scalar2>::value, bool>::type
    check_early_exit(Op const&, Scalar1 const& val1, Scalar2 const& val2) {
  return Op::check_early_exit(val1, val2);
}

template <typename Op, typename Scalar1, typename Scalar2>
constexpr DESUL_FUNCTION
    typename std::enable_if<!may_exit_early<Op, Scalar1, Scalar2>::value, bool>::type
    check_early_exit(Op const&, Scalar1 const&, Scalar2 const&) {
  return false;
}

template <class Scalar1, class Scalar2>
struct AddOper {
  DESUL_FORCEINLINE_FUNCTION
  static Scalar1 apply(const Scalar1& val1, const Scalar2& val2) { return val1 + val2; }
};

template <class Scalar1, class Scalar2>
struct SubOper {
  DESUL_FORCEINLINE_FUNCTION
  static Scalar1 apply(const Scalar1& val1, const Scalar2& val2) { return val1 - val2; }
};

template <class Scalar1, class Scalar2>
struct MulOper {
  DESUL_FORCEINLINE_FUNCTION
  static Scalar1 apply(const Scalar1& val1, const Scalar2& val2) { return val1 * val2; }
};

template <class Scalar1, class Scalar2>
struct DivOper {
  DESUL_FORCEINLINE_FUNCTION
  static Scalar1 apply(const Scalar1& val1, const Scalar2& val2) { return val1 / val2; }
};

template <class Scalar1, class Scalar2>
struct ModOper {
  DESUL_FORCEINLINE_FUNCTION
  static Scalar1 apply(const Scalar1& val1, const Scalar2& val2) { return val1 % val2; }
};

template <class Scalar1, class Scalar2>
struct AndOper {
  DESUL_FORCEINLINE_FUNCTION
  static Scalar1 apply(const Scalar1& val1, const Scalar2& val2) { return val1 & val2; }
};

template <class Scalar1, class Scalar2>
struct OrOper {
  DESUL_FORCEINLINE_FUNCTION
  static Scalar1 apply(const Scalar1& val1, const Scalar2& val2) { return val1 | val2; }
};

template <class Scalar1, class Scalar2>
struct XorOper {
  DESUL_FORCEINLINE_FUNCTION
  static Scalar1 apply(const Scalar1& val1, const Scalar2& val2) { return val1 ^ val2; }
};

template <class Scalar1, class Scalar2>
struct NandOper {
  DESUL_FORCEINLINE_FUNCTION
  static Scalar1 apply(const Scalar1& val1, const Scalar2& val2) {
    return ~(val1 & val2);
  }
};

template <class Scalar1, class Scalar2>
struct LShiftOper {
  DESUL_FORCEINLINE_FUNCTION
  static Scalar1 apply(const Scalar1& val1, const Scalar2& val2) {
    return val1 << val2;
  }
};

template <class Scalar1, class Scalar2>
struct RShiftOper {
  DESUL_FORCEINLINE_FUNCTION
  static Scalar1 apply(const Scalar1& val1, const Scalar2& val2) {
    return val1 >> val2;
  }
};

template <class Scalar1, class Scalar2>
struct IncModOper {
  DESUL_FORCEINLINE_FUNCTION
  static Scalar1 apply(const Scalar1& val1, const Scalar2& val2) {
    return ((val1 >= val2) ? Scalar1(0) : val1 + Scalar1(1));
  }
};

template <class Scalar1, class Scalar2>
struct DecModOper {
  DESUL_FORCEINLINE_FUNCTION
  static Scalar1 apply(const Scalar1& val1, const Scalar2& val2) {
    return (((val1 == Scalar1(0)) | (val1 > val2)) ? val2 : (val1 - Scalar1(1)));
  }
};

template <class Scalar1, class Scalar2>
struct StoreOper {
  DESUL_FORCEINLINE_FUNCTION
  static Scalar1 apply(const Scalar1&, const Scalar2& val2) { return val2; }
};

template <class Scalar1, class Scalar2>
struct LoadOper {
  DESUL_FORCEINLINE_FUNCTION
  static Scalar1 apply(const Scalar1& val1, const Scalar2&) { return val1; }
};

template <class Oper,
          typename T,
          class MemoryOrder,
          class MemoryScope,
          // equivalent to:
          //   requires atomic_always_lock_free(sizeof(T))
          std::enable_if_t<atomic_always_lock_free(sizeof(T)), int> = 0>
DESUL_INLINE_FUNCTION T atomic_fetch_oper(const Oper& op,
                                          T* const dest,
                                          dont_deduce_this_parameter_t<const T> val,
                                          MemoryOrder order,
                                          MemoryScope scope) {
  using cas_t = typename atomic_compare_exchange_type<sizeof(T)>::type;
  cas_t oldval = reinterpret_cast<cas_t&>(*dest);
  cas_t assume = oldval;

  do {
    if (Impl::check_early_exit(op, reinterpret_cast<T&>(oldval), val))
      return reinterpret_cast<T&>(oldval);
    assume = oldval;
    T newval = op.apply(reinterpret_cast<T&>(assume), val);
    oldval = desul::atomic_compare_exchange(reinterpret_cast<cas_t*>(dest),
                                            assume,
                                            reinterpret_cast<cas_t&>(newval),
                                            order,
                                            scope);
  } while (assume != oldval);

  return reinterpret_cast<T&>(oldval);
}

template <class Oper,
          typename T,
          class MemoryOrder,
          class MemoryScope,
          // equivalent to:
          //   requires atomic_always_lock_free(sizeof(T))
          std::enable_if_t<atomic_always_lock_free(sizeof(T)), int> = 0>
DESUL_INLINE_FUNCTION T atomic_oper_fetch(const Oper& op,
                                          T* const dest,
                                          dont_deduce_this_parameter_t<const T> val,
                                          MemoryOrder order,
                                          MemoryScope scope) {
  using cas_t = typename atomic_compare_exchange_type<sizeof(T)>::type;
  cas_t oldval = reinterpret_cast<cas_t&>(*dest);
  T newval = val;
  cas_t assume = oldval;
  do {
    if (Impl::check_early_exit(op, reinterpret_cast<T&>(oldval), val))
      return reinterpret_cast<T&>(oldval);
    assume = oldval;
    newval = op.apply(reinterpret_cast<T&>(assume), val);
    oldval = desul::atomic_compare_exchange(reinterpret_cast<cas_t*>(dest),
                                            assume,
                                            reinterpret_cast<cas_t&>(newval),
                                            order,
                                            scope);
  } while (assume != oldval);

  return newval;
}

template <class Oper,
          typename T,
          class MemoryOrder,
          class MemoryScope,
          // equivalent to:
          //   requires !atomic_always_lock_free(sizeof(T))
          std::enable_if_t<!atomic_always_lock_free(sizeof(T)), int> = 0>
DESUL_INLINE_FUNCTION T atomic_fetch_oper(const Oper& op,
                                          T* const dest,
                                          dont_deduce_this_parameter_t<const T> val,
                                          MemoryOrder /*order*/,
                                          MemoryScope scope) {
#if defined(DESUL_HAVE_FORWARD_PROGRESS)
  // Acquire a lock for the address
  while (!Impl::lock_address((void*)dest, scope)) {
  }

  atomic_thread_fence(MemoryOrderAcquire(), scope);
  T return_val = *dest;
  *dest = op.apply(return_val, val);
  atomic_thread_fence(MemoryOrderRelease(), scope);
  Impl::unlock_address((void*)dest, scope);
  return return_val;
#elif defined(DESUL_HAVE_GPU_LIKE_PROGRESS)
  // This is a way to avoid dead lock in a warp or wave front
  T return_val;
  int done = 0;
#ifdef __HIPCC__
  unsigned long long int active = DESUL_IMPL_BALLOT_MASK(1);
  unsigned long long int done_active = 0;
  while (active != done_active) {
    if (!done) {
      if (Impl::lock_address_hip((void*)dest, scope)) {
        atomic_thread_fence(MemoryOrderAcquire(), scope);
        return_val = *dest;
        *dest = op.apply(return_val, val);
        atomic_thread_fence(MemoryOrderRelease(), scope);
        Impl::unlock_address_hip((void*)dest, scope);
        done = 1;
      }
    }
    done_active = DESUL_IMPL_BALLOT_MASK(done);
  }
  return return_val;
// FIXME_SYCL not implemented
#elif defined(__SYCL_DEVICE_ONLY__)
  (void)op;
  (void)dest;
  (void)scope;
  (void)return_val;
  (void)done;

  assert(false);
  return val;
#else
  unsigned int mask = DESUL_IMPL_ACTIVEMASK;
  unsigned int active = DESUL_IMPL_BALLOT_MASK(mask, 1);
  unsigned int done_active = 0;
  while (active != done_active) {
    if (!done) {
      if (Impl::lock_address_cuda((void*)dest, scope)) {
        atomic_thread_fence(MemoryOrderAcquire(), scope);
        return_val = *dest;
        *dest = op.apply(return_val, val);
        atomic_thread_fence(MemoryOrderRelease(), scope);
        Impl::unlock_address_cuda((void*)dest, scope);
        done = 1;
      }
    }
    done_active = DESUL_IMPL_BALLOT_MASK(mask, done);
  }
  return return_val;
#endif
#else
  static_assert(false, "Unimplemented lock based atomic\n");
  return val;
#endif
}

template <class Oper,
          typename T,
          class MemoryOrder,
          class MemoryScope,
          // equivalent to:
          //   requires !atomic_always_lock_free(sizeof(T))
          std::enable_if_t<!atomic_always_lock_free(sizeof(T)), int> = 0>
DESUL_INLINE_FUNCTION T atomic_oper_fetch(const Oper& op,
                                          T* const dest,
                                          dont_deduce_this_parameter_t<const T> val,
                                          MemoryOrder /*order*/,
                                          MemoryScope scope) {
#if defined(DESUL_HAVE_FORWARD_PROGRESS)
  // Acquire a lock for the address
  while (!Impl::lock_address((void*)dest, scope)) {
  }

  atomic_thread_fence(MemoryOrderAcquire(), scope);
  T return_val = op.apply(*dest, val);
  *dest = return_val;
  atomic_thread_fence(MemoryOrderRelease(), scope);
  Impl::unlock_address((void*)dest, scope);
  return return_val;
#elif defined(DESUL_HAVE_GPU_LIKE_PROGRESS)
  // This is a way to avoid dead lock in a warp or wave front
  T return_val;
  int done = 0;
#ifdef __HIPCC__
  unsigned long long int active = DESUL_IMPL_BALLOT_MASK(1);
  unsigned long long int done_active = 0;
  while (active != done_active) {
    if (!done) {
      if (Impl::lock_address_hip((void*)dest, scope)) {
        atomic_thread_fence(MemoryOrderAcquire(), scope);
        return_val = op.apply(*dest, val);
        *dest = return_val;
        atomic_thread_fence(MemoryOrderRelease(), scope);
        Impl::unlock_address_hip((void*)dest, scope);
        done = 1;
      }
    }
    done_active = DESUL_IMPL_BALLOT_MASK(done);
  }
  return return_val;
  // FIXME_SYCL not implemented
#elif defined(__SYCL_DEVICE_ONLY__)
  (void)op;
  (void)dest;
  (void)scope;
  (void)done;

  assert(false);
  return val;
#else
  unsigned int mask = DESUL_IMPL_ACTIVEMASK;
  unsigned int active = DESUL_IMPL_BALLOT_MASK(mask, 1);
  unsigned int done_active = 0;
  while (active != done_active) {
    if (!done) {
      if (Impl::lock_address_cuda((void*)dest, scope)) {
        atomic_thread_fence(MemoryOrderAcquire(), scope);
        return_val = op.apply(*dest, val);
        *dest = return_val;
        atomic_thread_fence(MemoryOrderRelease(), scope);
        Impl::unlock_address_cuda((void*)dest, scope);
        done = 1;
      }
    }
    done_active = DESUL_IMPL_BALLOT_MASK(mask, done);
  }
  return return_val;
#endif
#else
  static_assert(false, "Unimplemented lock based atomic\n");
  return val;
#endif
}

template <class Oper, typename T, class MemoryOrder>
DESUL_INLINE_FUNCTION T atomic_fetch_oper(const Oper& op,
                                          T* const dest,
                                          dont_deduce_this_parameter_t<const T> val,
                                          MemoryOrder /*order*/,
                                          MemoryScopeCaller /*scope*/) {
  T oldval = *dest;
  *dest = op.apply(oldval, val);
  return oldval;
}

template <class Oper, typename T, class MemoryOrder>
DESUL_INLINE_FUNCTION T atomic_oper_fetch(const Oper& op,
                                          T* const dest,
                                          dont_deduce_this_parameter_t<const T> val,
                                          MemoryOrder /*order*/,
                                          MemoryScopeCaller /*scope*/) {
  T oldval = *dest;
  T newval = op.apply(oldval, val);
  *dest = newval;
  return newval;
}

}  // namespace Impl
}  // namespace desul

namespace desul {

// Fetch_Oper atomics: return value before operation
template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_fetch_add(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  return Impl::atomic_fetch_oper(Impl::AddOper<T, const T>(), dest, val, order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_fetch_sub(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  return Impl::atomic_fetch_oper(Impl::SubOper<T, const T>(), dest, val, order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_fetch_max(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  return Impl::atomic_fetch_oper(Impl::MaxOper<T, const T>(), dest, val, order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_fetch_min(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  return Impl::atomic_fetch_oper(Impl::MinOper<T, const T>(), dest, val, order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_fetch_mul(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  return Impl::atomic_fetch_oper(Impl::MulOper<T, const T>(), dest, val, order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_fetch_div(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  return Impl::atomic_fetch_oper(Impl::DivOper<T, const T>(), dest, val, order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_fetch_mod(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  return Impl::atomic_fetch_oper(Impl::ModOper<T, const T>(), dest, val, order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_fetch_and(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  return Impl::atomic_fetch_oper(Impl::AndOper<T, const T>(), dest, val, order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_fetch_or(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  return Impl::atomic_fetch_oper(Impl::OrOper<T, const T>(), dest, val, order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_fetch_xor(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  return Impl::atomic_fetch_oper(Impl::XorOper<T, const T>(), dest, val, order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_fetch_nand(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  return Impl::atomic_fetch_oper(Impl::NandOper<T, const T>(), dest, val, order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T atomic_fetch_lshift(T* const dest,
                                            const unsigned int val,
                                            MemoryOrder order,
                                            MemoryScope scope) {
  return Impl::atomic_fetch_oper(
      Impl::LShiftOper<T, const unsigned int>(), dest, val, order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T atomic_fetch_rshift(T* const dest,
                                            const unsigned int val,
                                            MemoryOrder order,
                                            MemoryScope scope) {
  return Impl::atomic_fetch_oper(
      Impl::RShiftOper<T, const unsigned int>(), dest, val, order, scope);
}

// Oper Fetch atomics: return value after operation
template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_add_fetch(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  return Impl::atomic_oper_fetch(Impl::AddOper<T, const T>(), dest, val, order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_sub_fetch(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  return Impl::atomic_oper_fetch(Impl::SubOper<T, const T>(), dest, val, order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_max_fetch(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  return Impl::atomic_oper_fetch(Impl::MaxOper<T, const T>(), dest, val, order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_min_fetch(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  return Impl::atomic_oper_fetch(Impl::MinOper<T, const T>(), dest, val, order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_mul_fetch(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  return Impl::atomic_oper_fetch(Impl::MulOper<T, const T>(), dest, val, order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_div_fetch(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  return Impl::atomic_oper_fetch(Impl::DivOper<T, const T>(), dest, val, order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_mod_fetch(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  return Impl::atomic_oper_fetch(Impl::ModOper<T, const T>(), dest, val, order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_and_fetch(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  return Impl::atomic_oper_fetch(Impl::AndOper<T, const T>(), dest, val, order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_or_fetch(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  return Impl::atomic_oper_fetch(Impl::OrOper<T, const T>(), dest, val, order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_xor_fetch(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  return Impl::atomic_oper_fetch(Impl::XorOper<T, const T>(), dest, val, order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_nand_fetch(T* const dest, const T val, MemoryOrder order, MemoryScope scope) {
  return Impl::atomic_oper_fetch(Impl::NandOper<T, const T>(), dest, val, order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T atomic_lshift_fetch(T* const dest,
                                            const unsigned int val,
                                            MemoryOrder order,
                                            MemoryScope scope) {
  return Impl::atomic_oper_fetch(
      Impl::LShiftOper<T, const unsigned int>(), dest, val, order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T atomic_rshift_fetch(T* const dest,
                                            const unsigned int val,
                                            MemoryOrder order,
                                            MemoryScope scope) {
  return Impl::atomic_oper_fetch(
      Impl::RShiftOper<T, const unsigned int>(), dest, val, order, scope);
}

// Other atomics

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T atomic_load(const T* const dest,
                                    MemoryOrder order,
                                    MemoryScope scope) {
  return Impl::atomic_fetch_oper(
      Impl::LoadOper<T, const T>(), const_cast<T*>(dest), T(), order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION void atomic_store(T* const dest,
                                        const T val,
                                        MemoryOrder order,
                                        MemoryScope scope) {
  (void)Impl::atomic_fetch_oper(Impl::StoreOper<T, const T>(), dest, val, order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION void atomic_add(T* const dest,
                                      const T val,
                                      MemoryOrder order,
                                      MemoryScope scope) {
  (void)atomic_fetch_add(dest, val, order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION void atomic_sub(T* const dest,
                                      const T val,
                                      MemoryOrder order,
                                      MemoryScope scope) {
  (void)atomic_fetch_sub(dest, val, order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION void atomic_mul(T* const dest,
                                      const T val,
                                      MemoryOrder order,
                                      MemoryScope scope) {
  (void)atomic_fetch_mul(dest, val, order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION void atomic_div(T* const dest,
                                      const T val,
                                      MemoryOrder order,
                                      MemoryScope scope) {
  (void)atomic_fetch_div(dest, val, order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION void atomic_min(T* const dest,
                                      const T val,
                                      MemoryOrder order,
                                      MemoryScope scope) {
  (void)atomic_fetch_min(dest, val, order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION void atomic_max(T* const dest,
                                      const T val,
                                      MemoryOrder order,
                                      MemoryScope scope) {
  (void)atomic_fetch_max(dest, val, order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T atomic_inc_fetch(T* const dest,
                                         MemoryOrder order,
                                         MemoryScope scope) {
  return atomic_add_fetch(dest, T(1), order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T atomic_dec_fetch(T* const dest,
                                         MemoryOrder order,
                                         MemoryScope scope) {
  return atomic_sub_fetch(dest, T(1), order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T atomic_fetch_inc(T* const dest,
                                         MemoryOrder order,
                                         MemoryScope scope) {
  return atomic_fetch_add(dest, T(1), order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_fetch_inc_mod(T* const dest, T val, MemoryOrder order, MemoryScope scope) {
  static_assert(std::is_unsigned<T>::value,
                "Signed types not supported by atomic_fetch_inc_mod.");
  return Impl::atomic_fetch_oper(
      Impl::IncModOper<T, const T>(), dest, val, order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T atomic_fetch_dec(T* const dest,
                                         MemoryOrder order,
                                         MemoryScope scope) {
  return atomic_fetch_sub(dest, T(1), order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION T
atomic_fetch_dec_mod(T* const dest, T val, MemoryOrder order, MemoryScope scope) {
  static_assert(std::is_unsigned<T>::value,
                "Signed types not supported by atomic_fetch_dec_mod.");
  return Impl::atomic_fetch_oper(
      Impl::DecModOper<T, const T>(), dest, val, order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION void atomic_inc(T* const dest,
                                      MemoryOrder order,
                                      MemoryScope scope) {
  return atomic_add(dest, T(1), order, scope);
}

template <typename T, class MemoryOrder, class MemoryScope>
DESUL_INLINE_FUNCTION void atomic_dec(T* const dest,
                                      MemoryOrder order,
                                      MemoryScope scope) {
  return atomic_sub(dest, T(1), order, scope);
}

// FIXME
template <typename T,
          class SuccessMemoryOrder,
          class FailureMemoryOrder,
          class MemoryScope>
DESUL_INLINE_FUNCTION bool atomic_compare_exchange_strong(
    T* const dest,
    T& expected,
    T desired,
    SuccessMemoryOrder success,
    FailureMemoryOrder /*failure*/,
    MemoryScope scope) {
  T const old = atomic_compare_exchange(dest, expected, desired, success, scope);
  if (old != expected) {
    expected = old;
    return false;
  } else {
    return true;
  }
}

template <typename T,
          class SuccessMemoryOrder,
          class FailureMemoryOrder,
          class MemoryScope>
DESUL_INLINE_FUNCTION bool atomic_compare_exchange_weak(T* const dest,
                                                        T& expected,
                                                        T desired,
                                                        SuccessMemoryOrder success,
                                                        FailureMemoryOrder failure,
                                                        MemoryScope scope) {
  return atomic_compare_exchange_strong(
      dest, expected, desired, success, failure, scope);
}

}  // namespace desul

#include <desul/atomics/CUDA.hpp>
#include <desul/atomics/GCC.hpp>
#include <desul/atomics/HIP.hpp>
#include <desul/atomics/OpenMP.hpp>
#include <desul/atomics/SYCL.hpp>
#if defined(__GNUC__) && (!defined(__clang__))
#pragma GCC diagnostic pop
#endif
#endif
