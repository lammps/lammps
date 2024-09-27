//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#include <Kokkos_Macros.hpp>
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#endif
#ifndef KOKKOS_DESUL_ATOMICS_WRAPPER_HPP_
#define KOKKOS_DESUL_ATOMICS_WRAPPER_HPP_
#include <Kokkos_Macros.hpp>
#include <desul/atomics.hpp>

#include <impl/Kokkos_Volatile_Load.hpp>

// clang-format off
namespace Kokkos {

// FIXME: These functions don't have any use/test in unit tests ...
// ==========================================================
inline const char* atomic_query_version() { return "KOKKOS_DESUL_ATOMICS"; }

#if defined(KOKKOS_COMPILER_GNU) && !defined(__PGIC__) && \
    !defined(__CUDA_ARCH__)

#define KOKKOS_NONTEMPORAL_PREFETCH_LOAD(addr) __builtin_prefetch(addr, 0, 0)
#define KOKKOS_NONTEMPORAL_PREFETCH_STORE(addr) __builtin_prefetch(addr, 1, 0)

#else

#define KOKKOS_NONTEMPORAL_PREFETCH_LOAD(addr) ((void)0)
#define KOKKOS_NONTEMPORAL_PREFETCH_STORE(addr) ((void)0)

#endif
// ============================================================

#ifdef KOKKOS_ENABLE_ATOMICS_BYPASS
#define KOKKOS_DESUL_MEM_SCOPE desul::MemoryScopeCaller()
#else
#define KOKKOS_DESUL_MEM_SCOPE desul::MemoryScopeDevice()
#endif

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_load(T* const dest) { return desul::atomic_load(dest, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
void atomic_store(T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_store(dest, val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
void atomic_assign(T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { atomic_store(dest,val); }

KOKKOS_INLINE_FUNCTION
void memory_fence() {
  desul::atomic_thread_fence(desul::MemoryOrderSeqCst(), KOKKOS_DESUL_MEM_SCOPE);
}

KOKKOS_INLINE_FUNCTION
void load_fence() { return desul::atomic_thread_fence(desul::MemoryOrderAcquire(), KOKKOS_DESUL_MEM_SCOPE); }

KOKKOS_INLINE_FUNCTION
void store_fence() { return desul::atomic_thread_fence(desul::MemoryOrderRelease(), KOKKOS_DESUL_MEM_SCOPE); }

// atomic_fetch_op
template<class T> KOKKOS_INLINE_FUNCTION
T atomic_fetch_add (T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_fetch_add (dest, val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_fetch_sub (T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_fetch_sub (dest, val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_fetch_max (T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_fetch_max (dest, val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_fetch_min (T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_fetch_min (dest, val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_fetch_mul (T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_fetch_mul (dest, val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_fetch_div (T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_fetch_div (dest, val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_fetch_mod (T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_fetch_mod (dest, val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_fetch_and (T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_fetch_and (dest, val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_fetch_or  (T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_fetch_or  (dest, val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_fetch_xor (T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_fetch_xor (dest, val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_fetch_nand(T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_fetch_nand(dest, val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_fetch_lshift(T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_fetch_lshift(dest, val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_fetch_rshift(T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_fetch_rshift(dest, val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_fetch_inc(T* const dest) { return desul::atomic_fetch_inc(dest, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_fetch_dec(T* const dest) { return desul::atomic_fetch_dec(dest, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }


// atomic_op_fetch
template<class T> KOKKOS_INLINE_FUNCTION
T atomic_add_fetch (T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_add_fetch (dest, val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_sub_fetch (T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_sub_fetch (dest, val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_max_fetch (T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_max_fetch (dest, val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_min_fetch (T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_min_fetch (dest, val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_mul_fetch (T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_mul_fetch (dest, val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_div_fetch (T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_div_fetch (dest, val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_mod_fetch (T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_mod_fetch (dest, val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_and_fetch (T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_and_fetch (dest, val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_or_fetch  (T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_or_fetch  (dest, val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_xor_fetch (T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_xor_fetch (dest, val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_nand_fetch(T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_nand_fetch(dest, val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_lshift_fetch(T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_lshift_fetch(dest, val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_rshift_fetch(T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_rshift_fetch(dest, val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_inc_fetch(T* const dest) { return desul::atomic_inc_fetch(dest, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_dec_fetch(T* const dest) { return desul::atomic_dec_fetch(dest, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }


// atomic_op
template<class T> KOKKOS_INLINE_FUNCTION
void atomic_add(T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_add (dest, val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
void atomic_sub(T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_sub (dest, val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
void atomic_mul(T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_mul (dest, val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
void atomic_div(T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_div (dest, val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
void atomic_min(T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_min (dest, val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
void atomic_max(T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_max (dest, val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

// FIXME: Desul doesn't have atomic_and yet so call fetch_and
template<class T> KOKKOS_INLINE_FUNCTION
void atomic_and(T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { (void) desul::atomic_fetch_and (dest, val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

// FIXME: Desul doesn't have atomic_or yet so call fetch_or
template<class T> KOKKOS_INLINE_FUNCTION
void atomic_or(T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val)  { (void) desul::atomic_fetch_or (dest, val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
void atomic_inc(T* const dest) { return desul::atomic_inc(dest, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
void atomic_dec(T* const dest) { return desul::atomic_dec(dest, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
void atomic_increment(T* const dest) { return desul::atomic_inc(dest, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
void atomic_decrement(T* const dest) { return desul::atomic_dec(dest, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

// Exchange

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_exchange(T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_exchange(dest, val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

template<class T> KOKKOS_INLINE_FUNCTION
bool atomic_compare_exchange_strong(T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> expected, desul::Impl::dont_deduce_this_parameter_t<const T> desired) {
  T expected_ref = expected;
  return desul::atomic_compare_exchange_strong(dest, expected_ref, desired,
                  desul::MemoryOrderRelaxed(), desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE);
}

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_compare_exchange(T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> compare, desul::Impl::dont_deduce_this_parameter_t<const T> desired) {
  return desul::atomic_compare_exchange(dest, compare, desired,
                  desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE);
}

namespace Impl {
  template<class T, class MemOrderSuccess, class MemOrderFailure> KOKKOS_INLINE_FUNCTION
  bool atomic_compare_exchange_strong(T* const dest, T& expected, const T desired, MemOrderSuccess succ, MemOrderFailure fail) {
    return desul::atomic_compare_exchange_strong(dest, expected, desired, succ, fail, KOKKOS_DESUL_MEM_SCOPE);
  }
  template<class T, class MemoryOrder>
  KOKKOS_INLINE_FUNCTION
  T atomic_load(const T* const src, MemoryOrder order) {
    return desul::atomic_load(src, order, KOKKOS_DESUL_MEM_SCOPE);
  }
  template<class T, class MemoryOrder>
  KOKKOS_INLINE_FUNCTION
  void atomic_store(T* const src, const T val, MemoryOrder order) {
    return desul::atomic_store(src, val, order, KOKKOS_DESUL_MEM_SCOPE);
  }
}  // namespace Impl

}  // namespace Kokkos

#undef KOKKOS_DESUL_MEM_SCOPE

// clang-format on
#endif
