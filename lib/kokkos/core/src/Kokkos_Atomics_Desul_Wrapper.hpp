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

#include <impl/Kokkos_Utilities.hpp>  // identity_type
#include <impl/Kokkos_Volatile_Load.hpp>

namespace Kokkos {

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
KOKKOS_DEPRECATED inline const char* atomic_query_version() {
  return "KOKKOS_DESUL_ATOMICS";
}
#endif

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
#if defined(KOKKOS_COMPILER_GNU) && !defined(__PGIC__) && \
    !defined(__CUDA_ARCH__)

#define KOKKOS_NONTEMPORAL_PREFETCH_LOAD(addr) __builtin_prefetch(addr, 0, 0)
#define KOKKOS_NONTEMPORAL_PREFETCH_STORE(addr) __builtin_prefetch(addr, 1, 0)

#else

#define KOKKOS_NONTEMPORAL_PREFETCH_LOAD(addr) ((void)0)
#define KOKKOS_NONTEMPORAL_PREFETCH_STORE(addr) ((void)0)

#endif
#endif
// ============================================================

#ifdef KOKKOS_ENABLE_ATOMICS_BYPASS
#define KOKKOS_DESUL_MEM_SCOPE desul::MemoryScopeCaller()
#else
#define KOKKOS_DESUL_MEM_SCOPE desul::MemoryScopeDevice()
#endif

namespace Impl {
template <class T>
using not_deduced_atomic_t =
    std::add_const_t<std::remove_volatile_t<type_identity_t<T>>>;

template <class T, class R>
using enable_if_atomic_t =
    std::enable_if_t<!std::is_reference_v<T> && !std::is_const_v<T>,
                     std::remove_volatile_t<R>>;
}  // namespace Impl

// clang-format off

// fences
KOKKOS_INLINE_FUNCTION void memory_fence() { desul::atomic_thread_fence(desul::MemoryOrderSeqCst(),  KOKKOS_DESUL_MEM_SCOPE); }
KOKKOS_INLINE_FUNCTION void load_fence()   { desul::atomic_thread_fence(desul::MemoryOrderAcquire(), KOKKOS_DESUL_MEM_SCOPE); }
KOKKOS_INLINE_FUNCTION void store_fence()  { desul::atomic_thread_fence(desul::MemoryOrderRelease(), KOKKOS_DESUL_MEM_SCOPE); }

// load/store
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T,    T> atomic_load (T const* ptr)                              { return desul::atomic_load (const_cast<std::remove_volatile_t<T>*>(ptr),      desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, void> atomic_store(T* ptr, Impl::not_deduced_atomic_t<T> val) { return desul::atomic_store(const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
template<class T> KOKKOS_DEPRECATED_WITH_COMMENT("Use atomic_store() instead!") KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, void> atomic_assign(T* ptr, Impl::not_deduced_atomic_t<T> val) { atomic_store(ptr, val); }
#endif

// atomic_fetch_op
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, T> atomic_fetch_add(T* ptr, Impl::not_deduced_atomic_t<T> val) { return desul::atomic_fetch_add(const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, T> atomic_fetch_sub(T* ptr, Impl::not_deduced_atomic_t<T> val) { return desul::atomic_fetch_sub(const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, T> atomic_fetch_max(T* ptr, Impl::not_deduced_atomic_t<T> val) { return desul::atomic_fetch_max(const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, T> atomic_fetch_min(T* ptr, Impl::not_deduced_atomic_t<T> val) { return desul::atomic_fetch_min(const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, T> atomic_fetch_mul(T* ptr, Impl::not_deduced_atomic_t<T> val) { return desul::atomic_fetch_mul(const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, T> atomic_fetch_div(T* ptr, Impl::not_deduced_atomic_t<T> val) { return desul::atomic_fetch_div(const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, T> atomic_fetch_mod(T* ptr, Impl::not_deduced_atomic_t<T> val) { return desul::atomic_fetch_mod(const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, T> atomic_fetch_and(T* ptr, Impl::not_deduced_atomic_t<T> val) { return desul::atomic_fetch_and(const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, T> atomic_fetch_or (T* ptr, Impl::not_deduced_atomic_t<T> val) { return desul::atomic_fetch_or (const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, T> atomic_fetch_xor(T* ptr, Impl::not_deduced_atomic_t<T> val) { return desul::atomic_fetch_xor(const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, T> atomic_fetch_nand(T* ptr, Impl::not_deduced_atomic_t<T> val) { return desul::atomic_fetch_nand(const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, T> atomic_fetch_lshift(T* ptr, Impl::not_deduced_atomic_t<T> val) { return desul::atomic_fetch_lshift(const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, T> atomic_fetch_rshift(T* ptr, Impl::not_deduced_atomic_t<T> val) { return desul::atomic_fetch_rshift(const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, T> atomic_fetch_inc(T* ptr) { return desul::atomic_fetch_inc(const_cast<std::remove_volatile_t<T>*>(ptr), desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, T> atomic_fetch_dec(T* ptr) { return desul::atomic_fetch_dec(const_cast<std::remove_volatile_t<T>*>(ptr), desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

// atomic_op_fetch
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, T> atomic_add_fetch(T* ptr, Impl::not_deduced_atomic_t<T> val) { return desul::atomic_add_fetch(const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, T> atomic_sub_fetch(T* ptr, Impl::not_deduced_atomic_t<T> val) { return desul::atomic_sub_fetch(const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, T> atomic_max_fetch(T* ptr, Impl::not_deduced_atomic_t<T> val) { return desul::atomic_max_fetch(const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, T> atomic_min_fetch(T* ptr, Impl::not_deduced_atomic_t<T> val) { return desul::atomic_min_fetch(const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, T> atomic_mul_fetch(T* ptr, Impl::not_deduced_atomic_t<T> val) { return desul::atomic_mul_fetch(const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, T> atomic_div_fetch(T* ptr, Impl::not_deduced_atomic_t<T> val) { return desul::atomic_div_fetch(const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, T> atomic_mod_fetch(T* ptr, Impl::not_deduced_atomic_t<T> val) { return desul::atomic_mod_fetch(const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, T> atomic_and_fetch(T* ptr, Impl::not_deduced_atomic_t<T> val) { return desul::atomic_and_fetch(const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, T> atomic_or_fetch (T* ptr, Impl::not_deduced_atomic_t<T> val) { return desul::atomic_or_fetch (const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, T> atomic_xor_fetch(T* ptr, Impl::not_deduced_atomic_t<T> val) { return desul::atomic_xor_fetch(const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, T> atomic_nand_fetch(T* ptr, Impl::not_deduced_atomic_t<T> val) { return desul::atomic_nand_fetch(const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, T> atomic_lshift_fetch(T* ptr, Impl::not_deduced_atomic_t<T> val) { return desul::atomic_lshift_fetch(const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, T> atomic_rshift_fetch(T* ptr, Impl::not_deduced_atomic_t<T> val) { return desul::atomic_rshift_fetch(const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, T> atomic_inc_fetch(T* ptr) { return desul::atomic_inc_fetch(const_cast<std::remove_volatile_t<T>*>(ptr), desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, T> atomic_dec_fetch(T* ptr) { return desul::atomic_dec_fetch(const_cast<std::remove_volatile_t<T>*>(ptr), desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }

// atomic_op
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, void> atomic_add(T* ptr, Impl::not_deduced_atomic_t<T> val) { desul::atomic_add(const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, void> atomic_sub(T* ptr, Impl::not_deduced_atomic_t<T> val) { desul::atomic_sub(const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, void> atomic_max(T* ptr, Impl::not_deduced_atomic_t<T> val) { desul::atomic_max(const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, void> atomic_min(T* ptr, Impl::not_deduced_atomic_t<T> val) { desul::atomic_min(const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, void> atomic_mul(T* ptr, Impl::not_deduced_atomic_t<T> val) { desul::atomic_mul(const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, void> atomic_div(T* ptr, Impl::not_deduced_atomic_t<T> val) { desul::atomic_div(const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, void> atomic_mod(T* ptr, Impl::not_deduced_atomic_t<T> val) { (void)desul::atomic_fetch_mod(const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, void> atomic_and(T* ptr, Impl::not_deduced_atomic_t<T> val) { (void)desul::atomic_fetch_and(const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, void> atomic_or (T* ptr, Impl::not_deduced_atomic_t<T> val) { (void)desul::atomic_fetch_or (const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, void> atomic_xor(T* ptr, Impl::not_deduced_atomic_t<T> val) { (void)desul::atomic_fetch_xor(const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T>  KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, void> atomic_nand(T* ptr, Impl::not_deduced_atomic_t<T> val) { (void)desul::atomic_nand_fetch(const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, void> atomic_lshift(T* ptr, Impl::not_deduced_atomic_t<T> val) { (void)desul::atomic_fetch_lshift(const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, void> atomic_rshift(T* ptr, Impl::not_deduced_atomic_t<T> val) { (void)desul::atomic_fetch_rshift(const_cast<std::remove_volatile_t<T>*>(ptr), val, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, void> atomic_inc(T* ptr) { desul::atomic_inc(const_cast<std::remove_volatile_t<T>*>(ptr), desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, void> atomic_dec(T* ptr) { desul::atomic_dec(const_cast<std::remove_volatile_t<T>*>(ptr), desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
template<class T> KOKKOS_DEPRECATED_WITH_COMMENT("Use atomic_inc() instead!") KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, void> atomic_increment(T* ptr) { atomic_inc(ptr); }
template<class T> KOKKOS_DEPRECATED_WITH_COMMENT("Use atomic_dec() instead!") KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, void> atomic_decrement(T* ptr) { atomic_dec(ptr); }
#endif

// exchange
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, T> atomic_exchange        (T* ptr, Impl::not_deduced_atomic_t<T> val)                                             { return desul::atomic_exchange        (const_cast<std::remove_volatile_t<T>*>(ptr), val,               desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
template<class T> KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, T> atomic_compare_exchange(T* ptr, Impl::not_deduced_atomic_t<T> expected, Impl::not_deduced_atomic_t<T> desired) { return desul::atomic_compare_exchange(const_cast<std::remove_volatile_t<T>*>(ptr), expected, desired, desul::MemoryOrderRelaxed(), KOKKOS_DESUL_MEM_SCOPE); }
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
template<class T> KOKKOS_DEPRECATED_WITH_COMMENT("Use atomic_compare_exchange() instead!") KOKKOS_FUNCTION Impl::enable_if_atomic_t<T, bool> atomic_compare_exchange_strong(T* ptr, Impl::not_deduced_atomic_t<T> expected, Impl::not_deduced_atomic_t<T> desired) { return expected == atomic_compare_exchange(ptr, expected, desired); }
#endif

// clang-format on
}  // namespace Kokkos

namespace Kokkos::Impl {

template <class T, class MemOrderSuccess, class MemOrderFailure>
KOKKOS_FUNCTION bool atomic_compare_exchange_strong(T* const dest, T& expected,
                                                    const T desired,
                                                    MemOrderSuccess succ,
                                                    MemOrderFailure fail) {
  return desul::atomic_compare_exchange_strong(dest, expected, desired, succ,
                                               fail, KOKKOS_DESUL_MEM_SCOPE);
}

template <class T, class MemoryOrder>
KOKKOS_FUNCTION T atomic_load(const T* const src, MemoryOrder order) {
  return desul::atomic_load(src, order, KOKKOS_DESUL_MEM_SCOPE);
}

template <class T, class MemoryOrder>
KOKKOS_FUNCTION void atomic_store(T* const src, const T val,
                                  MemoryOrder order) {
  return desul::atomic_store(src, val, order, KOKKOS_DESUL_MEM_SCOPE);
}

}  // namespace Kokkos::Impl

#undef KOKKOS_DESUL_MEM_SCOPE

#endif
