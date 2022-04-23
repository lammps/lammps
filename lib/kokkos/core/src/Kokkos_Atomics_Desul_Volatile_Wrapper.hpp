#ifndef KOKKOS_DESUL_ATOMICS_VOLATILE_WRAPPER_HPP_
#define KOKKOS_DESUL_ATOMICS_VOLATILE_WRAPPER_HPP_
#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_IMPL_DESUL_ATOMICS
#include <desul/atomics.hpp>

// clang-format off
namespace Kokkos { 

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_load(volatile T* const dest) { return desul::atomic_load(const_cast<T*>(dest), desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

template<class T> KOKKOS_INLINE_FUNCTION
void atomic_store(volatile T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_store(const_cast<T*>(dest), val, desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

// atomic_fetch_op
template<class T> KOKKOS_INLINE_FUNCTION
T atomic_fetch_add (volatile T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_fetch_add (const_cast<T*>(dest), val, desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

#ifdef DESUL_IMPL_ATOMIC_CUDA_USE_DOUBLE_ATOMICADD
KOKKOS_INLINE_FUNCTION
double atomic_fetch_add(volatile double* const dest, double val) {
  #ifdef __CUDA_ARCH__
  return atomicAdd(const_cast<double*>(dest),val);
  #else
  return desul::atomic_fetch_add (const_cast<double*>(dest), val, desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice());
  #endif
};

KOKKOS_INLINE_FUNCTION
double atomic_fetch_sub(volatile double* const dest, double val) {
  #ifdef __CUDA_ARCH__
  return atomicAdd(const_cast<double*>(dest),-val);
  #else
  return desul::atomic_fetch_sub (const_cast<double*>(dest), val, desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice());
  #endif
};
#endif

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_fetch_sub (volatile T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_fetch_sub (const_cast<T*>(dest), val, desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_fetch_max (volatile T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_fetch_max (const_cast<T*>(dest), val, desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_fetch_min (volatile T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_fetch_min (const_cast<T*>(dest), val, desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_fetch_mul (volatile T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_fetch_mul (const_cast<T*>(dest), val, desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_fetch_div (volatile T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_fetch_div (const_cast<T*>(dest), val, desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_fetch_mod (volatile T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_fetch_mod (const_cast<T*>(dest), val, desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_fetch_and (volatile T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_fetch_and (const_cast<T*>(dest), val, desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_fetch_or  (volatile T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_fetch_or  (const_cast<T*>(dest), val, desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_fetch_xor (volatile T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_fetch_xor (const_cast<T*>(dest), val, desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_fetch_nand(volatile T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_fetch_nand(const_cast<T*>(dest), val, desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_fetch_lshift(volatile T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_fetch_lshift(const_cast<T*>(dest), val, desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_fetch_rshift(volatile T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_fetch_rshift(const_cast<T*>(dest), val, desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_fetch_inc(volatile T* const dest) { return desul::atomic_fetch_inc(const_cast<T*>(dest),desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_fetch_dec(volatile T* const dest) { return desul::atomic_fetch_dec(const_cast<T*>(dest),desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }


// atomic_op_fetch
template<class T> KOKKOS_INLINE_FUNCTION
T atomic_add_fetch (volatile T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_add_fetch (const_cast<T*>(dest), val, desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_sub_fetch (volatile T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_sub_fetch (const_cast<T*>(dest), val, desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_max_fetch (volatile T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_max_fetch (const_cast<T*>(dest), val, desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_min_fetch (volatile T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_min_fetch (const_cast<T*>(dest), val, desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_mul_fetch (volatile T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_mul_fetch (const_cast<T*>(dest), val, desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_div_fetch (volatile T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_div_fetch (const_cast<T*>(dest), val, desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_mod_fetch (volatile T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_mod_fetch (const_cast<T*>(dest), val, desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_and_fetch (volatile T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_and_fetch (const_cast<T*>(dest), val, desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_or_fetch  (volatile T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_or_fetch  (const_cast<T*>(dest), val, desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_xor_fetch (volatile T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_xor_fetch (const_cast<T*>(dest), val, desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_nand_fetch(volatile T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_nand_fetch(const_cast<T*>(dest), val, desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_lshift_fetch(volatile T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_lshift_fetch(const_cast<T*>(dest), val, desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_rshift_fetch(volatile T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_rshift_fetch(const_cast<T*>(dest), val, desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_inc_fetch(volatile T* const dest) { return desul::atomic_inc_fetch(const_cast<T*>(dest),desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_dec_fetch(volatile T* const dest) { return desul::atomic_dec_fetch(const_cast<T*>(dest),desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }


// atomic_op
template<class T> KOKKOS_INLINE_FUNCTION
void atomic_add(volatile T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_add (const_cast<T*>(dest), val, desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

template<class T> KOKKOS_INLINE_FUNCTION
void atomic_sub(volatile T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_sub (const_cast<T*>(dest), val, desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

template<class T> KOKKOS_INLINE_FUNCTION
void atomic_mul(volatile T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_mul (const_cast<T*>(dest), val, desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

template<class T> KOKKOS_INLINE_FUNCTION
void atomic_div(volatile T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_div (const_cast<T*>(dest), val, desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

template<class T> KOKKOS_INLINE_FUNCTION
void atomic_min(volatile T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_min (const_cast<T*>(dest), val, desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

template<class T> KOKKOS_INLINE_FUNCTION
void atomic_max(volatile T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_max (const_cast<T*>(dest), val, desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

// FIXME: Desul doesn't have atomic_and yet so call fetch_and
template<class T> KOKKOS_INLINE_FUNCTION
void atomic_and(volatile T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { (void) desul::atomic_fetch_and (const_cast<T*>(dest), val, desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

// FIXME: Desul doesn't have atomic_or yet so call fetch_or
template<class T> KOKKOS_INLINE_FUNCTION
void atomic_or (volatile T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { (void) desul::atomic_fetch_or  (const_cast<T*>(dest), val, desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

template<class T> KOKKOS_INLINE_FUNCTION
void atomic_inc(volatile T* const dest) { return desul::atomic_inc(const_cast<T*>(dest),desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

template<class T> KOKKOS_INLINE_FUNCTION
void atomic_dec(volatile T* const dest) { return desul::atomic_dec(const_cast<T*>(dest),desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

template<class T> KOKKOS_INLINE_FUNCTION
void atomic_increment(volatile T* const dest) { return desul::atomic_inc(const_cast<T*>(dest),desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

template<class T> KOKKOS_INLINE_FUNCTION
void atomic_decrement(volatile T* const dest) { return desul::atomic_dec(const_cast<T*>(dest),desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

// Exchange

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_exchange(volatile T* const dest, desul::Impl::dont_deduce_this_parameter_t<const T> val) { return desul::atomic_exchange(const_cast<T*>(dest), val, desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice()); }

template<class T> KOKKOS_INLINE_FUNCTION
bool atomic_compare_exchange_strong(volatile T* const dest, T& expected, const T desired) {
  return desul::atomic_compare_exchange_strong(const_cast<T*>(dest),expected, desired,
                  desul::MemoryOrderRelaxed(), desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice());
}

template<class T> KOKKOS_INLINE_FUNCTION
T atomic_compare_exchange(volatile T* const dest, const T compare, const T desired) {
  return desul::atomic_compare_exchange(const_cast<T*>(dest),compare, desired,
                  desul::MemoryOrderRelaxed(), desul::MemoryScopeDevice());
}

}
// clang-format on
#endif  // KOKKOS_ENABLE_IMPL_DESUL_ATOMICS
#endif
