/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <Kokkos_Macros.hpp>
#if defined(KOKKOS_ATOMIC_HPP) && !defined(KOKKOS_ATOMIC_GENERIC_HPP)
#define KOKKOS_ATOMIC_GENERIC_HPP
#include <Kokkos_Macros.hpp>

#if defined(KOKKOS_ENABLE_CUDA)
#include <Cuda/Kokkos_Cuda_Version_9_8_Compatibility.hpp>
#endif

// Combination operands to be used in an Compare and Exchange based atomic
// operation
namespace Kokkos {
namespace Impl {

template <class Op, class Scalar1, class Scalar2, class Enable = bool>
struct _check_early_exit_impl {
  KOKKOS_FORCEINLINE_FUNCTION
  static constexpr bool check(Op const&, Scalar1 const&,
                              Scalar2 const&) noexcept {
    return false;
  }
};

template <class Op, class Scalar1, class Scalar2>
struct _check_early_exit_impl<
    Op, Scalar1, Scalar2,
    decltype(std::declval<Op const&>().check_early_exit(
        std::declval<Scalar1 const&>(), std::declval<Scalar2 const&>()))> {
  KOKKOS_FORCEINLINE_FUNCTION
  static constexpr bool check(Op const& op, Scalar1 const& v1,
                              Scalar2 const& v2) {
    return op.check_early_exit(v1, v2);
  }
};

template <class Op, class Scalar1, class Scalar2>
KOKKOS_FORCEINLINE_FUNCTION constexpr bool check_early_exit(
    Op const& op, Scalar1 const& v1, Scalar2 const& v2) noexcept {
  return _check_early_exit_impl<Op, Scalar1, Scalar2>::check(op, v1, v2);
}

template <class Scalar1, class Scalar2>
struct MaxOper {
  KOKKOS_FORCEINLINE_FUNCTION
  static Scalar1 apply(const Scalar1& val1, const Scalar2& val2) {
    return (val1 > val2 ? val1 : val2);
  }
  KOKKOS_FORCEINLINE_FUNCTION
  static constexpr bool check_early_exit(Scalar1 const& val1,
                                         Scalar2 const& val2) noexcept {
    return (val1 > val2);
  }
};

template <class Scalar1, class Scalar2>
struct MinOper {
  KOKKOS_FORCEINLINE_FUNCTION
  static Scalar1 apply(const Scalar1& val1, const Scalar2& val2) {
    return (val1 < val2 ? val1 : val2);
  }
  KOKKOS_FORCEINLINE_FUNCTION
  static constexpr bool check_early_exit(Scalar1 const& val1,
                                         Scalar2 const& val2) noexcept {
    return (val1 < val2);
  }
};

template <class Scalar1, class Scalar2>
struct AddOper {
  KOKKOS_FORCEINLINE_FUNCTION
  static Scalar1 apply(const Scalar1& val1, const Scalar2& val2) {
    return val1 + val2;
  }
};

template <class Scalar1, class Scalar2>
struct SubOper {
  KOKKOS_FORCEINLINE_FUNCTION
  static Scalar1 apply(const Scalar1& val1, const Scalar2& val2) {
    return val1 - val2;
  }
};

template <class Scalar1, class Scalar2>
struct MulOper {
  KOKKOS_FORCEINLINE_FUNCTION
  static Scalar1 apply(const Scalar1& val1, const Scalar2& val2) {
    return val1 * val2;
  }
};

template <class Scalar1, class Scalar2>
struct DivOper {
  KOKKOS_FORCEINLINE_FUNCTION
  static Scalar1 apply(const Scalar1& val1, const Scalar2& val2) {
    return val1 / val2;
  }
};

template <class Scalar1, class Scalar2>
struct ModOper {
  KOKKOS_FORCEINLINE_FUNCTION
  static Scalar1 apply(const Scalar1& val1, const Scalar2& val2) {
    return val1 % val2;
  }
};

template <class Scalar1, class Scalar2>
struct AndOper {
  KOKKOS_FORCEINLINE_FUNCTION
  static Scalar1 apply(const Scalar1& val1, const Scalar2& val2) {
    return val1 & val2;
  }
};

template <class Scalar1, class Scalar2>
struct OrOper {
  KOKKOS_FORCEINLINE_FUNCTION
  static Scalar1 apply(const Scalar1& val1, const Scalar2& val2) {
    return val1 | val2;
  }
};

template <class Scalar1, class Scalar2>
struct XorOper {
  KOKKOS_FORCEINLINE_FUNCTION
  static Scalar1 apply(const Scalar1& val1, const Scalar2& val2) {
    return val1 ^ val2;
  }
};

template <class Scalar1, class Scalar2>
struct LShiftOper {
  KOKKOS_FORCEINLINE_FUNCTION
  static Scalar1 apply(const Scalar1& val1, const Scalar2& val2) {
    return val1 << val2;
  }
};

template <class Scalar1, class Scalar2>
struct RShiftOper {
  KOKKOS_FORCEINLINE_FUNCTION
  static Scalar1 apply(const Scalar1& val1, const Scalar2& val2) {
    return val1 >> val2;
  }
};

template <class Oper, typename T>
KOKKOS_INLINE_FUNCTION T atomic_fetch_oper(
    const Oper& op, volatile T* const dest,
    typename std::enable_if<sizeof(T) != sizeof(int) &&
                                sizeof(T) == sizeof(unsigned long long int),
                            const T>::type val) {
  union U {
    unsigned long long int i;
    T t;
    KOKKOS_INLINE_FUNCTION U() {}
  } oldval, assume, newval;

  oldval.t = *dest;

  do {
    if (check_early_exit(op, oldval.t, val)) return oldval.t;
    assume.i = oldval.i;
    newval.t = op.apply(assume.t, val);
    oldval.i = Kokkos::atomic_compare_exchange((unsigned long long int*)dest,
                                               assume.i, newval.i);
  } while (assume.i != oldval.i);

  return oldval.t;
}

template <class Oper, typename T>
KOKKOS_INLINE_FUNCTION T atomic_oper_fetch(
    const Oper& op, volatile T* const dest,
    typename std::enable_if<sizeof(T) != sizeof(int) &&
                                sizeof(T) == sizeof(unsigned long long int),
                            const T>::type val) {
  union U {
    unsigned long long int i;
    T t;
    KOKKOS_INLINE_FUNCTION U() {}
  } oldval, assume, newval;

  oldval.t = *dest;

  do {
    if (check_early_exit(op, oldval.t, val)) return oldval.t;
    assume.i = oldval.i;
    newval.t = op.apply(assume.t, val);
    oldval.i = Kokkos::atomic_compare_exchange((unsigned long long int*)dest,
                                               assume.i, newval.i);
  } while (assume.i != oldval.i);

  return newval.t;
}

template <class Oper, typename T>
KOKKOS_INLINE_FUNCTION T atomic_fetch_oper(
    const Oper& op, volatile T* const dest,
    typename std::enable_if<sizeof(T) == sizeof(int), const T>::type val) {
  union U {
    int i;
    T t;
    KOKKOS_INLINE_FUNCTION U() {}
  } oldval, assume, newval;

  oldval.t = *dest;

  do {
    if (check_early_exit(op, oldval.t, val)) return oldval.t;
    assume.i = oldval.i;
    newval.t = op.apply(assume.t, val);
    oldval.i = Kokkos::atomic_compare_exchange((int*)dest, assume.i, newval.i);
  } while (assume.i != oldval.i);

  return oldval.t;
}

template <class Oper, typename T>
KOKKOS_INLINE_FUNCTION T atomic_oper_fetch(
    const Oper& op, volatile T* const dest,
    typename std::enable_if<sizeof(T) == sizeof(int), const T>::type val) {
  union U {
    int i;
    T t;
    KOKKOS_INLINE_FUNCTION U() {}
  } oldval, assume, newval;

  oldval.t = *dest;

  do {
    if (check_early_exit(op, oldval.t, val)) return oldval.t;
    assume.i = oldval.i;
    newval.t = op.apply(assume.t, val);
    oldval.i = Kokkos::atomic_compare_exchange((int*)dest, assume.i, newval.i);
  } while (assume.i != oldval.i);

  return newval.t;
}

template <class Oper, typename T>
KOKKOS_INLINE_FUNCTION T atomic_fetch_oper(
    const Oper& op, volatile T* const dest,
    typename std::enable_if<(sizeof(T) != 4) && (sizeof(T) != 8), const T>::type
        val) {
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
  while (!Impl::lock_address_host_space((void*)dest))
    ;
  Kokkos::memory_fence();
  T return_val = *dest;
  *dest        = op.apply(return_val, val);
  Kokkos::memory_fence();
  Impl::unlock_address_host_space((void*)dest);
  return return_val;
#elif defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_CUDA)
  // This is a way to (hopefully) avoid dead lock in a warp
  T return_val;
  int done                 = 0;
#ifdef KOKKOS_IMPL_CUDA_SYNCWARP_NEEDS_MASK
  unsigned int mask        = KOKKOS_IMPL_CUDA_ACTIVEMASK;
  unsigned int active      = KOKKOS_IMPL_CUDA_BALLOT_MASK(mask, 1);
#else
  unsigned int active = KOKKOS_IMPL_CUDA_BALLOT(1);
#endif
  unsigned int done_active = 0;
  while (active != done_active) {
    if (!done) {
      if (Impl::lock_address_cuda_space((void*)dest)) {
        Kokkos::memory_fence();
        return_val = *dest;
        *dest      = op.apply(return_val, val);
        Kokkos::memory_fence();
        Impl::unlock_address_cuda_space((void*)dest);
        done = 1;
      }
    }
#ifdef KOKKOS_IMPL_CUDA_SYNCWARP_NEEDS_MASK
    done_active = KOKKOS_IMPL_CUDA_BALLOT_MASK(mask, done);
#else
    done_active = KOKKOS_IMPL_CUDA_BALLOT(done);
#endif
  }
  return return_val;
#elif defined(__HIP_DEVICE_COMPILE__)
  T return_val             = *dest;
  int done                 = 0;
  unsigned int active      = __ballot(1);
  unsigned int done_active = 0;
  while (active != done_active) {
    if (!done) {
      if (Impl::lock_address_hip_space((void*)dest)) {
        return_val = *dest;
        *dest      = op.apply(return_val, val);
        Impl::unlock_address_hip_space((void*)dest);
        done = 1;
      }
    }
    done_active = __ballot(done);
  }
  return return_val;
#elif defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_SYCL)
  // FIXME_SYCL
  Kokkos::abort("Not implemented!");
  (void)op;
  (void)dest;
  (void)val;
  return 0;
#endif
}

template <class Oper, typename T>
KOKKOS_INLINE_FUNCTION T
atomic_oper_fetch(const Oper& op, volatile T* const dest,
                  typename std::enable_if<(sizeof(T) != 4) && (sizeof(T) != 8)
#if defined(KOKKOS_ENABLE_ASM) && \
    defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST)
                                              && (sizeof(T) != 16)
#endif
                                              ,
                                          const T>::type& val) {

#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
  while (!Impl::lock_address_host_space((void*)dest))
    ;
  Kokkos::memory_fence();
  T return_val = op.apply(*dest, val);
  *dest        = return_val;
  Kokkos::memory_fence();
  Impl::unlock_address_host_space((void*)dest);
  return return_val;
#elif defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_CUDA)
  T return_val;
  // This is a way to (hopefully) avoid dead lock in a warp
  int done                 = 0;
#ifdef KOKKOS_IMPL_CUDA_SYNCWARP_NEEDS_MASK
  unsigned int mask        = KOKKOS_IMPL_CUDA_ACTIVEMASK;
  unsigned int active      = KOKKOS_IMPL_CUDA_BALLOT_MASK(mask, 1);
#else
  unsigned int active = KOKKOS_IMPL_CUDA_BALLOT(1);
#endif
  unsigned int done_active = 0;
  while (active != done_active) {
    if (!done) {
      if (Impl::lock_address_cuda_space((void*)dest)) {
        Kokkos::memory_fence();
        return_val = op.apply(*dest, val);
        *dest      = return_val;
        Kokkos::memory_fence();
        Impl::unlock_address_cuda_space((void*)dest);
        done = 1;
      }
    }
#ifdef KOKKOS_IMPL_CUDA_SYNCWARP_NEEDS_MASK
    done_active = KOKKOS_IMPL_CUDA_BALLOT_MASK(mask, done);
#else
    done_active = KOKKOS_IMPL_CUDA_BALLOT(done);
#endif
  }
  return return_val;
#elif defined(__HIP_DEVICE_COMPILE__)
  T return_val;
  int done                 = 0;
  unsigned int active      = __ballot(1);
  unsigned int done_active = 0;
  while (active != done_active) {
    if (!done) {
      if (Impl::lock_address_hip_space((void*)dest)) {
        return_val = op.apply(*dest, val);
        *dest      = return_val;
        Impl::unlock_address_hip_space((void*)dest);
        done = 1;
      }
    }
    done_active = __ballot(done);
  }
  return return_val;
#elif defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_SYCL)
  // FIXME_SYCL
  std::abort();
  (void)op;
  (void)dest;
  (void)val;
  return 0;
#endif
}

}  // namespace Impl
}  // namespace Kokkos

namespace Kokkos {

// Fetch_Oper atomics: return value before operation
template <typename T>
KOKKOS_INLINE_FUNCTION T atomic_fetch_max(volatile T* const dest, const T val) {
  return Impl::atomic_fetch_oper(Impl::MaxOper<T, const T>(), dest, val);
}

template <typename T>
KOKKOS_INLINE_FUNCTION T atomic_fetch_min(volatile T* const dest, const T val) {
  return Impl::atomic_fetch_oper(Impl::MinOper<T, const T>(), dest, val);
}

template <typename T>
KOKKOS_INLINE_FUNCTION T atomic_fetch_mul(volatile T* const dest, const T val) {
  return Impl::atomic_fetch_oper(Impl::MulOper<T, const T>(), dest, val);
}

template <typename T>
KOKKOS_INLINE_FUNCTION T atomic_fetch_div(volatile T* const dest, const T val) {
  return Impl::atomic_fetch_oper(Impl::DivOper<T, const T>(), dest, val);
}

template <typename T>
KOKKOS_INLINE_FUNCTION T atomic_fetch_mod(volatile T* const dest, const T val) {
  return Impl::atomic_fetch_oper(Impl::ModOper<T, const T>(), dest, val);
}

#if !defined(KOKKOS_ENABLE_SERIAL_ATOMICS)

template <typename T>
KOKKOS_INLINE_FUNCTION T atomic_fetch_and(volatile T* const dest, const T val) {
  return Impl::atomic_fetch_oper(Impl::AndOper<T, const T>(), dest, val);
}

template <typename T>
KOKKOS_INLINE_FUNCTION T atomic_fetch_or(volatile T* const dest, const T val) {
  return Impl::atomic_fetch_oper(Impl::OrOper<T, const T>(), dest, val);
}

#endif

template <typename T>
KOKKOS_INLINE_FUNCTION T atomic_fetch_xor(volatile T* const dest, const T val) {
  return Impl::atomic_fetch_oper(Impl::XorOper<T, const T>(), dest, val);
}

template <typename T>
KOKKOS_INLINE_FUNCTION T atomic_fetch_lshift(volatile T* const dest,
                                             const unsigned int val) {
  return Impl::atomic_fetch_oper(Impl::LShiftOper<T, const unsigned int>(),
                                 dest, val);
}

template <typename T>
KOKKOS_INLINE_FUNCTION T atomic_fetch_rshift(volatile T* const dest,
                                             const unsigned int val) {
  return Impl::atomic_fetch_oper(Impl::RShiftOper<T, const unsigned int>(),
                                 dest, val);
}

// Oper Fetch atomics: return value after operation
template <typename T>
KOKKOS_INLINE_FUNCTION T atomic_max_fetch(volatile T* const dest, const T val) {
  return Impl::atomic_oper_fetch(Impl::MaxOper<T, const T>(), dest, val);
}

template <typename T>
KOKKOS_INLINE_FUNCTION T atomic_min_fetch(volatile T* const dest, const T val) {
  return Impl::atomic_oper_fetch(Impl::MinOper<T, const T>(), dest, val);
}

template <typename T>
KOKKOS_INLINE_FUNCTION T atomic_mul_fetch(volatile T* const dest, const T val) {
  return Impl::atomic_oper_fetch(Impl::MulOper<T, const T>(), dest, val);
}

template <typename T>
KOKKOS_INLINE_FUNCTION T atomic_div_fetch(volatile T* const dest, const T val) {
  return Impl::atomic_oper_fetch(Impl::DivOper<T, const T>(), dest, val);
}

template <typename T>
KOKKOS_INLINE_FUNCTION T atomic_mod_fetch(volatile T* const dest, const T val) {
  return Impl::atomic_oper_fetch(Impl::ModOper<T, const T>(), dest, val);
}

template <typename T>
KOKKOS_INLINE_FUNCTION T atomic_and_fetch(volatile T* const dest, const T val) {
  return Impl::atomic_oper_fetch(Impl::AndOper<T, const T>(), dest, val);
}

template <typename T>
KOKKOS_INLINE_FUNCTION T atomic_or_fetch(volatile T* const dest, const T val) {
  return Impl::atomic_oper_fetch(Impl::OrOper<T, const T>(), dest, val);
}

template <typename T>
KOKKOS_INLINE_FUNCTION T atomic_xor_fetch(volatile T* const dest, const T val) {
  return Impl::atomic_oper_fetch(Impl::XorOper<T, const T>(), dest, val);
}

template <typename T>
KOKKOS_INLINE_FUNCTION T atomic_lshift_fetch(volatile T* const dest,
                                             const unsigned int val) {
  return Impl::atomic_oper_fetch(Impl::LShiftOper<T, const unsigned int>(),
                                 dest, val);
}

template <typename T>
KOKKOS_INLINE_FUNCTION T atomic_rshift_fetch(volatile T* const dest,
                                             const unsigned int val) {
  return Impl::atomic_oper_fetch(Impl::RShiftOper<T, const unsigned int>(),
                                 dest, val);
}

#ifdef _WIN32
template <typename T>
KOKKOS_INLINE_FUNCTION T atomic_add_fetch(volatile T* const dest, const T val) {
  return Impl::atomic_oper_fetch(Impl::AddOper<T, const T>(), dest, val);
}

template <typename T>
KOKKOS_INLINE_FUNCTION T atomic_sub_fetch(volatile T* const dest, const T val) {
  return Impl::atomic_oper_fetch(Impl::SubOper<T, const T>(), dest, val);
}

template <typename T>
KOKKOS_INLINE_FUNCTION T atomic_fetch_add(volatile T* const dest, const T val) {
  return Impl::atomic_fetch_oper(Impl::AddOper<T, const T>(), dest, val);
}

template <typename T>
KOKKOS_INLINE_FUNCTION T atomic_fetch_sub(volatile T* const dest, const T val) {
  return Impl::atomic_fetch_oper(Impl::SubOper<T, const T>(), dest, val);
}
#endif

}  // namespace Kokkos
#endif
