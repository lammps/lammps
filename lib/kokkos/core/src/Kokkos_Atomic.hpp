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

/// \file Kokkos_Atomic.hpp
/// \brief Atomic functions
///
/// This header file defines prototypes for the following atomic functions:
///   - exchange
///   - compare and exchange
///   - add
///
/// Supported types include:
///   - signed and unsigned 4 and 8 byte integers
///   - float
///   - double
///
/// They are implemented through GCC compatible intrinsics, OpenMP
/// directives and native CUDA intrinsics.
///
/// Including this header file requires one of the following
/// compilers:
///   - NVCC (for CUDA device code only)
///   - GCC (for host code only)
///   - Intel (for host code only)
///   - A compiler that supports OpenMP 3.1 (for host code only)

#ifndef KOKKOS_ATOMIC_HPP
#define KOKKOS_ATOMIC_HPP
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_ATOMIC
#endif

#include <Kokkos_Macros.hpp>

#ifdef KOKKOS_ENABLE_IMPL_DESUL_ATOMICS
#include <Kokkos_Atomics_Desul_Wrapper.hpp>
#include <Kokkos_Atomics_Desul_Volatile_Wrapper.hpp>
#include <impl/Kokkos_Utilities.hpp>

// Helper functions for places where we really should have called SeqCst atomics
// anyway These can go away when we call desul unconditionally Non-Desul
// versions are below
namespace Kokkos {
namespace Impl {
using desul::MemoryOrderSeqCst;
using desul::MemoryScopeDevice;

template <class T>
KOKKOS_INLINE_FUNCTION void desul_atomic_dec(T* dest, MemoryOrderSeqCst,
                                             MemoryScopeDevice) {
  return desul::atomic_dec(const_cast<T*>(dest), desul::MemoryOrderSeqCst(),
                           desul::MemoryScopeDevice());
}

template <class T>
KOKKOS_INLINE_FUNCTION void desul_atomic_inc(T* dest, MemoryOrderSeqCst,
                                             MemoryScopeDevice) {
  return desul::atomic_inc(const_cast<T*>(dest), desul::MemoryOrderSeqCst(),
                           desul::MemoryScopeDevice());
}

template <class T>
KOKKOS_INLINE_FUNCTION T
desul_atomic_exchange(T* dest, const Kokkos::Impl::identity_t<T> val,
                      MemoryOrderSeqCst, MemoryScopeDevice) {
  return desul::atomic_exchange(const_cast<T*>(dest), val,
                                desul::MemoryOrderSeqCst(),
                                desul::MemoryScopeDevice());
}

template <class T>
KOKKOS_INLINE_FUNCTION T desul_atomic_compare_exchange(
    T* dest, Kokkos::Impl::identity_t<const T> compare,
    Kokkos::Impl::identity_t<const T> val, MemoryOrderSeqCst,
    MemoryScopeDevice) {
  return desul::atomic_compare_exchange(dest, compare, val,
                                        desul::MemoryOrderSeqCst(),
                                        desul::MemoryScopeDevice());
}

}  // namespace Impl
}  // namespace Kokkos
#else

#include <Kokkos_HostSpace.hpp>
#include <impl/Kokkos_Traits.hpp>

//----------------------------------------------------------------------------

// Need to fix this for pure clang on windows
#if defined(_WIN32)
#define KOKKOS_ENABLE_WINDOWS_ATOMICS

#if defined(KOKKOS_ENABLE_CUDA)
#define KOKKOS_ENABLE_CUDA_ATOMICS
#if defined(KOKKOS_COMPILER_CLANG)
#define KOKKOS_ENABLE_GNU_ATOMICS
#endif
#endif

#else  // _WIN32
#if defined(KOKKOS_ENABLE_CUDA)

// Compiling NVIDIA device code, must use Cuda atomics:

#define KOKKOS_ENABLE_CUDA_ATOMICS

#elif defined(KOKKOS_ENABLE_HIP)

#define KOKKOS_ENABLE_HIP_ATOMICS

#endif

#if !defined(KOKKOS_ENABLE_GNU_ATOMICS) &&    \
    !defined(KOKKOS_ENABLE_INTEL_ATOMICS) &&  \
    !defined(KOKKOS_ENABLE_OPENMP_ATOMICS) && \
    !defined(KOKKOS_ENABLE_STD_ATOMICS) &&    \
    !defined(KOKKOS_ENABLE_SERIAL_ATOMICS)

// Compiling for non-Cuda atomic implementation has not been pre-selected.
// Choose the best implementation for the detected compiler.
// Preference: GCC, INTEL, OMP31

#if defined(KOKKOS_INTERNAL_NOT_PARALLEL)

#define KOKKOS_ENABLE_SERIAL_ATOMICS

#elif defined(KOKKOS_COMPILER_GNU) || defined(KOKKOS_COMPILER_CLANG) || \
    (defined(KOKKOS_COMPILER_NVCC) || defined(KOKKOS_COMPILER_IBM))

#define KOKKOS_ENABLE_GNU_ATOMICS

#elif defined(KOKKOS_COMPILER_INTEL) || defined(KOKKOS_COMPILER_CRAYC)

#define KOKKOS_ENABLE_INTEL_ATOMICS

#elif defined(_OPENMP) && (201107 <= _OPENMP)

#define KOKKOS_ENABLE_OPENMP_ATOMICS

#else

#error "KOKKOS_ATOMICS_USE : Unsupported compiler"

#endif

#endif /* Not pre-selected atomic implementation */
#endif

#ifdef KOKKOS_ENABLE_CUDA
#include <Cuda/Kokkos_Cuda_Locks.hpp>
#endif

namespace Kokkos {
template <typename T>
KOKKOS_INLINE_FUNCTION void atomic_add(volatile T* const dest, const T src);

// Atomic increment
template <typename T>
KOKKOS_INLINE_FUNCTION void atomic_increment(volatile T* a);

template <typename T>
KOKKOS_INLINE_FUNCTION void atomic_decrement(volatile T* a);
}  // namespace Kokkos

namespace Kokkos {

inline const char* atomic_query_version() {
#if defined(KOKKOS_ENABLE_CUDA_ATOMICS)
  return "KOKKOS_ENABLE_CUDA_ATOMICS";
#elif defined(KOKKOS_ENABLE_GNU_ATOMICS)
  return "KOKKOS_ENABLE_GNU_ATOMICS";
#elif defined(KOKKOS_ENABLE_INTEL_ATOMICS)
  return "KOKKOS_ENABLE_INTEL_ATOMICS";
#elif defined(KOKKOS_ENABLE_OPENMP_ATOMICS)
  return "KOKKOS_ENABLE_OPENMP_ATOMICS";
#elif defined(KOKKOS_ENABLE_WINDOWS_ATOMICS)
  return "KOKKOS_ENABLE_WINDOWS_ATOMICS";
#elif defined(KOKKOS_ENABLE_SERIAL_ATOMICS)
  return "KOKKOS_ENABLE_SERIAL_ATOMICS";
#else
#error "No valid response for atomic_query_version!"
#endif
}

}  // namespace Kokkos

//----------------------------------------------------------------------------
// Atomic Memory Orders
//
// Implements Strongly-typed analogs of C++ standard memory orders
#include "impl/Kokkos_Atomic_Memory_Order.hpp"

#if defined(KOKKOS_ENABLE_HIP)
#include <HIP/Kokkos_HIP_Atomic.hpp>
#endif

#if defined(KOKKOS_ENABLE_WINDOWS_ATOMICS)
#include "impl/Kokkos_Atomic_Windows.hpp"
#endif
//----------------------------------------------------------------------------
// Atomic Assembly
//
// Implements CAS128-bit in assembly

#include "impl/Kokkos_Atomic_Assembly.hpp"

//----------------------------------------------------------------------------
// Memory fence
//
// All loads and stores from this thread will be globally consistent before
// continuing
//
// void memory_fence() {...};
#include "impl/Kokkos_Memory_Fence.hpp"

//----------------------------------------------------------------------------
// Atomic exchange
//
// template< typename T >
// T atomic_exchange( volatile T* const dest , const T val )
// { T tmp = *dest ; *dest = val ; return tmp ; }

#include "impl/Kokkos_Atomic_Exchange.hpp"

//----------------------------------------------------------------------------
// Atomic compare-and-exchange
//
// template<class T>
// bool atomic_compare_exchange_strong(volatile T* const dest, const T compare,
// const T val) { bool equal = compare == *dest ; if ( equal ) { *dest = val ; }
// return equal ; }

#include "impl/Kokkos_Atomic_Compare_Exchange_Strong.hpp"

#include "impl/Kokkos_Atomic_Generic.hpp"

//----------------------------------------------------------------------------
// Atomic fetch and add
//
// template<class T>
// T atomic_fetch_add(volatile T* const dest, const T val)
// { T tmp = *dest ; *dest += val ; return tmp ; }

#include "impl/Kokkos_Atomic_Fetch_Add.hpp"

//----------------------------------------------------------------------------
// Atomic increment
//
// template<class T>
// T atomic_increment(volatile T* const dest)
// { dest++; }

#include "impl/Kokkos_Atomic_Increment.hpp"

//----------------------------------------------------------------------------
// Atomic Decrement
//
// template<class T>
// T atomic_decrement(volatile T* const dest)
// { dest--; }

#include "impl/Kokkos_Atomic_Decrement.hpp"

//----------------------------------------------------------------------------
// Atomic fetch and sub
//
// template<class T>
// T atomic_fetch_sub(volatile T* const dest, const T val)
// { T tmp = *dest ; *dest -= val ; return tmp ; }

#include "impl/Kokkos_Atomic_Fetch_Sub.hpp"

//----------------------------------------------------------------------------
// Atomic fetch and or
//
// template<class T>
// T atomic_fetch_or(volatile T* const dest, const T val)
// { T tmp = *dest ; *dest = tmp | val ; return tmp ; }

#include "impl/Kokkos_Atomic_Fetch_Or.hpp"

//----------------------------------------------------------------------------
// Atomic fetch and and
//
// template<class T>
// T atomic_fetch_and(volatile T* const dest, const T val)
// { T tmp = *dest ; *dest = tmp & val ; return tmp ; }

#include "impl/Kokkos_Atomic_Fetch_And.hpp"

//----------------------------------------------------------------------------
// Atomic MinMax
//
// template<class T>
// T atomic_min(volatile T* const dest, const T val)
// { T tmp = *dest ; *dest = min(*dest, val); return tmp ; }
// template<class T>
// T atomic_max(volatile T* const dest, const T val)
// { T tmp = *dest ; *dest = max(*dest, val); return tmp ; }

#include "impl/Kokkos_Atomic_MinMax.hpp"

//----------------------------------------------------------------------------
// Provide volatile_load and safe_load
//
// T volatile_load(T const volatile * const ptr);
//
// T const& safe_load(T const * const ptr);
// XEON PHI
// T safe_load(T const * const ptr

#include "impl/Kokkos_Volatile_Load.hpp"

//----------------------------------------------------------------------------
// Provide atomic loads and stores with memory order semantics

#include "impl/Kokkos_Atomic_Load.hpp"
#include "impl/Kokkos_Atomic_Store.hpp"

// Generic functions using the above defined functions
#include "impl/Kokkos_Atomic_Generic_Secondary.hpp"
//----------------------------------------------------------------------------
// This atomic-style macro should be an inlined function, not a macro

#if defined(KOKKOS_COMPILER_GNU) && !defined(__PGIC__) && \
    !defined(__CUDA_ARCH__)

#define KOKKOS_NONTEMPORAL_PREFETCH_LOAD(addr) __builtin_prefetch(addr, 0, 0)
#define KOKKOS_NONTEMPORAL_PREFETCH_STORE(addr) __builtin_prefetch(addr, 1, 0)

#else

#define KOKKOS_NONTEMPORAL_PREFETCH_LOAD(addr) ((void)0)
#define KOKKOS_NONTEMPORAL_PREFETCH_STORE(addr) ((void)0)

#endif

//----------------------------------------------------------------------------

// Helper functions for places where we really should have called SeqCst atomics
// anyway These can go away when we call desul unconditionally
namespace Kokkos {
namespace Impl {
struct MemoryOrderSeqCst {};
struct MemoryScopeDevice {};

template <class T>
KOKKOS_INLINE_FUNCTION void desul_atomic_dec(T* dest, MemoryOrderSeqCst,
                                             MemoryScopeDevice) {
  return Kokkos::atomic_decrement(dest);
}

template <class T>
KOKKOS_INLINE_FUNCTION void desul_atomic_inc(T* dest, MemoryOrderSeqCst,
                                             MemoryScopeDevice) {
  return Kokkos::atomic_increment(dest);
}

template <class T>
KOKKOS_INLINE_FUNCTION T
desul_atomic_exchange(T* dest, Kokkos::Impl::identity_t<const T> val,
                      MemoryOrderSeqCst, MemoryScopeDevice) {
  return Kokkos::atomic_exchange(dest, val);
}

template <class T>
KOKKOS_INLINE_FUNCTION T desul_atomic_compare_exchange(
    T* dest, Kokkos::Impl::identity_t<const T> compare,
    Kokkos::Impl::identity_t<const T> val, MemoryOrderSeqCst,
    MemoryScopeDevice) {
  return Kokkos::atomic_compare_exchange(dest, compare, val);
}

}  // namespace Impl
}  // namespace Kokkos

#endif /* !KOKKOS_ENABLE_IMPL_DESUL_ATOMICS */
#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_ATOMIC
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_ATOMIC
#endif
#endif /* KOKKOS_ATOMIC_HPP */
