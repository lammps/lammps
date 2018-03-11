/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/


#ifndef KOKKOS_SPINWAIT_HPP
#define KOKKOS_SPINWAIT_HPP

#include <Kokkos_Macros.hpp>
#include <Kokkos_Atomic.hpp>

#include <cstdint>

#include <type_traits>

namespace Kokkos {
namespace Impl {

#if defined( KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST )

enum class WaitMode : int {
    ACTIVE   // Used for tight loops to keep threads active longest
  , PASSIVE  // Used to quickly yield the thread to quite down the system
};


void host_thread_yield( const uint32_t i , const WaitMode mode );


template <typename T>
typename std::enable_if< std::is_integral<T>::value, void>::type
spinwait_while_equal( T const volatile & flag, const T value )
{
  Kokkos::store_fence();
  uint32_t i = 0 ;
  while( value == flag ) {
    host_thread_yield(++i, WaitMode::ACTIVE);
  }
  Kokkos::load_fence();
}

template <typename T>
typename std::enable_if< std::is_integral<T>::value, void>::type
yield_while_equal( T const volatile & flag, const T value )
{
  Kokkos::store_fence();
  uint32_t i = 0 ;
  while( value == flag ) {
    host_thread_yield(++i, WaitMode::PASSIVE);
  }
  Kokkos::load_fence();
}

template <typename T>
typename std::enable_if< std::is_integral<T>::value, void>::type
spinwait_until_equal( T const volatile & flag, const T value )
{
  Kokkos::store_fence();
  uint32_t i = 0 ;
  while( value != flag ) {
    host_thread_yield(++i, WaitMode::ACTIVE);
  }
  Kokkos::load_fence();
}

template <typename T>
typename std::enable_if< std::is_integral<T>::value, void>::type
yield_until_equal( T const volatile & flag, const T value )
{
  Kokkos::store_fence();
  uint32_t i = 0 ;
  while( value != flag ) {
    host_thread_yield(++i, WaitMode::PASSIVE);
  }
  Kokkos::load_fence();
}

#else

template <typename T>
KOKKOS_INLINE_FUNCTION
typename std::enable_if< std::is_integral<T>::value, void>::type
spinwait_while_equal( T const volatile & flag, const T value ) {}

template <typename T>
KOKKOS_INLINE_FUNCTION
typename std::enable_if< std::is_integral<T>::value, void>::type
yield_while_equal( T const volatile & flag, const T value ) {}

template <typename T>
KOKKOS_INLINE_FUNCTION
typename std::enable_if< std::is_integral<T>::value, void>::type
spinwait_until_equal( T const volatile & flag, const T value ) {}

template <typename T>
KOKKOS_INLINE_FUNCTION
typename std::enable_if< std::is_integral<T>::value, void>::type
yield_until_equal( T const volatile & flag, const T value ) {}

#endif

} /* namespace Impl */
} /* namespace Kokkos */

#endif /* #ifndef KOKKOS_SPINWAIT_HPP */

