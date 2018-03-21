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
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#if defined( KOKKOS_ATOMIC_HPP ) && ! defined( KOKKOS_VOLATILE_LOAD_HPP )
#define KOKKOS_VOLATILE_LOAD_HPP

#if defined( __GNUC__ ) /* GNU C   */ || \
    defined( __GNUG__ ) /* GNU C++ */ || \
    defined( __clang__ )

#define KOKKOS_IMPL_MAY_ALIAS __attribute__((__may_alias__))

#else

#define KOKKOS_IMPL_MAY_ALIAS

#endif

namespace Kokkos {

//----------------------------------------------------------------------------

template <typename T>
KOKKOS_FORCEINLINE_FUNCTION
T volatile_load(T const volatile * const src_ptr)
{
  typedef uint64_t KOKKOS_IMPL_MAY_ALIAS T64;
  typedef uint32_t KOKKOS_IMPL_MAY_ALIAS T32;
  typedef uint16_t KOKKOS_IMPL_MAY_ALIAS T16;
  typedef uint8_t  KOKKOS_IMPL_MAY_ALIAS T8;

  enum {
    NUM_8  = sizeof(T),
    NUM_16 = NUM_8 / 2,
    NUM_32 = NUM_8 / 4,
    NUM_64 = NUM_8 / 8
  };

  union {
    T   const volatile * const ptr;
    T64 const volatile * const ptr64;
    T32 const volatile * const ptr32;
    T16 const volatile * const ptr16;
    T8  const volatile * const ptr8;
  } src = {src_ptr};

  T result;

  union {
    T   * const ptr;
    T64 * const ptr64;
    T32 * const ptr32;
    T16 * const ptr16;
    T8  * const ptr8;
  } dst = {&result};

  for (int i=0; i < NUM_64; ++i) {
    dst.ptr64[i] = src.ptr64[i];
  }

  if ( NUM_64*2 < NUM_32 ) {
    dst.ptr32[NUM_64*2] = src.ptr32[NUM_64*2];
  }

  if ( NUM_32*2 < NUM_16 ) {
    dst.ptr16[NUM_32*2] = src.ptr16[NUM_32*2];
  }

  if ( NUM_16*2 < NUM_8 ) {
    dst.ptr8[NUM_16*2] = src.ptr8[NUM_16*2];
  }

  return result;
}

template <typename T>
KOKKOS_FORCEINLINE_FUNCTION
void volatile_store(T volatile * const dst_ptr, T const volatile * const src_ptr)
{
  typedef uint64_t KOKKOS_IMPL_MAY_ALIAS T64;
  typedef uint32_t KOKKOS_IMPL_MAY_ALIAS T32;
  typedef uint16_t KOKKOS_IMPL_MAY_ALIAS T16;
  typedef uint8_t  KOKKOS_IMPL_MAY_ALIAS T8;

  enum {
    NUM_8  = sizeof(T),
    NUM_16 = NUM_8 / 2,
    NUM_32 = NUM_8 / 4,
    NUM_64 = NUM_8 / 8
  };

  union {
    T   const volatile * const ptr;
    T64 const volatile * const ptr64;
    T32 const volatile * const ptr32;
    T16 const volatile * const ptr16;
    T8  const volatile * const ptr8;
  } src = {src_ptr};

  union {
    T   volatile * const ptr;
    T64 volatile * const ptr64;
    T32 volatile * const ptr32;
    T16 volatile * const ptr16;
    T8  volatile * const ptr8;
  } dst = {dst_ptr};

  for (int i=0; i < NUM_64; ++i) {
    dst.ptr64[i] = src.ptr64[i];
  }

  if ( NUM_64*2 < NUM_32 ) {
    dst.ptr32[NUM_64*2] = src.ptr32[NUM_64*2];
  }

  if ( NUM_32*2 < NUM_16 ) {
    dst.ptr16[NUM_32*2] = src.ptr16[NUM_32*2];
  }

  if ( NUM_16*2 < NUM_8 ) {
    dst.ptr8[NUM_16*2] = src.ptr8[NUM_16*2];
  }
}

template <typename T>
KOKKOS_FORCEINLINE_FUNCTION
void volatile_store(T volatile * const dst_ptr, T const * const src_ptr)
{
  typedef uint64_t KOKKOS_IMPL_MAY_ALIAS T64;
  typedef uint32_t KOKKOS_IMPL_MAY_ALIAS T32;
  typedef uint16_t KOKKOS_IMPL_MAY_ALIAS T16;
  typedef uint8_t  KOKKOS_IMPL_MAY_ALIAS T8;

  enum {
    NUM_8  = sizeof(T),
    NUM_16 = NUM_8 / 2,
    NUM_32 = NUM_8 / 4,
    NUM_64 = NUM_8 / 8
  };

  union {
    T   const * const ptr;
    T64 const * const ptr64;
    T32 const * const ptr32;
    T16 const * const ptr16;
    T8  const * const ptr8;
  } src = {src_ptr};

  union {
    T   volatile * const ptr;
    T64 volatile * const ptr64;
    T32 volatile * const ptr32;
    T16 volatile * const ptr16;
    T8  volatile * const ptr8;
  } dst = {dst_ptr};

  for (int i=0; i < NUM_64; ++i) {
    dst.ptr64[i] = src.ptr64[i];
  }

  if ( NUM_64*2 < NUM_32 ) {
    dst.ptr32[NUM_64*2] = src.ptr32[NUM_64*2];
  }

  if ( NUM_32*2 < NUM_16 ) {
    dst.ptr16[NUM_32*2] = src.ptr16[NUM_32*2];
  }

  if ( NUM_16*2 < NUM_8 ) {
    dst.ptr8[NUM_16*2] = src.ptr8[NUM_16*2];
  }
}

template <typename T>
KOKKOS_FORCEINLINE_FUNCTION
void volatile_store(T volatile * dst_ptr, T const volatile & src)
{ volatile_store(dst_ptr, &src); }

template <typename T>
KOKKOS_FORCEINLINE_FUNCTION
void volatile_store(T volatile * dst_ptr, T const & src)
{ volatile_store(dst_ptr, &src); }

template <typename T>
KOKKOS_FORCEINLINE_FUNCTION
T safe_load(T const * const ptr)
{
#if !defined( __MIC__ )
  return *ptr;
#else
  return volatile_load(ptr);
#endif
}

} // namespace kokkos

#undef KOKKOS_IMPL_MAY_ALIAS

#endif

