/*
//@HEADER
// ************************************************************************
//
//                             Kokkos
//         Manycore Performance-Portable Multidimensional Arrays
//
//              Copyright (2012) Sandia Corporation
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
// Questions?  Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/
#if defined( KOKKOS_ATOMIC_HPP ) && ! defined( KOKKOS_ATOMIC_ASSEMBLY_X86_HPP )
#define KOKKOS_ATOMIC_ASSEMBLY_X86_HPP
namespace Kokkos {

#ifndef __CUDA_ARCH__
template<>
KOKKOS_INLINE_FUNCTION
void atomic_increment<char>(volatile char* a) {
  __asm__ __volatile__(
    "lock incb %0"
    : /* no output registers */
    : "m" (a[0])
    : "memory"
  );
}

template<>
KOKKOS_INLINE_FUNCTION
void atomic_increment<short>(volatile short* a) {
  __asm__ __volatile__(
    "lock incw %0"
    : /* no output registers */
    : "m" (a[0])
    : "memory"
  );
}

template<>
KOKKOS_INLINE_FUNCTION
void atomic_increment<int>(volatile int* a) {
  __asm__ __volatile__(
    "lock incl %0"
    : /* no output registers */
    : "m" (a[0])
    : "memory"
  );
}

template<>
KOKKOS_INLINE_FUNCTION
void atomic_increment<long long int>(volatile long long int* a) {
  __asm__ __volatile__(
    "lock incq %0"
    : /* no output registers */
    : "m" (a[0])
    : "memory"
  );
}

template<>
KOKKOS_INLINE_FUNCTION
void atomic_decrement<char>(volatile char* a) {
  __asm__ __volatile__(
    "lock decb %0"
    : /* no output registers */
    : "m" (a[0])
    : "memory"
  );
}

template<>
KOKKOS_INLINE_FUNCTION
void atomic_decrement<short>(volatile short* a) {
  __asm__ __volatile__(
    "lock decw %0"
    : /* no output registers */
    : "m" (a[0])
    : "memory"
  );
}

template<>
KOKKOS_INLINE_FUNCTION
void atomic_decrement<int>(volatile int* a) {
  __asm__ __volatile__(
    "lock decl %0"
    : /* no output registers */
    : "m" (a[0])
    : "memory"
  );
}

template<>
KOKKOS_INLINE_FUNCTION
void atomic_decrement<long long int>(volatile long long int* a) {
  __asm__ __volatile__(
    "lock decq %0"
    : /* no output registers */
    : "m" (a[0])
    : "memory"
  );
}
#endif
}

#endif
