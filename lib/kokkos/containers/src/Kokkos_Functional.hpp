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

#ifndef KOKKOS_FUNCTIONAL_HPP
#define KOKKOS_FUNCTIONAL_HPP

#include <Kokkos_Macros.hpp>
#include <impl/Kokkos_Functional_impl.hpp>

namespace Kokkos {

// These should work for most types

template <typename T>
struct pod_hash
{
  typedef T argument_type;
  typedef T first_argument_type;
  typedef uint32_t second_argument_type;
  typedef uint32_t result_type;

  KOKKOS_FORCEINLINE_FUNCTION
  uint32_t operator()(T const & t) const
  { return Impl::MurmurHash3_x86_32( &t, sizeof(T), 0); }

  KOKKOS_FORCEINLINE_FUNCTION
  uint32_t operator()(T const & t, uint32_t seed) const
  { return Impl::MurmurHash3_x86_32( &t, sizeof(T), seed); }
};

template <typename T>
struct pod_equal_to
{
  typedef T first_argument_type;
  typedef T second_argument_type;
  typedef bool result_type;

  KOKKOS_FORCEINLINE_FUNCTION
  bool operator()(T const & a, T const & b) const
  { return Impl::bitwise_equal(&a,&b); }
};

template <typename T>
struct pod_not_equal_to
{
  typedef T first_argument_type;
  typedef T second_argument_type;
  typedef bool result_type;

  KOKKOS_FORCEINLINE_FUNCTION
  bool operator()(T const & a, T const & b) const
  { return !Impl::bitwise_equal(&a,&b); }
};

template <typename T>
struct equal_to
{
  typedef T first_argument_type;
  typedef T second_argument_type;
  typedef bool result_type;

  KOKKOS_FORCEINLINE_FUNCTION
  bool operator()(T const & a, T const & b) const
  { return a == b; }
};

template <typename T>
struct not_equal_to
{
  typedef T first_argument_type;
  typedef T second_argument_type;
  typedef bool result_type;

  KOKKOS_FORCEINLINE_FUNCTION
  bool operator()(T const & a, T const & b) const
  { return a != b; }
};


template <typename T>
struct greater
{
  typedef T first_argument_type;
  typedef T second_argument_type;
  typedef bool result_type;

  KOKKOS_FORCEINLINE_FUNCTION
  bool operator()(T const & a, T const & b) const
  { return a > b; }
};


template <typename T>
struct less
{
  typedef T first_argument_type;
  typedef T second_argument_type;
  typedef bool result_type;

  KOKKOS_FORCEINLINE_FUNCTION
  bool operator()(T const & a, T const & b) const
  { return a < b; }
};

template <typename T>
struct greater_equal
{
  typedef T first_argument_type;
  typedef T second_argument_type;
  typedef bool result_type;

  KOKKOS_FORCEINLINE_FUNCTION
  bool operator()(T const & a, T const & b) const
  { return a >= b; }
};


template <typename T>
struct less_equal
{
  typedef T first_argument_type;
  typedef T second_argument_type;
  typedef bool result_type;

  KOKKOS_FORCEINLINE_FUNCTION
  bool operator()(T const & a, T const & b) const
  { return a <= b; }
};

} // namespace Kokkos


#endif //KOKKOS_FUNCTIONAL_HPP

