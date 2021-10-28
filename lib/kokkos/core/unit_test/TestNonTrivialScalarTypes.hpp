
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

#ifndef TESTNONTRIVIALSCALARTYPES_HPP_
#define TESTNONTRIVIALSCALARTYPES_HPP_

#include <Kokkos_Core.hpp>

#include <impl/Kokkos_Timer.hpp>
#include <iostream>
#include <cstdlib>
#include <cstdint>
#include <cinttypes>

namespace Test {

struct my_complex {
  double re, im;
  int dummy;

  KOKKOS_INLINE_FUNCTION
  my_complex() {
    re    = 0.0;
    im    = 0.0;
    dummy = 0;
  }

  KOKKOS_INLINE_FUNCTION
  my_complex(const my_complex &src) {
    re    = src.re;
    im    = src.im;
    dummy = src.dummy;
  }

  KOKKOS_INLINE_FUNCTION
  my_complex &operator=(const my_complex &src) {
    re    = src.re;
    im    = src.im;
    dummy = src.dummy;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  my_complex &operator=(const volatile my_complex &src) {
    re    = src.re;
    im    = src.im;
    dummy = src.dummy;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  volatile my_complex &operator=(const my_complex &src) volatile {
    re    = src.re;
    im    = src.im;
    dummy = src.dummy;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  volatile my_complex &operator=(const volatile my_complex &src) volatile {
    re    = src.re;
    im    = src.im;
    dummy = src.dummy;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  my_complex(const volatile my_complex &src) {
    re    = src.re;
    im    = src.im;
    dummy = src.dummy;
  }

  KOKKOS_INLINE_FUNCTION
  my_complex(const double &val) {
    re    = val;
    im    = 0.0;
    dummy = 0;
  }

  KOKKOS_INLINE_FUNCTION
  my_complex &operator+=(const my_complex &src) {
    re += src.re;
    im += src.im;
    dummy += src.dummy;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  void operator+=(const volatile my_complex &src) volatile {
    re += src.re;
    im += src.im;
    dummy += src.dummy;
  }

  KOKKOS_INLINE_FUNCTION
  my_complex operator+(const my_complex &src) {
    my_complex tmp = *this;
    tmp.re += src.re;
    tmp.im += src.im;
    tmp.dummy += src.dummy;
    return tmp;
  }

  KOKKOS_INLINE_FUNCTION
  my_complex operator+(const volatile my_complex &src) volatile {
    my_complex tmp = *this;
    tmp.re += src.re;
    tmp.im += src.im;
    tmp.dummy += src.dummy;
    return tmp;
  }

  KOKKOS_INLINE_FUNCTION
  my_complex &operator*=(const my_complex &src) {
    double re_tmp = re * src.re - im * src.im;
    double im_tmp = re * src.im + im * src.re;
    re            = re_tmp;
    im            = im_tmp;
    dummy *= src.dummy;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  void operator*=(const volatile my_complex &src) volatile {
    double re_tmp = re * src.re - im * src.im;
    double im_tmp = re * src.im + im * src.re;
    re            = re_tmp;
    im            = im_tmp;
    dummy *= src.dummy;
  }

  KOKKOS_INLINE_FUNCTION
  bool operator==(const my_complex &src) const {
    return (re == src.re) && (im == src.im) && (dummy == src.dummy);
  }

  KOKKOS_INLINE_FUNCTION
  bool operator!=(const my_complex &src) const {
    return (re != src.re) || (im != src.im) || (dummy != src.dummy);
  }

  KOKKOS_INLINE_FUNCTION
  bool operator!=(const double &val) const {
    return (re != val) || (im != 0) || (dummy != 0);
  }

  KOKKOS_INLINE_FUNCTION
  my_complex &operator=(const int &val) {
    re    = val;
    im    = 0.0;
    dummy = 0;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  my_complex &operator=(const double &val) {
    re    = val;
    im    = 0.0;
    dummy = 0;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  operator double() { return re; }
};

template <class scalar_t, int N>
struct array_reduce {
  scalar_t data[N];
  KOKKOS_INLINE_FUNCTION
  array_reduce() {
    for (int i = 0; i < N; i++) data[i] = scalar_t();
  }
  KOKKOS_INLINE_FUNCTION
  array_reduce(const array_reduce &rhs) {
    for (int i = 0; i < N; i++) data[i] = rhs.data[i];
  }
  KOKKOS_INLINE_FUNCTION
  array_reduce(const scalar_t value) {
    for (int i = 0; i < N; i++) data[i] = scalar_t(value);
  }

  KOKKOS_INLINE_FUNCTION
  array_reduce &operator=(const array_reduce &src) {
    for (int i = 0; i < N; i++) data[i] = src.data[i];
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  array_reduce &operator=(const volatile array_reduce &src) {
    for (int i = 0; i < N; i++) data[i] = src.data[i];
    return *this;
  }

  KOKKOS_INLINE_FUNCTION  // add operator
      array_reduce &
      operator=(const scalar_t val) {
    for (int i = 0; i < N; i++) data[i] = val;
    return *this;
  }
  KOKKOS_INLINE_FUNCTION  // add operator
      array_reduce &
      operator=(const int val) {
    for (int i = 0; i < N; i++) data[i] = val;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION  // add operator
      array_reduce &
      operator+=(const array_reduce &src) {
    for (int i = 0; i < N; i++) data[i] += src.data[i];
    return *this;
  }
  KOKKOS_INLINE_FUNCTION  // volatile add operator
      void
      operator+=(const volatile array_reduce &src) volatile {
    for (int i = 0; i < N; i++) data[i] += src.data[i];
  }
  KOKKOS_INLINE_FUNCTION  // add operator
      array_reduce
      operator+(const array_reduce &src) const {
    array_reduce result(*this);
    for (int i = 0; i < N; i++) result.data[i] += src.data[i];
    return result;
  }
  KOKKOS_INLINE_FUNCTION  // add operator
      array_reduce
      operator-(const array_reduce &src) const {
    array_reduce result(*this);
    for (int i = 0; i < N; i++) result.data[i] -= src.data[i];
    return result;
  }
  KOKKOS_INLINE_FUNCTION  // add operator
      array_reduce &
      operator*=(const array_reduce &src) {
    for (int i = 0; i < N; i++) data[i] *= src.data[i];
    return *this;
  }
  KOKKOS_INLINE_FUNCTION  // volatile add operator
      void
      operator*=(const volatile array_reduce &src) volatile {
    for (int i = 0; i < N; i++) data[i] *= src.data[i];
  }
  KOKKOS_INLINE_FUNCTION  // add operator
      array_reduce
      operator*(const array_reduce &src) const {
    array_reduce result(*this);
    for (int i = 0; i < N; i++) result.data[i] *= src.data[i];
    return result;
  }
  KOKKOS_INLINE_FUNCTION
  bool operator==(const array_reduce &src) const {
    bool equal = true;
    for (int i = 0; i < N; i++) equal = equal && (data[i] == src.data[i]);
    return equal;
  }
  KOKKOS_INLINE_FUNCTION
  bool operator!=(const array_reduce &src) const {
    bool equal = true;
    for (int i = 0; i < N; i++) equal = equal && (data[i] == src.data[i]);
    return !equal;
  }
  KOKKOS_INLINE_FUNCTION
  explicit operator double() const {
    double lsum = 0.0;
    for (int i = 0; i < N; i++) lsum += data[i];
    return lsum;
  }
};
}  // namespace Test

namespace Kokkos {
template <>
struct reduction_identity<Test::my_complex> {
  using t_red_ident = reduction_identity<double>;
  KOKKOS_FORCEINLINE_FUNCTION static Test::my_complex sum() {
    return Test::my_complex(t_red_ident::sum());
  }
  KOKKOS_FORCEINLINE_FUNCTION static Test::my_complex prod() {
    return Test::my_complex(t_red_ident::prod());
  }
};

template <class scalar_t, int N>
struct reduction_identity<Test::array_reduce<scalar_t, N>> {
  using t_red_ident = reduction_identity<scalar_t>;
  KOKKOS_FORCEINLINE_FUNCTION static Test::array_reduce<scalar_t, N> sum() {
    return Test::array_reduce<scalar_t, N>(t_red_ident::sum());
  }
  KOKKOS_FORCEINLINE_FUNCTION static Test::array_reduce<scalar_t, N> prod() {
    return Test::array_reduce<scalar_t, N>(t_red_ident::prod());
  }
};
}  // namespace Kokkos
#endif  // TESTNONTRIVIALSCALARTYPES_HPP_
