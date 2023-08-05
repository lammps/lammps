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

#ifndef TESTNONTRIVIALSCALARTYPES_HPP_
#define TESTNONTRIVIALSCALARTYPES_HPP_

#include <Kokkos_Core.hpp>

#include <Kokkos_Timer.hpp>
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
  my_complex operator+(const my_complex &src) {
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

struct point_t {
  uint8_t x, y, z;

  KOKKOS_FUNCTION
  point_t() : x(1), y(1), z(1){};

  KOKKOS_FUNCTION
  point_t(const point_t &val) : x(val.x), y(val.y), z(val.z){};

  KOKKOS_FUNCTION
  point_t(const int rhs) { x = y = z = static_cast<uint8_t>(rhs); }

  KOKKOS_FUNCTION
  explicit operator int() const { return static_cast<int>(x + y + z); }

  KOKKOS_FUNCTION
  bool operator==(const point_t rhs) const {
    return (x == rhs.x && y == rhs.y && z == rhs.z);
  }

  KOKKOS_FUNCTION
  void operator=(point_t rhs) {
    x = rhs.x;
    y = rhs.y;
    z = rhs.z;
  }

  KOKKOS_FUNCTION
  point_t operator+=(const point_t rhs) {
    x += rhs.x;
    y += rhs.y;
    z += rhs.z;
    return *this;
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

template <>
struct reduction_identity<Test::point_t> {
  KOKKOS_FORCEINLINE_FUNCTION constexpr static uint8_t sum() noexcept {
    return 0;
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static uint8_t prod() noexcept {
    return 1;
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static uint8_t max() noexcept {
    return 0xff;
  }
  KOKKOS_FORCEINLINE_FUNCTION constexpr static uint8_t min() noexcept {
    return 0x0;
  }
};
}  // namespace Kokkos
#endif  // TESTNONTRIVIALSCALARTYPES_HPP_
