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

#ifndef KOKKOS_HALF_FLOATING_POINT_WRAPPER_HPP_
#define KOKKOS_HALF_FLOATING_POINT_WRAPPER_HPP_

#include <Kokkos_Macros.hpp>

#include <type_traits>
#include <iosfwd>  // istream & ostream for extraction and insertion ops
#include <string>

#ifdef KOKKOS_IMPL_HALF_TYPE_DEFINED

// KOKKOS_HALF_IS_FULL_TYPE_ON_ARCH: A macro to select which
// floating_pointer_wrapper operator paths should be used. For CUDA, let the
// compiler conditionally select when device ops are used For SYCL, we have a
// full half type on both host and device
#if defined(__CUDA_ARCH__) || defined(KOKKOS_ENABLE_SYCL)
#define KOKKOS_HALF_IS_FULL_TYPE_ON_ARCH
#endif

/************************* BEGIN forward declarations *************************/
namespace Kokkos {
namespace Experimental {
namespace Impl {
template <class FloatType>
class floating_point_wrapper;
}

// Declare half_t (binary16)
using half_t = Kokkos::Experimental::Impl::floating_point_wrapper<
    Kokkos::Impl::half_impl_t ::type>;
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(float val);
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(bool val);
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(double val);
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(short val);
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(int val);
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(long val);
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(long long val);
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(unsigned short val);
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(unsigned int val);
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(unsigned long val);
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(unsigned long long val);
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(half_t);

template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, float>::value, T>
    cast_from_half(half_t);
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, bool>::value, T>
    cast_from_half(half_t);
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, double>::value, T>
    cast_from_half(half_t);
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, short>::value, T>
    cast_from_half(half_t);
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, int>::value, T>
    cast_from_half(half_t);
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, long>::value, T>
    cast_from_half(half_t);
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, long long>::value, T>
    cast_from_half(half_t);
template <class T>
KOKKOS_INLINE_FUNCTION
    std::enable_if_t<std::is_same<T, unsigned short>::value, T>
        cast_from_half(half_t);
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, unsigned int>::value, T>
    cast_from_half(half_t);
template <class T>
KOKKOS_INLINE_FUNCTION
    std::enable_if_t<std::is_same<T, unsigned long>::value, T>
        cast_from_half(half_t);
template <class T>
KOKKOS_INLINE_FUNCTION
    std::enable_if_t<std::is_same<T, unsigned long long>::value, T>
        cast_from_half(half_t);

// declare bhalf_t
#ifdef KOKKOS_IMPL_BHALF_TYPE_DEFINED
using bhalf_t = Kokkos::Experimental::Impl::floating_point_wrapper<
    Kokkos::Impl ::bhalf_impl_t ::type>;

KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(float val);
KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(bool val);
KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(double val);
KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(short val);
KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(int val);
KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(long val);
KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(long long val);
KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(unsigned short val);
KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(unsigned int val);
KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(unsigned long val);
KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(unsigned long long val);
KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(bhalf_t val);

template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, float>::value, T>
    cast_from_bhalf(bhalf_t);
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, bool>::value, T>
    cast_from_bhalf(bhalf_t);
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, double>::value, T>
    cast_from_bhalf(bhalf_t);
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, short>::value, T>
    cast_from_bhalf(bhalf_t);
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, int>::value, T>
    cast_from_bhalf(bhalf_t);
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, long>::value, T>
    cast_from_bhalf(bhalf_t);
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, long long>::value, T>
    cast_from_bhalf(bhalf_t);
template <class T>
KOKKOS_INLINE_FUNCTION
    std::enable_if_t<std::is_same<T, unsigned short>::value, T>
        cast_from_bhalf(bhalf_t);
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_same<T, unsigned int>::value, T>
    cast_from_bhalf(bhalf_t);
template <class T>
KOKKOS_INLINE_FUNCTION
    std::enable_if_t<std::is_same<T, unsigned long>::value, T>
        cast_from_bhalf(bhalf_t);
template <class T>
KOKKOS_INLINE_FUNCTION
    std::enable_if_t<std::is_same<T, unsigned long long>::value, T>
        cast_from_bhalf(bhalf_t);
#endif  // KOKKOS_IMPL_BHALF_TYPE_DEFINED

template <class T>
static KOKKOS_INLINE_FUNCTION Kokkos::Experimental::half_t cast_to_wrapper(
    T x, const volatile Kokkos::Impl::half_impl_t::type&);

#ifdef KOKKOS_IMPL_BHALF_TYPE_DEFINED
template <class T>
static KOKKOS_INLINE_FUNCTION Kokkos::Experimental::bhalf_t cast_to_wrapper(
    T x, const volatile Kokkos::Impl::bhalf_impl_t::type&);
#endif  // KOKKOS_IMPL_BHALF_TYPE_DEFINED

template <class T>
static KOKKOS_INLINE_FUNCTION T
cast_from_wrapper(const Kokkos::Experimental::half_t& x);

#ifdef KOKKOS_IMPL_BHALF_TYPE_DEFINED
template <class T>
static KOKKOS_INLINE_FUNCTION T
cast_from_wrapper(const Kokkos::Experimental::bhalf_t& x);
#endif  // KOKKOS_IMPL_BHALF_TYPE_DEFINED
/************************** END forward declarations **************************/

namespace Impl {
template <class FloatType>
class alignas(FloatType) floating_point_wrapper {
 public:
  using impl_type = FloatType;

 private:
  impl_type val;
  using fixed_width_integer_type = std::conditional_t<
      sizeof(impl_type) == 2, uint16_t,
      std::conditional_t<
          sizeof(impl_type) == 4, uint32_t,
          std::conditional_t<sizeof(impl_type) == 8, uint64_t, void>>>;
  static_assert(!std::is_void<fixed_width_integer_type>::value,
                "Invalid impl_type");

 public:
  // In-class initialization and defaulted default constructors not used
  // since Cuda supports half precision initialization via the below constructor
  KOKKOS_FUNCTION
  floating_point_wrapper() : val(0.0F) {}

// Copy constructors
// Getting "C2580: multiple versions of a defaulted special
// member function are not allowed" with VS 16.11.3 and CUDA 11.4.2
#if defined(_WIN32) && defined(KOKKOS_ENABLE_CUDA)
  KOKKOS_FUNCTION
  floating_point_wrapper(const floating_point_wrapper& rhs) : val(rhs.val) {}

  KOKKOS_FUNCTION
  floating_point_wrapper& operator=(const floating_point_wrapper& rhs) {
    val = rhs.val;
    return *this;
  }
#else
  KOKKOS_DEFAULTED_FUNCTION
  floating_point_wrapper(const floating_point_wrapper&) noexcept = default;

  KOKKOS_DEFAULTED_FUNCTION
  floating_point_wrapper& operator=(const floating_point_wrapper&) noexcept =
      default;
#endif

  KOKKOS_INLINE_FUNCTION
  floating_point_wrapper(const volatile floating_point_wrapper& rhs) {
#if defined(KOKKOS_HALF_IS_FULL_TYPE_ON_ARCH) && !defined(KOKKOS_ENABLE_SYCL)
    val = rhs.val;
#else
    const volatile fixed_width_integer_type* rv_ptr =
        reinterpret_cast<const volatile fixed_width_integer_type*>(&rhs.val);
    const fixed_width_integer_type rv_val = *rv_ptr;
    val       = reinterpret_cast<const impl_type&>(rv_val);
#endif  // KOKKOS_HALF_IS_FULL_TYPE_ON_ARCH
  }

  // Don't support implicit conversion back to impl_type.
  // impl_type is a storage only type on host.
  KOKKOS_FUNCTION
  explicit operator impl_type() const { return val; }
  KOKKOS_FUNCTION
  explicit operator float() const { return cast_from_wrapper<float>(*this); }
  KOKKOS_FUNCTION
  explicit operator bool() const { return cast_from_wrapper<bool>(*this); }
  KOKKOS_FUNCTION
  explicit operator double() const { return cast_from_wrapper<double>(*this); }
  KOKKOS_FUNCTION
  explicit operator short() const { return cast_from_wrapper<short>(*this); }
  KOKKOS_FUNCTION
  explicit operator int() const { return cast_from_wrapper<int>(*this); }
  KOKKOS_FUNCTION
  explicit operator long() const { return cast_from_wrapper<long>(*this); }
  KOKKOS_FUNCTION
  explicit operator long long() const {
    return cast_from_wrapper<long long>(*this);
  }
  KOKKOS_FUNCTION
  explicit operator unsigned short() const {
    return cast_from_wrapper<unsigned short>(*this);
  }
  KOKKOS_FUNCTION
  explicit operator unsigned int() const {
    return cast_from_wrapper<unsigned int>(*this);
  }
  KOKKOS_FUNCTION
  explicit operator unsigned long() const {
    return cast_from_wrapper<unsigned long>(*this);
  }
  KOKKOS_FUNCTION
  explicit operator unsigned long long() const {
    return cast_from_wrapper<unsigned long long>(*this);
  }

  /**
   * Conversion constructors.
   *
   * Support implicit conversions from impl_type, float, double ->
   * floating_point_wrapper. Mixed precision expressions require upcasting which
   * is done in the
   * "// Binary Arithmetic" operator overloads below.
   *
   * Support implicit conversions from integral types -> floating_point_wrapper.
   * Expressions involving floating_point_wrapper with integral types require
   * downcasting the integral types to floating_point_wrapper. Existing operator
   * overloads can handle this with the addition of the below implicit
   * conversion constructors.
   */
  KOKKOS_FUNCTION
  constexpr floating_point_wrapper(impl_type rhs) : val(rhs) {}
  KOKKOS_FUNCTION
  floating_point_wrapper(float rhs) : val(cast_to_wrapper(rhs, val).val) {}
  KOKKOS_FUNCTION
  floating_point_wrapper(double rhs) : val(cast_to_wrapper(rhs, val).val) {}
  KOKKOS_FUNCTION
  explicit floating_point_wrapper(bool rhs)
      : val(cast_to_wrapper(rhs, val).val) {}
  KOKKOS_FUNCTION
  floating_point_wrapper(short rhs) : val(cast_to_wrapper(rhs, val).val) {}
  KOKKOS_FUNCTION
  floating_point_wrapper(int rhs) : val(cast_to_wrapper(rhs, val).val) {}
  KOKKOS_FUNCTION
  floating_point_wrapper(long rhs) : val(cast_to_wrapper(rhs, val).val) {}
  KOKKOS_FUNCTION
  floating_point_wrapper(long long rhs) : val(cast_to_wrapper(rhs, val).val) {}
  KOKKOS_FUNCTION
  floating_point_wrapper(unsigned short rhs)
      : val(cast_to_wrapper(rhs, val).val) {}
  KOKKOS_FUNCTION
  floating_point_wrapper(unsigned int rhs)
      : val(cast_to_wrapper(rhs, val).val) {}
  KOKKOS_FUNCTION
  floating_point_wrapper(unsigned long rhs)
      : val(cast_to_wrapper(rhs, val).val) {}
  KOKKOS_FUNCTION
  floating_point_wrapper(unsigned long long rhs)
      : val(cast_to_wrapper(rhs, val).val) {}

  // Unary operators
  KOKKOS_FUNCTION
  floating_point_wrapper operator+() const {
    floating_point_wrapper tmp = *this;
#ifdef KOKKOS_HALF_IS_FULL_TYPE_ON_ARCH
    tmp.val = +tmp.val;
#else
    tmp.val   = cast_to_wrapper(+cast_from_wrapper<float>(tmp), val).val;
#endif
    return tmp;
  }

  KOKKOS_FUNCTION
  floating_point_wrapper operator-() const {
    floating_point_wrapper tmp = *this;
#ifdef KOKKOS_HALF_IS_FULL_TYPE_ON_ARCH
    tmp.val = -tmp.val;
#else
    tmp.val   = cast_to_wrapper(-cast_from_wrapper<float>(tmp), val).val;
#endif
    return tmp;
  }

  // Prefix operators
  KOKKOS_FUNCTION
  floating_point_wrapper& operator++() {
#ifdef KOKKOS_HALF_IS_FULL_TYPE_ON_ARCH
    val = val + impl_type(1.0F);  // cuda has no operator++ for __nv_bfloat
#else
    float tmp = cast_from_wrapper<float>(*this);
    ++tmp;
    val       = cast_to_wrapper(tmp, val).val;
#endif
    return *this;
  }

  KOKKOS_FUNCTION
  floating_point_wrapper& operator--() {
#ifdef KOKKOS_HALF_IS_FULL_TYPE_ON_ARCH
    val = val - impl_type(1.0F);  // cuda has no operator-- for __nv_bfloat
#else
    float tmp = cast_from_wrapper<float>(*this);
    --tmp;
    val = cast_to_wrapper(tmp, val).val;
#endif
    return *this;
  }

  // Postfix operators
  KOKKOS_FUNCTION
  floating_point_wrapper operator++(int) {
    floating_point_wrapper tmp = *this;
    operator++();
    return tmp;
  }

  KOKKOS_FUNCTION
  floating_point_wrapper operator--(int) {
    floating_point_wrapper tmp = *this;
    operator--();
    return tmp;
  }

  // Binary operators
  KOKKOS_FUNCTION
  floating_point_wrapper& operator=(impl_type rhs) {
    val = rhs;
    return *this;
  }

  template <class T>
  KOKKOS_FUNCTION floating_point_wrapper& operator=(T rhs) {
    val = cast_to_wrapper(rhs, val).val;
    return *this;
  }

  template <class T>
  KOKKOS_FUNCTION void operator=(T rhs) volatile {
    impl_type new_val = cast_to_wrapper(rhs, val).val;
    volatile fixed_width_integer_type* val_ptr =
        reinterpret_cast<volatile fixed_width_integer_type*>(
            const_cast<impl_type*>(&val));
    *val_ptr = reinterpret_cast<fixed_width_integer_type&>(new_val);
  }

  // Compound operators
  KOKKOS_FUNCTION
  floating_point_wrapper& operator+=(floating_point_wrapper rhs) {
#ifdef KOKKOS_HALF_IS_FULL_TYPE_ON_ARCH
    val = val + rhs.val;  // cuda has no operator+= for __nv_bfloat
#else
    val = cast_to_wrapper(
              cast_from_wrapper<float>(*this) + cast_from_wrapper<float>(rhs),
              val)
              .val;
#endif
    return *this;
  }

  KOKKOS_FUNCTION
  void operator+=(const volatile floating_point_wrapper& rhs) volatile {
    floating_point_wrapper tmp_rhs = rhs;
    floating_point_wrapper tmp_lhs = *this;

    tmp_lhs += tmp_rhs;
    *this = tmp_lhs;
  }

  // Compound operators: upcast overloads for +=
  template <class T>
  KOKKOS_FUNCTION friend std::enable_if_t<
      std::is_same<T, float>::value || std::is_same<T, double>::value, T>
  operator+=(T& lhs, floating_point_wrapper rhs) {
    lhs += static_cast<T>(rhs);
    return lhs;
  }

  KOKKOS_FUNCTION
  floating_point_wrapper& operator+=(float rhs) {
    float result = static_cast<float>(val) + rhs;
    val          = static_cast<impl_type>(result);
    return *this;
  }

  KOKKOS_FUNCTION
  floating_point_wrapper& operator+=(double rhs) {
    double result = static_cast<double>(val) + rhs;
    val           = static_cast<impl_type>(result);
    return *this;
  }

  KOKKOS_FUNCTION
  floating_point_wrapper& operator-=(floating_point_wrapper rhs) {
#ifdef KOKKOS_HALF_IS_FULL_TYPE_ON_ARCH
    val = val - rhs.val;  // cuda has no operator-= for __nv_bfloat
#else
    val = cast_to_wrapper(
              cast_from_wrapper<float>(*this) - cast_from_wrapper<float>(rhs),
              val)
              .val;
#endif
    return *this;
  }

  KOKKOS_FUNCTION
  void operator-=(const volatile floating_point_wrapper& rhs) volatile {
    floating_point_wrapper tmp_rhs = rhs;
    floating_point_wrapper tmp_lhs = *this;

    tmp_lhs -= tmp_rhs;
    *this = tmp_lhs;
  }

  // Compund operators: upcast overloads for -=
  template <class T>
  KOKKOS_FUNCTION friend std::enable_if_t<
      std::is_same<T, float>::value || std::is_same<T, double>::value, T>
  operator-=(T& lhs, floating_point_wrapper rhs) {
    lhs -= static_cast<T>(rhs);
    return lhs;
  }

  KOKKOS_FUNCTION
  floating_point_wrapper& operator-=(float rhs) {
    float result = static_cast<float>(val) - rhs;
    val          = static_cast<impl_type>(result);
    return *this;
  }

  KOKKOS_FUNCTION
  floating_point_wrapper& operator-=(double rhs) {
    double result = static_cast<double>(val) - rhs;
    val           = static_cast<impl_type>(result);
    return *this;
  }

  KOKKOS_FUNCTION
  floating_point_wrapper& operator*=(floating_point_wrapper rhs) {
#ifdef KOKKOS_HALF_IS_FULL_TYPE_ON_ARCH
    val = val * rhs.val;  // cuda has no operator*= for __nv_bfloat
#else
    val = cast_to_wrapper(
              cast_from_wrapper<float>(*this) * cast_from_wrapper<float>(rhs),
              val)
              .val;
#endif
    return *this;
  }

  KOKKOS_FUNCTION
  void operator*=(const volatile floating_point_wrapper& rhs) volatile {
    floating_point_wrapper tmp_rhs = rhs;
    floating_point_wrapper tmp_lhs = *this;

    tmp_lhs *= tmp_rhs;
    *this = tmp_lhs;
  }

  // Compund operators: upcast overloads for *=
  template <class T>
  KOKKOS_FUNCTION friend std::enable_if_t<
      std::is_same<T, float>::value || std::is_same<T, double>::value, T>
  operator*=(T& lhs, floating_point_wrapper rhs) {
    lhs *= static_cast<T>(rhs);
    return lhs;
  }

  KOKKOS_FUNCTION
  floating_point_wrapper& operator*=(float rhs) {
    float result = static_cast<float>(val) * rhs;
    val          = static_cast<impl_type>(result);
    return *this;
  }

  KOKKOS_FUNCTION
  floating_point_wrapper& operator*=(double rhs) {
    double result = static_cast<double>(val) * rhs;
    val           = static_cast<impl_type>(result);
    return *this;
  }

  KOKKOS_FUNCTION
  floating_point_wrapper& operator/=(floating_point_wrapper rhs) {
#ifdef KOKKOS_HALF_IS_FULL_TYPE_ON_ARCH
    val = val / rhs.val;  // cuda has no operator/= for __nv_bfloat
#else
    val = cast_to_wrapper(
              cast_from_wrapper<float>(*this) / cast_from_wrapper<float>(rhs),
              val)
              .val;
#endif
    return *this;
  }

  KOKKOS_FUNCTION
  void operator/=(const volatile floating_point_wrapper& rhs) volatile {
    floating_point_wrapper tmp_rhs = rhs;
    floating_point_wrapper tmp_lhs = *this;

    tmp_lhs /= tmp_rhs;
    *this = tmp_lhs;
  }

  // Compund operators: upcast overloads for /=
  template <class T>
  KOKKOS_FUNCTION friend std::enable_if_t<
      std::is_same<T, float>::value || std::is_same<T, double>::value, T>
  operator/=(T& lhs, floating_point_wrapper rhs) {
    lhs /= static_cast<T>(rhs);
    return lhs;
  }

  KOKKOS_FUNCTION
  floating_point_wrapper& operator/=(float rhs) {
    float result = static_cast<float>(val) / rhs;
    val          = static_cast<impl_type>(result);
    return *this;
  }

  KOKKOS_FUNCTION
  floating_point_wrapper& operator/=(double rhs) {
    double result = static_cast<double>(val) / rhs;
    val           = static_cast<impl_type>(result);
    return *this;
  }

  // Binary Arithmetic
  KOKKOS_FUNCTION
  friend floating_point_wrapper operator+(floating_point_wrapper lhs,
                                          floating_point_wrapper rhs) {
#ifdef KOKKOS_HALF_IS_FULL_TYPE_ON_ARCH
    lhs += rhs;
#else
    lhs.val = cast_to_wrapper(
                  cast_from_wrapper<float>(lhs) + cast_from_wrapper<float>(rhs),
                  lhs.val)
                  .val;
#endif
    return lhs;
  }

  // Binary Arithmetic upcast operators for +
  template <class T>
  KOKKOS_FUNCTION friend std::enable_if_t<
      std::is_same<T, float>::value || std::is_same<T, double>::value, T>
  operator+(floating_point_wrapper lhs, T rhs) {
    return T(lhs) + rhs;
  }

  template <class T>
  KOKKOS_FUNCTION friend std::enable_if_t<
      std::is_same<T, float>::value || std::is_same<T, double>::value, T>
  operator+(T lhs, floating_point_wrapper rhs) {
    return lhs + T(rhs);
  }

  KOKKOS_FUNCTION
  friend floating_point_wrapper operator-(floating_point_wrapper lhs,
                                          floating_point_wrapper rhs) {
#ifdef KOKKOS_HALF_IS_FULL_TYPE_ON_ARCH
    lhs -= rhs;
#else
    lhs.val = cast_to_wrapper(
                  cast_from_wrapper<float>(lhs) - cast_from_wrapper<float>(rhs),
                  lhs.val)
                  .val;
#endif
    return lhs;
  }

  // Binary Arithmetic upcast operators for -
  template <class T>
  KOKKOS_FUNCTION friend std::enable_if_t<
      std::is_same<T, float>::value || std::is_same<T, double>::value, T>
  operator-(floating_point_wrapper lhs, T rhs) {
    return T(lhs) - rhs;
  }

  template <class T>
  KOKKOS_FUNCTION friend std::enable_if_t<
      std::is_same<T, float>::value || std::is_same<T, double>::value, T>
  operator-(T lhs, floating_point_wrapper rhs) {
    return lhs - T(rhs);
  }

  KOKKOS_FUNCTION
  friend floating_point_wrapper operator*(floating_point_wrapper lhs,
                                          floating_point_wrapper rhs) {
#ifdef KOKKOS_HALF_IS_FULL_TYPE_ON_ARCH
    lhs *= rhs;
#else
    lhs.val = cast_to_wrapper(
                  cast_from_wrapper<float>(lhs) * cast_from_wrapper<float>(rhs),
                  lhs.val)
                  .val;
#endif
    return lhs;
  }

  // Binary Arithmetic upcast operators for *
  template <class T>
  KOKKOS_FUNCTION friend std::enable_if_t<
      std::is_same<T, float>::value || std::is_same<T, double>::value, T>
  operator*(floating_point_wrapper lhs, T rhs) {
    return T(lhs) * rhs;
  }

  template <class T>
  KOKKOS_FUNCTION friend std::enable_if_t<
      std::is_same<T, float>::value || std::is_same<T, double>::value, T>
  operator*(T lhs, floating_point_wrapper rhs) {
    return lhs * T(rhs);
  }

  KOKKOS_FUNCTION
  friend floating_point_wrapper operator/(floating_point_wrapper lhs,
                                          floating_point_wrapper rhs) {
#ifdef KOKKOS_HALF_IS_FULL_TYPE_ON_ARCH
    lhs /= rhs;
#else
    lhs.val = cast_to_wrapper(
                  cast_from_wrapper<float>(lhs) / cast_from_wrapper<float>(rhs),
                  lhs.val)
                  .val;
#endif
    return lhs;
  }

  // Binary Arithmetic upcast operators for /
  template <class T>
  KOKKOS_FUNCTION friend std::enable_if_t<
      std::is_same<T, float>::value || std::is_same<T, double>::value, T>
  operator/(floating_point_wrapper lhs, T rhs) {
    return T(lhs) / rhs;
  }

  template <class T>
  KOKKOS_FUNCTION friend std::enable_if_t<
      std::is_same<T, float>::value || std::is_same<T, double>::value, T>
  operator/(T lhs, floating_point_wrapper rhs) {
    return lhs / T(rhs);
  }

  // Logical operators
  KOKKOS_FUNCTION
  bool operator!() const {
#ifdef KOKKOS_HALF_IS_FULL_TYPE_ON_ARCH
    return static_cast<bool>(!val);
#else
    return !cast_from_wrapper<float>(*this);
#endif
  }

  // NOTE: Loses short-circuit evaluation
  KOKKOS_FUNCTION
  bool operator&&(floating_point_wrapper rhs) const {
#ifdef KOKKOS_HALF_IS_FULL_TYPE_ON_ARCH
    return static_cast<bool>(val && rhs.val);
#else
    return cast_from_wrapper<float>(*this) && cast_from_wrapper<float>(rhs);
#endif
  }

  // NOTE: Loses short-circuit evaluation
  KOKKOS_FUNCTION
  bool operator||(floating_point_wrapper rhs) const {
#ifdef KOKKOS_HALF_IS_FULL_TYPE_ON_ARCH
    return static_cast<bool>(val || rhs.val);
#else
    return cast_from_wrapper<float>(*this) || cast_from_wrapper<float>(rhs);
#endif
  }

  // Comparison operators
  KOKKOS_FUNCTION
  bool operator==(floating_point_wrapper rhs) const {
#ifdef KOKKOS_HALF_IS_FULL_TYPE_ON_ARCH
    return static_cast<bool>(val == rhs.val);
#else
    return cast_from_wrapper<float>(*this) == cast_from_wrapper<float>(rhs);
#endif
  }

  KOKKOS_FUNCTION
  bool operator!=(floating_point_wrapper rhs) const {
#ifdef KOKKOS_HALF_IS_FULL_TYPE_ON_ARCH
    return static_cast<bool>(val != rhs.val);
#else
    return cast_from_wrapper<float>(*this) != cast_from_wrapper<float>(rhs);
#endif
  }

  KOKKOS_FUNCTION
  bool operator<(floating_point_wrapper rhs) const {
#ifdef KOKKOS_HALF_IS_FULL_TYPE_ON_ARCH
    return static_cast<bool>(val < rhs.val);
#else
    return cast_from_wrapper<float>(*this) < cast_from_wrapper<float>(rhs);
#endif
  }

  KOKKOS_FUNCTION
  bool operator>(floating_point_wrapper rhs) const {
#ifdef KOKKOS_HALF_IS_FULL_TYPE_ON_ARCH
    return static_cast<bool>(val > rhs.val);
#else
    return cast_from_wrapper<float>(*this) > cast_from_wrapper<float>(rhs);
#endif
  }

  KOKKOS_FUNCTION
  bool operator<=(floating_point_wrapper rhs) const {
#ifdef KOKKOS_HALF_IS_FULL_TYPE_ON_ARCH
    return static_cast<bool>(val <= rhs.val);
#else
    return cast_from_wrapper<float>(*this) <= cast_from_wrapper<float>(rhs);
#endif
  }

  KOKKOS_FUNCTION
  bool operator>=(floating_point_wrapper rhs) const {
#ifdef KOKKOS_HALF_IS_FULL_TYPE_ON_ARCH
    return static_cast<bool>(val >= rhs.val);
#else
    return cast_from_wrapper<float>(*this) >= cast_from_wrapper<float>(rhs);
#endif
  }

  KOKKOS_FUNCTION
  friend bool operator==(const volatile floating_point_wrapper& lhs,
                         const volatile floating_point_wrapper& rhs) {
    floating_point_wrapper tmp_lhs = lhs, tmp_rhs = rhs;
    return tmp_lhs == tmp_rhs;
  }

  KOKKOS_FUNCTION
  friend bool operator!=(const volatile floating_point_wrapper& lhs,
                         const volatile floating_point_wrapper& rhs) {
    floating_point_wrapper tmp_lhs = lhs, tmp_rhs = rhs;
    return tmp_lhs != tmp_rhs;
  }

  KOKKOS_FUNCTION
  friend bool operator<(const volatile floating_point_wrapper& lhs,
                        const volatile floating_point_wrapper& rhs) {
    floating_point_wrapper tmp_lhs = lhs, tmp_rhs = rhs;
    return tmp_lhs < tmp_rhs;
  }

  KOKKOS_FUNCTION
  friend bool operator>(const volatile floating_point_wrapper& lhs,
                        const volatile floating_point_wrapper& rhs) {
    floating_point_wrapper tmp_lhs = lhs, tmp_rhs = rhs;
    return tmp_lhs > tmp_rhs;
  }

  KOKKOS_FUNCTION
  friend bool operator<=(const volatile floating_point_wrapper& lhs,
                         const volatile floating_point_wrapper& rhs) {
    floating_point_wrapper tmp_lhs = lhs, tmp_rhs = rhs;
    return tmp_lhs <= tmp_rhs;
  }

  KOKKOS_FUNCTION
  friend bool operator>=(const volatile floating_point_wrapper& lhs,
                         const volatile floating_point_wrapper& rhs) {
    floating_point_wrapper tmp_lhs = lhs, tmp_rhs = rhs;
    return tmp_lhs >= tmp_rhs;
  }

  // Insertion and extraction operators
  friend std::ostream& operator<<(std::ostream& os,
                                  const floating_point_wrapper& x) {
    const std::string out = std::to_string(static_cast<double>(x));
    os << out;
    return os;
  }

  friend std::istream& operator>>(std::istream& is, floating_point_wrapper& x) {
    std::string in;
    is >> in;
    x = std::stod(in);
    return is;
  }
};
}  // namespace Impl

// Declare wrapper overloads now that floating_point_wrapper is declared
template <class T>
static KOKKOS_INLINE_FUNCTION Kokkos::Experimental::half_t cast_to_wrapper(
    T x, const volatile Kokkos::Impl::half_impl_t::type&) {
  return Kokkos::Experimental::cast_to_half(x);
}

#ifdef KOKKOS_IMPL_BHALF_TYPE_DEFINED
template <class T>
static KOKKOS_INLINE_FUNCTION Kokkos::Experimental::bhalf_t cast_to_wrapper(
    T x, const volatile Kokkos::Impl::bhalf_impl_t::type&) {
  return Kokkos::Experimental::cast_to_bhalf(x);
}
#endif  // KOKKOS_IMPL_BHALF_TYPE_DEFINED

template <class T>
static KOKKOS_INLINE_FUNCTION T
cast_from_wrapper(const Kokkos::Experimental::half_t& x) {
  return Kokkos::Experimental::cast_from_half<T>(x);
}

#ifdef KOKKOS_IMPL_BHALF_TYPE_DEFINED
template <class T>
static KOKKOS_INLINE_FUNCTION T
cast_from_wrapper(const Kokkos::Experimental::bhalf_t& x) {
  return Kokkos::Experimental::cast_from_bhalf<T>(x);
}
#endif  // KOKKOS_IMPL_BHALF_TYPE_DEFINED

}  // namespace Experimental
}  // namespace Kokkos

#endif  // KOKKOS_IMPL_HALF_TYPE_DEFINED

// If none of the above actually did anything and defined a half precision type
// define a fallback implementation here using float
#ifndef KOKKOS_IMPL_HALF_TYPE_DEFINED
#define KOKKOS_IMPL_HALF_TYPE_DEFINED
#define KOKKOS_HALF_T_IS_FLOAT true
namespace Kokkos {
namespace Impl {
struct half_impl_t {
  using type = float;
};
}  // namespace Impl
namespace Experimental {

using half_t = Kokkos::Impl::half_impl_t::type;

// cast_to_half
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(float val) { return half_t(val); }
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(bool val) { return half_t(val); }
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(double val) { return half_t(val); }
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(short val) { return half_t(val); }
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(unsigned short val) { return half_t(val); }
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(int val) { return half_t(val); }
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(unsigned int val) { return half_t(val); }
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(long val) { return half_t(val); }
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(unsigned long val) { return half_t(val); }
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(long long val) { return half_t(val); }
KOKKOS_INLINE_FUNCTION
half_t cast_to_half(unsigned long long val) { return half_t(val); }

// cast_from_half
// Using an explicit list here too, since the other ones are explicit and for
// example don't include char
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<
    std::is_same<T, float>::value || std::is_same<T, bool>::value ||
        std::is_same<T, double>::value || std::is_same<T, short>::value ||
        std::is_same<T, unsigned short>::value || std::is_same<T, int>::value ||
        std::is_same<T, unsigned int>::value || std::is_same<T, long>::value ||
        std::is_same<T, unsigned long>::value ||
        std::is_same<T, long long>::value ||
        std::is_same<T, unsigned long long>::value,
    T>
cast_from_half(half_t val) {
  return T(val);
}

}  // namespace Experimental
}  // namespace Kokkos

#else
#define KOKKOS_HALF_T_IS_FLOAT false
#endif  // KOKKOS_IMPL_HALF_TYPE_DEFINED

#ifndef KOKKOS_IMPL_BHALF_TYPE_DEFINED
#define KOKKOS_IMPL_BHALF_TYPE_DEFINED
#define KOKKOS_BHALF_T_IS_FLOAT true
namespace Kokkos {
namespace Impl {
struct bhalf_impl_t {
  using type = float;
};
}  // namespace Impl

namespace Experimental {

using bhalf_t = Kokkos::Impl::bhalf_impl_t::type;

// cast_to_bhalf
KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(float val) { return bhalf_t(val); }
KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(bool val) { return bhalf_t(val); }
KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(double val) { return bhalf_t(val); }
KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(short val) { return bhalf_t(val); }
KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(unsigned short val) { return bhalf_t(val); }
KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(int val) { return bhalf_t(val); }
KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(unsigned int val) { return bhalf_t(val); }
KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(long val) { return bhalf_t(val); }
KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(unsigned long val) { return bhalf_t(val); }
KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(long long val) { return bhalf_t(val); }
KOKKOS_INLINE_FUNCTION
bhalf_t cast_to_bhalf(unsigned long long val) { return bhalf_t(val); }

// cast_from_bhalf
template <class T>
KOKKOS_INLINE_FUNCTION std::enable_if_t<
    std::is_same<T, float>::value || std::is_same<T, bool>::value ||
        std::is_same<T, double>::value || std::is_same<T, short>::value ||
        std::is_same<T, unsigned short>::value || std::is_same<T, int>::value ||
        std::is_same<T, unsigned int>::value || std::is_same<T, long>::value ||
        std::is_same<T, unsigned long>::value ||
        std::is_same<T, long long>::value ||
        std::is_same<T, unsigned long long>::value,
    T>
cast_from_bhalf(bhalf_t val) {
  return T(val);
}
}  // namespace Experimental
}  // namespace Kokkos
#else
#define KOKKOS_BHALF_T_IS_FLOAT false
#endif  // KOKKOS_IMPL_BHALF_TYPE_DEFINED

#endif  // KOKKOS_HALF_FLOATING_POINT_WRAPPER_HPP_
