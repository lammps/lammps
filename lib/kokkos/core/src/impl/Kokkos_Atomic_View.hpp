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
#ifndef KOKKOS_ATOMIC_VIEW_HPP
#define KOKKOS_ATOMIC_VIEW_HPP

#include <Kokkos_Macros.hpp>
#include <Kokkos_Atomic.hpp>

namespace Kokkos {
namespace Impl {

// The following tag is used to prevent an implicit call of the constructor when
// trying to assign a literal 0 int ( = 0 );
struct AtomicViewConstTag {};

template <class ViewTraits>
class AtomicDataElement {
 public:
  using value_type           = typename ViewTraits::value_type;
  using const_value_type     = typename ViewTraits::const_value_type;
  using non_const_value_type = typename ViewTraits::non_const_value_type;
  value_type* const ptr;

  KOKKOS_INLINE_FUNCTION
  AtomicDataElement(value_type* ptr_, AtomicViewConstTag) : ptr(ptr_) {}

  KOKKOS_INLINE_FUNCTION
  const_value_type operator=(const_value_type& val) const {
    Kokkos::atomic_store(ptr, val);
    return val;
  }

  KOKKOS_INLINE_FUNCTION
  void inc() const { Kokkos::atomic_increment(ptr); }

  KOKKOS_INLINE_FUNCTION
  void dec() const { Kokkos::atomic_decrement(ptr); }

  KOKKOS_INLINE_FUNCTION
  const_value_type operator++() const {
    const_value_type tmp =
        Kokkos::atomic_fetch_add(ptr, non_const_value_type(1));
    return tmp + 1;
  }

  KOKKOS_INLINE_FUNCTION
  const_value_type operator--() const {
    const_value_type tmp =
        Kokkos::atomic_fetch_sub(ptr, non_const_value_type(1));
    return tmp - 1;
  }

  KOKKOS_INLINE_FUNCTION
  const_value_type operator++(int) const {
    return Kokkos::atomic_fetch_add(ptr, non_const_value_type(1));
  }

  KOKKOS_INLINE_FUNCTION
  const_value_type operator--(int) const {
    return Kokkos::atomic_fetch_sub(ptr, non_const_value_type(1));
  }

  KOKKOS_INLINE_FUNCTION
  const_value_type operator+=(const_value_type& val) const {
    const_value_type tmp = Kokkos::atomic_fetch_add(ptr, val);
    return tmp + val;
  }

  KOKKOS_INLINE_FUNCTION
  const_value_type operator-=(const_value_type& val) const {
    const_value_type tmp = Kokkos::atomic_fetch_sub(ptr, val);
    return tmp - val;
  }

  KOKKOS_INLINE_FUNCTION
  const_value_type operator*=(const_value_type& val) const {
    return Kokkos::atomic_mul_fetch(ptr, val);
  }

  KOKKOS_INLINE_FUNCTION
  const_value_type operator/=(const_value_type& val) const {
    return Kokkos::atomic_div_fetch(ptr, val);
  }

  KOKKOS_INLINE_FUNCTION
  const_value_type operator%=(const_value_type& val) const {
    return Kokkos::atomic_mod_fetch(ptr, val);
  }

  KOKKOS_INLINE_FUNCTION
  const_value_type operator&=(const_value_type& val) const {
    return Kokkos::atomic_and_fetch(ptr, val);
  }

  KOKKOS_INLINE_FUNCTION
  const_value_type operator^=(const_value_type& val) const {
    return Kokkos::atomic_xor_fetch(ptr, val);
  }

  KOKKOS_INLINE_FUNCTION
  const_value_type operator|=(const_value_type& val) const {
    return Kokkos::atomic_or_fetch(ptr, val);
  }

  KOKKOS_INLINE_FUNCTION
  const_value_type operator<<=(const_value_type& val) const {
    return Kokkos::atomic_lshift_fetch(ptr, val);
  }

  KOKKOS_INLINE_FUNCTION
  const_value_type operator>>=(const_value_type& val) const {
    return Kokkos::atomic_rshift_fetch(ptr, val);
  }

  KOKKOS_INLINE_FUNCTION
  const_value_type operator+(const_value_type& val) const { return *ptr + val; }

  KOKKOS_INLINE_FUNCTION
  const_value_type operator-(const_value_type& val) const { return *ptr - val; }

  KOKKOS_INLINE_FUNCTION
  const_value_type operator*(const_value_type& val) const { return *ptr * val; }

  KOKKOS_INLINE_FUNCTION
  const_value_type operator/(const_value_type& val) const { return *ptr / val; }

  KOKKOS_INLINE_FUNCTION
  const_value_type operator%(const_value_type& val) const { return *ptr ^ val; }

  KOKKOS_INLINE_FUNCTION
  const_value_type operator!() const { return !*ptr; }

  KOKKOS_INLINE_FUNCTION
  const_value_type operator&&(const_value_type& val) const {
    return *ptr && val;
  }

  KOKKOS_INLINE_FUNCTION
  const_value_type operator||(const_value_type& val) const {
    return *ptr | val;
  }

  KOKKOS_INLINE_FUNCTION
  const_value_type operator&(const_value_type& val) const { return *ptr & val; }

  KOKKOS_INLINE_FUNCTION
  const_value_type operator|(const_value_type& val) const { return *ptr | val; }

  KOKKOS_INLINE_FUNCTION
  const_value_type operator^(const_value_type& val) const { return *ptr ^ val; }

  KOKKOS_INLINE_FUNCTION
  const_value_type operator~() const { return ~*ptr; }

  KOKKOS_INLINE_FUNCTION
  const_value_type operator<<(const unsigned int& val) const {
    return *ptr << val;
  }

  KOKKOS_INLINE_FUNCTION
  const_value_type operator>>(const unsigned int& val) const {
    return *ptr >> val;
  }

  KOKKOS_INLINE_FUNCTION
  bool operator==(const AtomicDataElement& val) const { return *ptr == val; }

  KOKKOS_INLINE_FUNCTION
  bool operator!=(const AtomicDataElement& val) const { return *ptr != val; }

  KOKKOS_INLINE_FUNCTION
  bool operator>=(const_value_type& val) const { return *ptr >= val; }

  KOKKOS_INLINE_FUNCTION
  bool operator<=(const_value_type& val) const { return *ptr <= val; }

  KOKKOS_INLINE_FUNCTION
  bool operator<(const_value_type& val) const { return *ptr < val; }

  KOKKOS_INLINE_FUNCTION
  bool operator>(const_value_type& val) const { return *ptr > val; }

  KOKKOS_INLINE_FUNCTION
  operator value_type() const { return Kokkos::atomic_load(ptr); }
};

template <class ViewTraits>
class AtomicViewDataHandle {
 public:
  typename ViewTraits::value_type* ptr;

  KOKKOS_INLINE_FUNCTION
  AtomicViewDataHandle() : ptr(nullptr) {}

  KOKKOS_INLINE_FUNCTION
  AtomicViewDataHandle(typename ViewTraits::value_type* ptr_) : ptr(ptr_) {}

  template <class iType>
  KOKKOS_INLINE_FUNCTION AtomicDataElement<ViewTraits> operator[](
      const iType& i) const {
    return AtomicDataElement<ViewTraits>(ptr + i, AtomicViewConstTag());
  }

  KOKKOS_INLINE_FUNCTION
  operator typename ViewTraits::value_type*() const { return ptr; }
};

}  // namespace Impl
}  // namespace Kokkos

#endif
