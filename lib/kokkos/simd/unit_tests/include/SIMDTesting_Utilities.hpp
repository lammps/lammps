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

#ifndef KOKKOS_SIMD_TESTING_UTILITIES_HPP
#define KOKKOS_SIMD_TESTING_UTILITIES_HPP

#include <gtest/gtest.h>
#include <Kokkos_SIMD.hpp>
#include <SIMDTesting_Ops.hpp>

class gtest_checker {
 public:
  void truth(bool x) const { EXPECT_TRUE(x); }
  template <class T>
  void equality(T const& a, T const& b) const {
    EXPECT_EQ(a, b);
  }
};

class kokkos_checker {
 public:
  KOKKOS_INLINE_FUNCTION void truth(bool x) const {
    if (!x) Kokkos::abort("SIMD unit test truth condition failed on device");
  }
  template <class T>
  KOKKOS_INLINE_FUNCTION void equality(T const& a, T const& b) const {
    if (a != b)
      Kokkos::abort("SIMD unit test equality condition failed on device");
  }
};

template <class T, class Abi>
inline void host_check_equality(
    Kokkos::Experimental::simd<T, Abi> const& expected_result,
    Kokkos::Experimental::simd<T, Abi> const& computed_result,
    std::size_t nlanes) {
  gtest_checker checker;
  for (std::size_t i = 0; i < nlanes; ++i) {
    checker.equality(expected_result[i], computed_result[i]);
  }
  using mask_type = typename Kokkos::Experimental::simd<T, Abi>::mask_type;
  mask_type mask(false);
  for (std::size_t i = 0; i < nlanes; ++i) {
    mask[i] = true;
  }
  checker.equality((expected_result == computed_result) && mask, mask);
}

template <class T, class Abi>
KOKKOS_INLINE_FUNCTION void device_check_equality(
    Kokkos::Experimental::simd<T, Abi> const& expected_result,
    Kokkos::Experimental::simd<T, Abi> const& computed_result,
    std::size_t nlanes) {
  kokkos_checker checker;
  for (std::size_t i = 0; i < nlanes; ++i) {
    checker.equality(expected_result[i], computed_result[i]);
  }
  using mask_type = typename Kokkos::Experimental::simd<T, Abi>::mask_type;
  mask_type mask(false);
  for (std::size_t i = 0; i < nlanes; ++i) {
    mask[i] = true;
  }
  checker.equality((expected_result == computed_result) && mask, mask);
}

template <typename T, typename Abi>
KOKKOS_INLINE_FUNCTION void check_equality(
    Kokkos::Experimental::simd<T, Abi> const& expected_result,
    Kokkos::Experimental::simd<T, Abi> const& computed_result,
    std::size_t nlanes) {
  KOKKOS_IF_ON_HOST(
      (host_check_equality(expected_result, computed_result, nlanes);))
  KOKKOS_IF_ON_DEVICE(
      (device_check_equality(expected_result, computed_result, nlanes);))
}

class load_element_aligned {
 public:
  template <class T, class Abi>
  bool host_load(T const* mem, std::size_t n,
                 Kokkos::Experimental::simd<T, Abi>& result) const {
    if (n < result.size()) return false;
    result.copy_from(mem, Kokkos::Experimental::simd_flag_default);
    return true;
  }
  template <class T, class Abi>
  KOKKOS_INLINE_FUNCTION bool device_load(
      T const* mem, std::size_t n,
      Kokkos::Experimental::simd<T, Abi>& result) const {
    if (n < result.size()) return false;
    result.copy_from(mem, Kokkos::Experimental::simd_flag_default);
    return true;
  }
};

class load_vector_aligned {
 public:
  template <class T, class Abi>
  bool host_load(T const* mem, std::size_t n,
                 Kokkos::Experimental::simd<T, Abi>& result) const {
    if (n < result.size()) return false;
    result.copy_from(mem, Kokkos::Experimental::simd_flag_aligned);
    return true;
  }
  template <class T, class Abi>
  KOKKOS_INLINE_FUNCTION bool device_load(
      T const* mem, std::size_t n,
      Kokkos::Experimental::simd<T, Abi>& result) const {
    if (n < result.size()) return false;
    result.copy_from(mem, Kokkos::Experimental::simd_flag_aligned);
    return true;
  }
};

class load_masked {
 public:
  template <class T, class Abi>
  bool host_load(T const* mem, std::size_t n,
                 Kokkos::Experimental::simd<T, Abi>& result) const {
    using mask_type = typename Kokkos::Experimental::simd<T, Abi>::mask_type;
    mask_type mask(false);
    for (std::size_t i = 0; i < n; ++i) {
      mask[i] = true;
    }
    result = T(0);
    where(mask, result).copy_from(mem, Kokkos::Experimental::simd_flag_default);
    return true;
  }
  template <class T, class Abi>
  KOKKOS_INLINE_FUNCTION bool device_load(
      T const* mem, std::size_t n,
      Kokkos::Experimental::simd<T, Abi>& result) const {
    using mask_type = typename Kokkos::Experimental::simd<T, Abi>::mask_type;
    mask_type mask(false);
    for (std::size_t i = 0; i < n; ++i) {
      mask[i] = true;
    }
    where(mask, result).copy_from(mem, Kokkos::Experimental::simd_flag_default);
    where(!mask, result) = T(0);
    return true;
  }
};

class load_as_scalars {
 public:
  template <class T, class Abi>
  bool host_load(T const* mem, std::size_t n,
                 Kokkos::Experimental::simd<T, Abi>& result) const {
    for (std::size_t i = 0; i < n; ++i) {
      result[i] = mem[i];
    }
    for (std::size_t i = n; i < result.size(); ++i) {
      result[i] = T(0);
    }
    return true;
  }
  template <class T, class Abi>
  KOKKOS_INLINE_FUNCTION bool device_load(
      T const* mem, std::size_t n,
      Kokkos::Experimental::simd<T, Abi>& result) const {
    for (std::size_t i = 0; i < n; ++i) {
      result[i] = mem[i];
    }
    for (std::size_t i = n; i < result.size(); ++i) {
      result[i] = T(0);
    }
    return true;
  }
};

// Simple check to loosely test that T is a complete type.
// Some capabilities are only defined for specific data type and abi pairs (i.e.
// extended vector width); this is used to exclude pairs that
// are not defined from being tested.
template <typename T, typename = void>
constexpr bool is_type_v = false;

template <typename T>
constexpr bool is_type_v<T, decltype(void(sizeof(T)))> = true;

#endif
