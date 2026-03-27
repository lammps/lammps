// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_SIMD_TESTING_UTILITIES_HPP
#define KOKKOS_SIMD_TESTING_UTILITIES_HPP

#include <gtest/gtest.h>
#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.simd;
#else
#include <Kokkos_SIMD.hpp>
#endif
#include <SIMDTesting_Ops.hpp>

class gtest_checker {
 public:
  void truth(bool x) const { EXPECT_TRUE(x); }
  template <class T>
  void equality(T const& a, T const& b) const {
    if constexpr (std::is_same_v<T, double>) {
      EXPECT_DOUBLE_EQ(a, b);
    } else if constexpr (std::is_same_v<T, float>) {
      EXPECT_FLOAT_EQ(a, b);
    } else {
      EXPECT_EQ(a, b);
    }
  }
};

class kokkos_checker {
 public:
  KOKKOS_INLINE_FUNCTION void truth(bool x) const {
    if (!x) Kokkos::abort("SIMD unit test truth condition failed on device");
  }
  template <class T>
  KOKKOS_INLINE_FUNCTION void equality(T const& a, T const& b) const {
#if defined(KOKKOS_IMPL_32BIT)
    // This is needed to work around a bug where the comparison fails because it
    // is done on the x87 fpu (which is the default for 32 bit gcc) in long
    // double and a and b end up being different in long double but have the
    // same value when casted to float or double. (see
    // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=323#c109)
    T const volatile* const volatile pa = &a;
    T const volatile* const volatile pb = &b;
    if (*pa != *pb)
      Kokkos::abort("SIMD unit test equality condition failed on device");
#else
    if (a != b)
      Kokkos::abort("SIMD unit test equality condition failed on device");
#endif
  }
};

template <class T, class Abi>
inline void host_check_equality(
    Kokkos::Experimental::basic_simd<T, Abi> const& expected_result,
    Kokkos::Experimental::basic_simd<T, Abi> const& computed_result,
    std::size_t nlanes) {
  gtest_checker checker;
  for (std::size_t i = 0; i < nlanes; ++i) {
    checker.equality(expected_result[i], computed_result[i]);
  }
}

template <class T, class Abi>
KOKKOS_INLINE_FUNCTION void device_check_equality(
    Kokkos::Experimental::basic_simd<T, Abi> const& expected_result,
    Kokkos::Experimental::basic_simd<T, Abi> const& computed_result,
    std::size_t nlanes) {
  kokkos_checker checker;
  for (std::size_t i = 0; i < nlanes; ++i) {
    checker.equality(expected_result[i], computed_result[i]);
  }
}

template <typename T, typename Abi>
KOKKOS_INLINE_FUNCTION void check_equality(
    Kokkos::Experimental::basic_simd<T, Abi> const& expected_result,
    Kokkos::Experimental::basic_simd<T, Abi> const& computed_result,
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
                 Kokkos::Experimental::basic_simd<T, Abi>& result) const {
    if (n < result.size()) return false;
    using simd_type = Kokkos::Experimental::basic_simd<T, Abi>;
    result          = Kokkos::Experimental::simd_unchecked_load<simd_type>(
        mem, Kokkos::Experimental::simd_flag_default);
    return true;
  }
  template <class T, class Abi>
  KOKKOS_INLINE_FUNCTION bool device_load(
      T const* mem, std::size_t n,
      Kokkos::Experimental::basic_simd<T, Abi>& result) const {
    if (n < result.size()) return false;
    using simd_type = Kokkos::Experimental::basic_simd<T, Abi>;
    result          = Kokkos::Experimental::simd_unchecked_load<simd_type>(
        mem, Kokkos::Experimental::simd_flag_default);
    return true;
  }
};

class load_vector_aligned {
 public:
  template <class T, class Abi>
  bool host_load(T const* mem, std::size_t n,
                 Kokkos::Experimental::basic_simd<T, Abi>& result) const {
    if (n < result.size()) return false;
    using simd_type = Kokkos::Experimental::basic_simd<T, Abi>;
    result          = Kokkos::Experimental::simd_unchecked_load<simd_type>(
        mem, Kokkos::Experimental::simd_flag_aligned);
    return true;
  }
  template <class T, class Abi>
  KOKKOS_INLINE_FUNCTION bool device_load(
      T const* mem, std::size_t n,
      Kokkos::Experimental::basic_simd<T, Abi>& result) const {
    if (n < result.size()) return false;
    using simd_type = Kokkos::Experimental::basic_simd<T, Abi>;
    result          = Kokkos::Experimental::simd_unchecked_load<simd_type>(
        mem, Kokkos::Experimental::simd_flag_aligned);
    return true;
  }
};

class load_masked {
 public:
  template <class T, class Abi>
  bool host_load(T const* mem, std::size_t n,
                 Kokkos::Experimental::basic_simd<T, Abi>& result) const {
    using mask_type =
        typename Kokkos::Experimental::basic_simd<T, Abi>::mask_type;
    mask_type mask(KOKKOS_LAMBDA(std::size_t i) { return i < n; });
    result =
        simd_partial_load(mem, mask, Kokkos::Experimental::simd_flag_default);
    return true;
  }
  template <class T, class Abi>
  KOKKOS_INLINE_FUNCTION bool device_load(
      T const* mem, std::size_t n,
      Kokkos::Experimental::basic_simd<T, Abi>& result) const {
    using mask_type =
        typename Kokkos::Experimental::basic_simd<T, Abi>::mask_type;
    mask_type mask(KOKKOS_LAMBDA(std::size_t i) { return i < n; });
    result =
        simd_partial_load(mem, mask, Kokkos::Experimental::simd_flag_default);
    return true;
  }
};

class load_as_scalars {
 public:
  template <class T, class Abi>
  bool host_load(T const* mem, std::size_t n,
                 Kokkos::Experimental::basic_simd<T, Abi>& result) const {
    Kokkos::Experimental::basic_simd<T, Abi> init(
        KOKKOS_LAMBDA(std::size_t i) { return (i < n) ? mem[i] : T(0); });
    result = init;

    return true;
  }
  template <class T, class Abi>
  KOKKOS_INLINE_FUNCTION bool device_load(
      T const* mem, std::size_t n,
      Kokkos::Experimental::basic_simd<T, Abi>& result) const {
    Kokkos::Experimental::basic_simd<T, Abi> init(
        KOKKOS_LAMBDA(std::size_t i) { return (i < n) ? mem[i] : T(0); });

    result = init;
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

// We consider a fully-implemented 'simd' type is always accompanied by the
// same-capability 'simd_mask' type
template <typename DataType, typename Abi>
constexpr bool is_simd_avail_v =
    is_type_v<Kokkos::Experimental::basic_simd<DataType, Abi>> &&
    is_type_v<Kokkos::Experimental::basic_simd_mask<DataType, Abi>>;

#endif
