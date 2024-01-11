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

#ifndef TESTHALFCONVERSION_HPP_
#define TESTHALFCONVERSION_HPP_
namespace Test {

template <class T>
void test_half_conversion_type() {
  double epsilon                 = KOKKOS_HALF_T_IS_FLOAT ? 0.0000003 : 0.0003;
  T base                         = static_cast<T>(3.3);
  Kokkos::Experimental::half_t a = Kokkos::Experimental::cast_to_half(base);
  T b                            = Kokkos::Experimental::cast_from_half<T>(a);
  ASSERT_LT((double(b - base) / double(base)), epsilon);

#ifdef KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA
  Kokkos::View<T> b_v("b_v");
  Kokkos::parallel_for(
      "TestHalfConversion", 1, KOKKOS_LAMBDA(int) {
        Kokkos::Experimental::half_t d_a =
            Kokkos::Experimental::cast_to_half(base);
        b_v() = Kokkos::Experimental::cast_from_half<T>(d_a);
      });

  Kokkos::deep_copy(b, b_v);
  ASSERT_LT((double(b - base) / double(base)), epsilon);
#endif  // KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA
}

template <class T>
void test_bhalf_conversion_type() {
  double epsilon = KOKKOS_BHALF_T_IS_FLOAT ? 0.0000003 : 0.0003;
  T base         = static_cast<T>(3.3);
  Kokkos::Experimental::bhalf_t a = Kokkos::Experimental::cast_to_bhalf(base);
  T b                             = Kokkos::Experimental::cast_from_bhalf<T>(a);
  ASSERT_LT((double(b - base) / double(base)), epsilon);

#ifdef KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA
  Kokkos::View<T> b_v("b_v");
  Kokkos::parallel_for(
      "TestHalfConversion", 1, KOKKOS_LAMBDA(int) {
        Kokkos::Experimental::bhalf_t d_a =
            Kokkos::Experimental::cast_to_bhalf(base);
        b_v() = Kokkos::Experimental::cast_from_bhalf<T>(d_a);
      });

  Kokkos::deep_copy(b, b_v);
  ASSERT_LT((double(b - base) / double(base)), epsilon);
#endif  // KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA
}

void test_half_conversion() {
  test_half_conversion_type<float>();
  test_half_conversion_type<double>();
  test_half_conversion_type<short>();
  test_half_conversion_type<int>();
  test_half_conversion_type<long>();
  test_half_conversion_type<long long>();
  test_half_conversion_type<unsigned short>();
  test_half_conversion_type<unsigned int>();
  test_half_conversion_type<unsigned long>();
  test_half_conversion_type<unsigned long long>();
}

void test_bhalf_conversion() {
  test_bhalf_conversion_type<float>();
  test_bhalf_conversion_type<double>();
  test_bhalf_conversion_type<short>();
  test_bhalf_conversion_type<int>();
  test_bhalf_conversion_type<long>();
  test_bhalf_conversion_type<long long>();
  test_bhalf_conversion_type<unsigned short>();
  test_bhalf_conversion_type<unsigned int>();
  test_bhalf_conversion_type<unsigned long>();
  test_bhalf_conversion_type<unsigned long long>();
}

TEST(TEST_CATEGORY, half_conversion) { test_half_conversion(); }

TEST(TEST_CATEGORY, bhalf_conversion) { test_bhalf_conversion(); }

}  // namespace Test
#endif
