// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef TESTHALFCONVERSION_HPP_
#define TESTHALFCONVERSION_HPP_

#include <impl/Kokkos_Half_FloatingPointWrapper.hpp>

namespace Test {

template <class T>
void test_half_conversion_type() {
  // When truncating mantissa to 10bits (like f16), 3.3 becomes 3.298828125
  // 3.3 - 3.298828125 < 1.1719e-3, so conversion error should be smaller
  double epsilon                 = KOKKOS_HALF_T_IS_FLOAT ? 3e-7 : 1.1719e-3;
  T base                         = static_cast<T>(3.3);
  Kokkos::Experimental::half_t a = Kokkos::Experimental::cast_to_half(base);
  T b                            = Kokkos::Experimental::cast_from_half<T>(a);
  ASSERT_NEAR(b, base, epsilon);

  Kokkos::View<T> b_v("b_v");
  Kokkos::parallel_for(
      "TestHalfConversion", 1, KOKKOS_LAMBDA(int) {
        Kokkos::Experimental::half_t d_a =
            Kokkos::Experimental::cast_to_half(base);
        b_v() = Kokkos::Experimental::cast_from_half<T>(d_a);
      });

  Kokkos::deep_copy(b, b_v);
  ASSERT_NEAR(b, base, epsilon);
}

template <class T>
void test_bhalf_conversion_type() {
  // When truncating mantissa to 7bits (like b16), 3.3 becomes 3.296875
  // 3.3 - 3.296875 < 3.125e-3, so conversion error should be smaller
  double epsilon                  = KOKKOS_BHALF_T_IS_FLOAT ? 3e-7 : 3.125e-3;
  T base                          = static_cast<T>(3.3);
  Kokkos::Experimental::bhalf_t a = Kokkos::Experimental::cast_to_bhalf(base);
  T b                             = Kokkos::Experimental::cast_from_bhalf<T>(a);
  ASSERT_NEAR(b, base, epsilon);

  Kokkos::View<T> b_v("b_v");
  Kokkos::parallel_for(
      "TestHalfConversion", 1, KOKKOS_LAMBDA(int) {
        Kokkos::Experimental::bhalf_t d_a =
            Kokkos::Experimental::cast_to_bhalf(base);
        b_v() = Kokkos::Experimental::cast_from_bhalf<T>(d_a);
      });

  Kokkos::deep_copy(b, b_v);
  ASSERT_NEAR(b, base, epsilon);
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
