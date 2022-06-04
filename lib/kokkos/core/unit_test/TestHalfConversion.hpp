
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
