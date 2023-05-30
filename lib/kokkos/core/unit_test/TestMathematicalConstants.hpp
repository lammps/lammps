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

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

template <class T>
KOKKOS_FUNCTION T *take_address_of(T &arg) {
  return &arg;
}

template <class T>
KOKKOS_FUNCTION void take_by_value(T) {}

#if defined(KOKKOS_ENABLE_CXX17)
#define DEFINE_MATH_CONSTANT_TRAIT(TRAIT)                          \
  template <class T>                                               \
  struct TRAIT {                                                   \
    static constexpr T value = Kokkos::Experimental::TRAIT##_v<T>; \
  }
#else
#define DEFINE_MATH_CONSTANT_TRAIT(TRAIT)                                    \
  template <class>                                                           \
  struct TRAIT;                                                              \
  template <>                                                                \
  struct TRAIT<float> {                                                      \
    static constexpr float value = Kokkos::Experimental::TRAIT##_v<float>;   \
  };                                                                         \
  template <>                                                                \
  struct TRAIT<double> {                                                     \
    static constexpr double value = Kokkos::Experimental::TRAIT##_v<double>; \
  };                                                                         \
  template <>                                                                \
  struct TRAIT<long double> {                                                \
    static constexpr long double value =                                     \
        Kokkos::Experimental::TRAIT##_v<long double>;                        \
  };                                                                         \
  constexpr float TRAIT<float>::value;                                       \
  constexpr double TRAIT<double>::value;                                     \
  constexpr long double TRAIT<long double>::value
#endif

DEFINE_MATH_CONSTANT_TRAIT(e);
DEFINE_MATH_CONSTANT_TRAIT(log2e);
DEFINE_MATH_CONSTANT_TRAIT(log10e);
DEFINE_MATH_CONSTANT_TRAIT(pi);
DEFINE_MATH_CONSTANT_TRAIT(inv_pi);
DEFINE_MATH_CONSTANT_TRAIT(inv_sqrtpi);
DEFINE_MATH_CONSTANT_TRAIT(ln2);
DEFINE_MATH_CONSTANT_TRAIT(ln10);
DEFINE_MATH_CONSTANT_TRAIT(sqrt2);
DEFINE_MATH_CONSTANT_TRAIT(sqrt3);
DEFINE_MATH_CONSTANT_TRAIT(inv_sqrt3);
DEFINE_MATH_CONSTANT_TRAIT(egamma);
DEFINE_MATH_CONSTANT_TRAIT(phi);

template <class Space, class Trait>
struct TestMathematicalConstants {
  using T = std::decay_t<decltype(Trait::value)>;

  TestMathematicalConstants() { run(); }

  void run() const {
    int errors = 0;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<Space, Trait>(0, 1), *this,
                            errors);
    ASSERT_EQ(errors, 0);
    (void)take_address_of(Trait::value);  // use on host
  }

  KOKKOS_FUNCTION void operator()(Trait, int, int &) const { use_on_device(); }

  KOKKOS_FUNCTION void use_on_device() const {
#if defined(KOKKOS_COMPILER_NVCC) || defined(KOKKOS_ENABLE_OPENMPTARGET)
    take_by_value(Trait::value);
#else
    (void)take_address_of(Trait::value);
#endif
  }
};

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) || \
    defined(KOKKOS_ENABLE_SYCL) || defined(KOKKOS_ENABLE_OPENMPTARGET)
#define TEST_MATH_CONSTANT(TRAIT)                               \
  TEST(TEST_CATEGORY, mathematical_constants_##TRAIT) {         \
    TestMathematicalConstants<TEST_EXECSPACE, TRAIT<float>>();  \
    TestMathematicalConstants<TEST_EXECSPACE, TRAIT<double>>(); \
  }
#else
#define TEST_MATH_CONSTANT(TRAIT)                                    \
  TEST(TEST_CATEGORY, mathematical_constants_##TRAIT) {              \
    TestMathematicalConstants<TEST_EXECSPACE, TRAIT<float>>();       \
    TestMathematicalConstants<TEST_EXECSPACE, TRAIT<double>>();      \
    TestMathematicalConstants<TEST_EXECSPACE, TRAIT<long double>>(); \
  }
#endif

TEST_MATH_CONSTANT(e)
TEST_MATH_CONSTANT(log2e)
TEST_MATH_CONSTANT(log10e)
TEST_MATH_CONSTANT(pi)
TEST_MATH_CONSTANT(inv_pi)
TEST_MATH_CONSTANT(inv_sqrtpi)
TEST_MATH_CONSTANT(ln2)
TEST_MATH_CONSTANT(ln10)
TEST_MATH_CONSTANT(sqrt2)
TEST_MATH_CONSTANT(sqrt3)
TEST_MATH_CONSTANT(inv_sqrt3)
TEST_MATH_CONSTANT(egamma)
TEST_MATH_CONSTANT(phi)
