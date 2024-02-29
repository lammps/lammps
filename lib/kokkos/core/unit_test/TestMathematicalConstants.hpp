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

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

template <class T>
KOKKOS_FUNCTION T *take_address_of(T &arg) {
  return &arg;
}

template <class T>
KOKKOS_FUNCTION void take_by_value(T) {}

#define DEFINE_MATH_CONSTANT_TRAIT(TRAIT)                     \
  template <class T>                                          \
  struct TRAIT {                                              \
    static constexpr T value = Kokkos::numbers::TRAIT##_v<T>; \
  }

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
#if defined(KOKKOS_COMPILER_NVCC) || defined(KOKKOS_ENABLE_OPENMPTARGET) || \
    defined(KOKKOS_ENABLE_OPENACC) ||                                       \
    defined(KOKKOS_COMPILER_NVHPC)  // FIXME_NVHPC 23.7
    take_by_value(Trait::value);
#else
    (void)take_address_of(Trait::value);
#endif
  }
};

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) ||          \
    defined(KOKKOS_ENABLE_SYCL) || defined(KOKKOS_ENABLE_OPENMPTARGET) || \
    defined(KOKKOS_ENABLE_OPENACC)
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
