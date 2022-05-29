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

#ifndef KOKKOS_TEST_VECTOR_HPP
#define KOKKOS_TEST_VECTOR_HPP

#include <gtest/gtest.h>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <Kokkos_Vector.hpp>

namespace Test {

namespace Impl {

template <typename Scalar, class Device>
struct test_vector_insert {
  using scalar_type     = Scalar;
  using execution_space = Device;

  template <typename Vector>
  void run_test(Vector& a) {
    auto n = a.size();

    auto it = a.begin();
    if (n > 0) {
      ASSERT_EQ(a.data(), &a[0]);
    }
    it += 15;
    ASSERT_EQ(*it, scalar_type(1));

    auto it_return = a.insert(it, scalar_type(3));
    ASSERT_EQ(a.size(), n + 1);
    ASSERT_EQ(std::distance(it_return, a.begin() + 15), 0);

    it = a.begin();
    it += 17;
// Looks like some std::vector implementations do not have the restriction
// right on the overload taking three iterators, and thus the following call
// will hit that overload and then fail to compile.
#if defined(KOKKOS_COMPILER_INTEL)
// And at least GCC 4.8.4 doesn't implement vector insert correct for C++11
// Return type is void ...
#if (__GNUC__ < 5)
    a.insert(it, typename Vector::size_type(n + 5), scalar_type(5));
    it_return = a.begin() + 17;
#else
    it_return = a.insert(it, typename Vector::size_type(n + 5), scalar_type(5));
#endif
#else
#if (__GNUC__ < 5)
    a.insert(it, n + 5, scalar_type(5));
    it_return = a.begin() + 17;
#else
    it_return = a.insert(it, n + 5, scalar_type(5));
#endif
#endif

    ASSERT_EQ(a.size(), n + 1 + n + 5);
    ASSERT_EQ(std::distance(it_return, a.begin() + 17), 0u);

    Vector b;

// Looks like some std::vector implementations do not have the restriction
// right on the overload taking three iterators, and thus the following call
// will hit that overload and then fail to compile.
#if defined(KOKKOS_COMPILER_INTEL)
    b.insert(b.begin(), typename Vector::size_type(7), 9);
#else
    b.insert(b.begin(), 7, 9);
#endif
    ASSERT_EQ(b.size(), 7u);
    ASSERT_EQ(b[0], scalar_type(9));

    it = a.begin();
    it += 27 + n;
#if (__GNUC__ < 5)
    a.insert(it, b.begin(), b.end());
    it_return = a.begin() + (27 + n);
#else
    it_return = a.insert(it, b.begin(), b.end());
#endif
    ASSERT_EQ(a.size(), n + 1 + n + 5 + 7);
    ASSERT_EQ(std::distance(it_return, a.begin() + 27 + n), 0u);

    // Testing insert at end via all three function interfaces
    a.insert(a.end(), 11);
#if defined(KOKKOS_COMPILER_INTEL)
    a.insert(a.end(), typename Vector::size_type(2), 12);
#else
    a.insert(a.end(), 2, 12);
#endif
    a.insert(a.end(), b.begin(), b.end());
  }

  template <typename Vector>
  void check_test(Vector& a, int n) {
    for (int i = 0; i < (int)a.size(); i++) {
      if (i == 15)
        ASSERT_EQ(a[i], scalar_type(3));
      else if (i > 16 && i < 16 + 6 + n)
        ASSERT_EQ(a[i], scalar_type(5));
      else if (i > 26 + n && i < 34 + n)
        ASSERT_EQ(a[i], scalar_type(9));
      else if (i == (int)a.size() - 10)
        ASSERT_EQ(a[i], scalar_type(11));
      else if ((i == (int)a.size() - 9) || (i == (int)a.size() - 8))
        ASSERT_EQ(a[i], scalar_type(12));
      else if (i > (int)a.size() - 8)
        ASSERT_EQ(a[i], scalar_type(9));
      else
        ASSERT_EQ(a[i], scalar_type(1));
    }
  }

  test_vector_insert(unsigned int size) {
    {
      std::vector<Scalar> a(size, scalar_type(1));
      run_test(a);
      check_test(a, size);
    }
    {
      Kokkos::vector<Scalar, Device> a(size, scalar_type(1));
      a.sync_device();
      run_test(a);
      a.sync_host();
      check_test(a, size);
    }
    {
      Kokkos::vector<Scalar, Device> a(size, scalar_type(1));
      a.sync_host();
      run_test(a);
      check_test(a, size);
    }
  }
};

template <typename Scalar, class Device>
struct test_vector_allocate {
  using self_type = test_vector_allocate<Scalar, Device>;

  using scalar_type     = Scalar;
  using execution_space = Device;

  bool result = false;

  template <typename Vector>
  Scalar run_me(unsigned int n) {
    {
      Vector v1;
      if (v1.is_allocated() == true) return false;

      v1 = Vector(n, 1);
      Vector v2(v1);
      Vector v3(n, 1);

      if (v1.is_allocated() == false) return false;
      if (v2.is_allocated() == false) return false;
      if (v3.is_allocated() == false) return false;
    }
    return true;
  }

  test_vector_allocate(unsigned int size) {
    result = run_me<Kokkos::vector<Scalar, Device> >(size);
  }
};

template <typename Scalar, class Device>
struct test_vector_combinations {
  using self_type = test_vector_combinations<Scalar, Device>;

  using scalar_type     = Scalar;
  using execution_space = Device;

  Scalar reference;
  Scalar result;

  template <typename Vector>
  Scalar run_me(unsigned int n) {
    Vector a(n, 1);

    a.push_back(2);
    a.resize(n + 4);
    a[n + 1] = 3;
    a[n + 2] = 4;
    a[n + 3] = 5;

    Scalar temp1 = a[2];
    Scalar temp2 = a[n];
    Scalar temp3 = a[n + 1];

    a.assign(n + 2, -1);

    a[2]     = temp1;
    a[n]     = temp2;
    a[n + 1] = temp3;

    Scalar test1 = 0;
    for (unsigned int i = 0; i < a.size(); i++) test1 += a[i];

    a.assign(n + 1, -2);
    Scalar test2 = 0;
    for (unsigned int i = 0; i < a.size(); i++) test2 += a[i];

    a.reserve(n + 10);

    Scalar test3 = 0;
    for (unsigned int i = 0; i < a.size(); i++) test3 += a[i];

    return (test1 * test2 + test3) * test2 + test1 * test3;
  }

  test_vector_combinations(unsigned int size) {
    reference = run_me<std::vector<Scalar> >(size);
    result    = run_me<Kokkos::vector<Scalar, Device> >(size);
  }
};

}  // namespace Impl

template <typename Scalar, typename Device>
void test_vector_combinations(unsigned int size) {
  Impl::test_vector_combinations<Scalar, Device> test(size);
  ASSERT_EQ(test.reference, test.result);
}

template <typename Scalar, typename Device>
void test_vector_allocate(unsigned int size) {
  Impl::test_vector_allocate<Scalar, Device> test(size);
  ASSERT_TRUE(test.result);
}

TEST(TEST_CATEGORY, vector_combination) {
  test_vector_allocate<int, TEST_EXECSPACE>(10);
  test_vector_combinations<int, TEST_EXECSPACE>(10);
  test_vector_combinations<int, TEST_EXECSPACE>(3057);
}

TEST(TEST_CATEGORY, vector_insert) {
  Impl::test_vector_insert<int, TEST_EXECSPACE>(3057);
}

}  // namespace Test

#endif  // KOKKOS_TEST_UNORDERED_MAP_HPP
