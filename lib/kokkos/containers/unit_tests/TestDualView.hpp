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

#ifndef KOKKOS_TEST_DUALVIEW_HPP
#define KOKKOS_TEST_DUALVIEW_HPP

#include <gtest/gtest.h>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <impl/Kokkos_Timer.hpp>
#include <Kokkos_DualView.hpp>

namespace Test {

namespace Impl {
template <typename Scalar, class Device>
struct test_dualview_alloc {
  using scalar_type     = Scalar;
  using execution_space = Device;

  template <typename ViewType>
  bool run_me(unsigned int n, unsigned int m) {
    if (n < 10) n = 10;
    if (m < 3) m = 3;

    {
      ViewType b1;
      if (b1.is_allocated() == true) return false;

      b1 = ViewType("B1", n, m);
      ViewType b2(b1);
      ViewType b3("B3", n, m);

      if (b1.is_allocated() == false) return false;
      if (b2.is_allocated() == false) return false;
      if (b3.is_allocated() == false) return false;
    }
    return true;
  }

  bool result = false;

  test_dualview_alloc(unsigned int size) {
    result = run_me<Kokkos::DualView<Scalar**, Kokkos::LayoutLeft, Device> >(
        size, 3);
  }
};

template <typename Scalar, class Device>
struct test_dualview_combinations {
  using self_type = test_dualview_combinations<Scalar, Device>;

  using scalar_type     = Scalar;
  using execution_space = Device;

  Scalar reference;
  Scalar result;

  template <typename ViewType>
  Scalar run_me(unsigned int n, unsigned int m, bool with_init) {
    if (n < 10) n = 10;
    if (m < 3) m = 3;

    ViewType a;

    if (with_init) {
      a = ViewType("A", n, m);
    } else {
      a = ViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing, "A"), n, m);
    }
    Kokkos::deep_copy(a.d_view, 1);

    a.template modify<typename ViewType::execution_space>();
    a.template sync<typename ViewType::host_mirror_space>();

    a.h_view(5, 1) = 3;
    a.h_view(6, 1) = 4;
    a.h_view(7, 2) = 5;
    a.template modify<typename ViewType::host_mirror_space>();
    ViewType b = Kokkos::subview(a, std::pair<unsigned int, unsigned int>(6, 9),
                                 std::pair<unsigned int, unsigned int>(0, 1));
    a.template sync<typename ViewType::execution_space>();
    b.template modify<typename ViewType::execution_space>();

    Kokkos::deep_copy(b.d_view, 2);

    a.template sync<typename ViewType::host_mirror_space>();
    Scalar count = 0;
    for (unsigned int i = 0; i < a.d_view.extent(0); i++)
      for (unsigned int j = 0; j < a.d_view.extent(1); j++)
        count += a.h_view(i, j);
    return count - a.d_view.extent(0) * a.d_view.extent(1) - 2 - 4 - 3 * 2;
  }

  test_dualview_combinations(unsigned int size, bool with_init) {
    result = run_me<Kokkos::DualView<Scalar**, Kokkos::LayoutLeft, Device> >(
        size, 3, with_init);
  }
};

template <typename Scalar, class ViewType>
struct SumViewEntriesFunctor {
  using value_type = Scalar;

  ViewType fv;

  SumViewEntriesFunctor(const ViewType& fv_) : fv(fv_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, value_type& total) const {
    for (size_t j = 0; j < fv.extent(1); ++j) {
      total += fv(i, j);
    }
  }
};

template <typename Scalar, class Device>
struct test_dual_view_deep_copy {
  using scalar_type     = Scalar;
  using execution_space = Device;

  template <typename ViewType>
  void run_me(int n, const int m, const bool use_templ_sync) {
    ViewType a, b;
    if (n >= 0) {
      a = ViewType("A", n, m);
      b = ViewType("B", n, m);
    } else {
      n = 0;
    }
    const scalar_type sum_total = scalar_type(n * m);

    Kokkos::deep_copy(a.d_view, 1);

    if (use_templ_sync) {
      a.template modify<typename ViewType::execution_space>();
      a.template sync<typename ViewType::host_mirror_space>();
    } else {
      a.modify_device();
      a.sync_host();
    }

    // Check device view is initialized as expected
    scalar_type a_d_sum = 0;
    // Execute on the execution_space associated with t_dev's memory space
    using t_dev_exec_space =
        typename ViewType::t_dev::memory_space::execution_space;
    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<t_dev_exec_space>(0, n),
        SumViewEntriesFunctor<scalar_type, typename ViewType::t_dev>(a.d_view),
        a_d_sum);
    ASSERT_EQ(a_d_sum, sum_total);

    // Check host view is synced as expected
    scalar_type a_h_sum = 0;
    for (size_t i = 0; i < a.h_view.extent(0); ++i)
      for (size_t j = 0; j < a.h_view.extent(1); ++j) {
        a_h_sum += a.h_view(i, j);
      }

    ASSERT_EQ(a_h_sum, sum_total);

    // Test deep_copy
    Kokkos::deep_copy(b, a);
    if (use_templ_sync) {
      b.template sync<typename ViewType::host_mirror_space>();
    } else {
      b.sync_host();
    }

    // Perform same checks on b as done on a
    // Check device view is initialized as expected
    scalar_type b_d_sum = 0;
    // Execute on the execution_space associated with t_dev's memory space
    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<t_dev_exec_space>(0, n),
        SumViewEntriesFunctor<scalar_type, typename ViewType::t_dev>(b.d_view),
        b_d_sum);
    ASSERT_EQ(b_d_sum, sum_total);

    // Check host view is synced as expected
    scalar_type b_h_sum = 0;
    for (size_t i = 0; i < b.h_view.extent(0); ++i)
      for (size_t j = 0; j < b.h_view.extent(1); ++j) {
        b_h_sum += b.h_view(i, j);
      }

    ASSERT_EQ(b_h_sum, sum_total);

  }  // end run_me

  test_dual_view_deep_copy() {
    run_me<Kokkos::DualView<Scalar**, Kokkos::LayoutLeft, Device> >(10, 5,
                                                                    true);
    run_me<Kokkos::DualView<Scalar**, Kokkos::LayoutLeft, Device> >(10, 5,
                                                                    false);
    // Test zero length but allocated (a.d_view.data!=nullptr but
    // a.d_view.span()==0)
    run_me<Kokkos::DualView<Scalar**, Kokkos::LayoutLeft, Device> >(0, 5, true);
    run_me<Kokkos::DualView<Scalar**, Kokkos::LayoutLeft, Device> >(0, 5,
                                                                    false);

    // Test default constructed view
    run_me<Kokkos::DualView<Scalar**, Kokkos::LayoutLeft, Device> >(-1, 5,
                                                                    true);
    run_me<Kokkos::DualView<Scalar**, Kokkos::LayoutLeft, Device> >(-1, 5,
                                                                    false);
  }
};

template <typename Scalar, class Device>
struct test_dualview_resize {
  using scalar_type     = Scalar;
  using execution_space = Device;

  template <typename ViewType>
  void run_me() {
    const unsigned int n      = 10;
    const unsigned int m      = 5;
    const unsigned int factor = 2;

    ViewType a("A", n, m);
    Kokkos::deep_copy(a.d_view, 1);

    /* Covers case "Resize on Device" */
    a.modify_device();
    Kokkos::resize(a, factor * n, factor * m);
    ASSERT_EQ(a.extent(0), n * factor);
    ASSERT_EQ(a.extent(1), m * factor);

    Kokkos::deep_copy(a.d_view, 1);
    a.sync_host();

    // Check device view is initialized as expected
    scalar_type a_d_sum = 0;
    // Execute on the execution_space associated with t_dev's memory space
    using t_dev_exec_space =
        typename ViewType::t_dev::memory_space::execution_space;
    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<t_dev_exec_space>(0, a.d_view.extent(0)),
        SumViewEntriesFunctor<scalar_type, typename ViewType::t_dev>(a.d_view),
        a_d_sum);

    // Check host view is synced as expected
    scalar_type a_h_sum = 0;
    for (size_t i = 0; i < a.h_view.extent(0); ++i)
      for (size_t j = 0; j < a.h_view.extent(1); ++j) {
        a_h_sum += a.h_view(i, j);
      }

    // Check
    ASSERT_EQ(a_h_sum, a_d_sum);
    ASSERT_EQ(a_h_sum, a.extent(0) * a.extent(1));

    /* Covers case "Resize on Host" */
    a.modify_host();

    Kokkos::resize(a, n / factor, m / factor);
    ASSERT_EQ(a.extent(0), n / factor);
    ASSERT_EQ(a.extent(1), m / factor);

    a.sync_device();

    // Check device view is initialized as expected
    a_d_sum = 0;
    // Execute on the execution_space associated with t_dev's memory space
    using t_dev_exec_space =
        typename ViewType::t_dev::memory_space::execution_space;
    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<t_dev_exec_space>(0, a.d_view.extent(0)),
        SumViewEntriesFunctor<scalar_type, typename ViewType::t_dev>(a.d_view),
        a_d_sum);

    // Check host view is synced as expected
    a_h_sum = 0;
    for (size_t i = 0; i < a.h_view.extent(0); ++i)
      for (size_t j = 0; j < a.h_view.extent(1); ++j) {
        a_h_sum += a.h_view(i, j);
      }

    // Check
    ASSERT_EQ(a_h_sum, a.extent(0) * a.extent(1));
    ASSERT_EQ(a_h_sum, a_d_sum);

  }  // end run_me

  test_dualview_resize() {
    run_me<Kokkos::DualView<Scalar**, Kokkos::LayoutLeft, Device> >();
  }
};

template <typename Scalar, class Device>
struct test_dualview_realloc {
  using scalar_type     = Scalar;
  using execution_space = Device;

  template <typename ViewType>
  void run_me() {
    const unsigned int n = 10;
    const unsigned int m = 5;

    ViewType a("A", n, m);
    Kokkos::realloc(a, n, m);

    Kokkos::deep_copy(a.d_view, 1);
    a.modify_device();
    a.sync_host();

    // Check device view is initialized as expected
    scalar_type a_d_sum = 0;
    // Execute on the execution_space associated with t_dev's memory space
    using t_dev_exec_space =
        typename ViewType::t_dev::memory_space::execution_space;
    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<t_dev_exec_space>(0, a.d_view.extent(0)),
        SumViewEntriesFunctor<scalar_type, typename ViewType::t_dev>(a.d_view),
        a_d_sum);

    // Check host view is synced as expected
    scalar_type a_h_sum = 0;
    for (size_t i = 0; i < a.h_view.extent(0); ++i)
      for (size_t j = 0; j < a.h_view.extent(1); ++j) {
        a_h_sum += a.h_view(i, j);
      }

    // Check
    ASSERT_EQ(a_h_sum, a.extent(0) * a.extent(1));
    ASSERT_EQ(a_h_sum, a_d_sum);
  }  // end run_me

  test_dualview_realloc() {
    run_me<Kokkos::DualView<Scalar**, Kokkos::LayoutLeft, Device> >();
  }
};

}  // namespace Impl

template <typename Scalar, typename Device>
void test_dualview_combinations(unsigned int size, bool with_init) {
  Impl::test_dualview_combinations<Scalar, Device> test(size, with_init);
  ASSERT_EQ(test.result, 0);
}

template <typename Scalar, typename Device>
void test_dualview_alloc(unsigned int size) {
  Impl::test_dualview_alloc<Scalar, Device> test(size);
  ASSERT_TRUE(test.result);
}

template <typename Scalar, typename Device>
void test_dualview_deep_copy() {
  Impl::test_dual_view_deep_copy<Scalar, Device>();
}

template <typename Scalar, typename Device>
void test_dualview_realloc() {
  Impl::test_dualview_realloc<Scalar, Device>();
}

template <typename Scalar, typename Device>
void test_dualview_resize() {
  Impl::test_dualview_resize<Scalar, Device>();
}

// FIXME_SYCL requires MDRange policy
#ifndef KOKKOS_ENABLE_SYCL
TEST(TEST_CATEGORY, dualview_combination) {
  test_dualview_combinations<int, TEST_EXECSPACE>(10, true);
}
#endif

TEST(TEST_CATEGORY, dualview_alloc) {
  test_dualview_alloc<int, TEST_EXECSPACE>(10);
}

// FIXME_SYCL requires MDRange policy
#ifndef KOKKOS_ENABLE_SYCL
TEST(TEST_CATEGORY, dualview_combinations_without_init) {
  test_dualview_combinations<int, TEST_EXECSPACE>(10, false);
}

TEST(TEST_CATEGORY, dualview_deep_copy) {
  test_dualview_deep_copy<int, TEST_EXECSPACE>();
  test_dualview_deep_copy<double, TEST_EXECSPACE>();
}

TEST(TEST_CATEGORY, dualview_realloc) {
  test_dualview_realloc<int, TEST_EXECSPACE>();
}

TEST(TEST_CATEGORY, dualview_resize) {
  test_dualview_resize<int, TEST_EXECSPACE>();
}
#endif

}  // namespace Test

#endif  // KOKKOS_TEST_DUALVIEW_HPP
