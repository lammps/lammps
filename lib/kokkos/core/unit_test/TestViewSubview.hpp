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
#ifndef TESTVIEWSUBVIEW_HPP_
#define TESTVIEWSUBVIEW_HPP_
#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>
#include <sstream>
#include <iostream>
#include <type_traits>

// TODO @refactoring move this to somewhere common

//------------------------------------------------------------------------------

template <class...>
struct _kokkos____________________static_test_failure_____;

template <class...>
struct static_predicate_message {};

//------------------------------------------------------------------------------

template <class, template <class...> class, class...>
struct static_assert_predicate_true_impl;

template <template <class...> class predicate, class... message, class... args>
struct static_assert_predicate_true_impl<
    std::enable_if_t<predicate<args...>::type::value>, predicate,
    static_predicate_message<message...>, args...> {
  using type = int;
};

template <template <class...> class predicate, class... message, class... args>
struct static_assert_predicate_true_impl<
    std::enable_if_t<!predicate<args...>::type::value>, predicate,
    static_predicate_message<message...>, args...> {
  using type = typename _kokkos____________________static_test_failure_____<
      message...>::type;
};

template <template <class...> class predicate, class... args>
struct static_assert_predicate_true
    : static_assert_predicate_true_impl<void, predicate,
                                        static_predicate_message<>, args...> {};

template <template <class...> class predicate, class... message, class... args>
struct static_assert_predicate_true<
    predicate, static_predicate_message<message...>, args...>
    : static_assert_predicate_true_impl<
          void, predicate, static_predicate_message<message...>, args...> {};

//------------------------------------------------------------------------------

// error "messages"
struct _kokkos__________types_should_be_the_same_____expected_type__ {};
struct _kokkos__________actual_type_was__ {};
template <class Expected, class Actual>
struct static_expect_same {
  using type = typename static_assert_predicate_true<
      std::is_same,
      static_predicate_message<
          _kokkos__________types_should_be_the_same_____expected_type__,
          Expected, _kokkos__________actual_type_was__, Actual>,
      Expected, Actual>::type;
};

//------------------------------------------------------------------------------

namespace TestViewSubview {

template <class Layout, class Space>
struct getView {
  static Kokkos::View<double**, Layout, Space> get(int n, int m) {
    return Kokkos::View<double**, Layout, Space>("G", n, m);
  }
};

template <class Space>
struct getView<Kokkos::LayoutStride, Space> {
  static Kokkos::View<double**, Kokkos::LayoutStride, Space> get(int n, int m) {
    const int rank       = 2;
    const int order[]    = {0, 1};
    const unsigned dim[] = {unsigned(n), unsigned(m)};
    Kokkos::LayoutStride stride =
        Kokkos::LayoutStride::order_dimensions(rank, order, dim);

    return Kokkos::View<double**, Kokkos::LayoutStride, Space>("G", stride);
  }
};

template <class ViewType, class Space>
struct fill_1D {
  using execution_space = typename Space::execution_space;
  using size_type       = typename ViewType::size_type;

  ViewType a;
  double val;

  fill_1D(ViewType a_, double val_) : a(a_), val(val_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const { a(i) = val; }
};

template <class ViewType, class Space>
struct fill_2D {
  using execution_space = typename Space::execution_space;
  using size_type       = typename ViewType::size_type;

  ViewType a;
  double val;

  fill_2D(ViewType a_, double val_) : a(a_), val(val_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const {
    for (int j = 0; j < static_cast<int>(a.extent(1)); j++) {
      a(i, j) = val;
    }
  }
};

template <class Layout, class Space>
void test_auto_1d() {
  using mv_type         = Kokkos::View<double**, Layout, Space>;
  using execution_space = typename Space::execution_space;
  using size_type       = typename mv_type::size_type;

  const double ZERO = 0.0;
  const double ONE  = 1.0;
  const double TWO  = 2.0;

  const size_type numRows = 10;
  const size_type numCols = 3;

  mv_type X = getView<Layout, Space>::get(numRows, numCols);
  typename mv_type::HostMirror X_h = Kokkos::create_mirror_view(X);

  fill_2D<mv_type, Space> f1(X, ONE);
#if (HIP_VERSION_MAJOR == 5) && (HIP_VERSION_MINOR == 3)
  using Property =
      Kokkos::Experimental::WorkItemProperty::ImplForceGlobalLaunch_t;
#else
  using Property = Kokkos::Experimental::WorkItemProperty::None_t;
#endif
  Kokkos::parallel_for(
      Kokkos::RangePolicy<execution_space, Property>(0, X.extent(0)), f1);
  Kokkos::fence();
  Kokkos::deep_copy(X_h, X);
  for (size_type j = 0; j < numCols; ++j) {
    for (size_type i = 0; i < numRows; ++i) {
      ASSERT_EQ(X_h(i, j), ONE);
    }
  }

  fill_2D<mv_type, Space> f2(X, 0.0);
  Kokkos::parallel_for(
      Kokkos::RangePolicy<execution_space, Property>(0, X.extent(0)), f2);
  Kokkos::fence();
  Kokkos::deep_copy(X_h, X);
  for (size_type j = 0; j < numCols; ++j) {
    for (size_type i = 0; i < numRows; ++i) {
      ASSERT_EQ(X_h(i, j), ZERO);
    }
  }

  fill_2D<mv_type, Space> f3(X, TWO);
  Kokkos::parallel_for(
      Kokkos::RangePolicy<execution_space, Property>(0, X.extent(0)), f3);
  Kokkos::fence();
  Kokkos::deep_copy(X_h, X);
  for (size_type j = 0; j < numCols; ++j) {
    for (size_type i = 0; i < numRows; ++i) {
      ASSERT_EQ(X_h(i, j), TWO);
    }
  }

  for (size_type j = 0; j < numCols; ++j) {
    auto X_j = Kokkos::subview(X, Kokkos::ALL, j);

    fill_1D<decltype(X_j), Space> f4(X_j, ZERO);
    Kokkos::parallel_for(
        Kokkos::RangePolicy<execution_space, Property>(0, X_j.extent(0)), f4);
    Kokkos::fence();
    Kokkos::deep_copy(X_h, X);
    for (size_type i = 0; i < numRows; ++i) {
      ASSERT_EQ(X_h(i, j), ZERO);
    }

    for (size_type jj = 0; jj < numCols; ++jj) {
      auto X_jj = Kokkos::subview(X, Kokkos::ALL, jj);
      fill_1D<decltype(X_jj), Space> f5(X_jj, ONE);
      Kokkos::parallel_for(
          Kokkos::RangePolicy<execution_space, Property>(0, X_jj.extent(0)),
          f5);
      Kokkos::fence();
      Kokkos::deep_copy(X_h, X);
      for (size_type i = 0; i < numRows; ++i) {
        ASSERT_EQ(X_h(i, jj), ONE);
      }
    }
  }
}

template <class LD, class LS, class Space>
void test_1d_strided_assignment_impl(bool a, bool b, bool c, bool d, int n,
                                     int m) {
  Kokkos::View<double**, LS, Space> l2d("l2d", n, m);

  int col = n > 2 ? 2 : 0;
  int row = m > 2 ? 2 : 0;

  if (Kokkos::SpaceAccessibility<Kokkos::HostSpace,
                                 typename Space::memory_space>::accessible) {
    if (a) {
      Kokkos::View<double*, LD, Space> l1da =
          Kokkos::subview(l2d, Kokkos::ALL, row);
      ASSERT_EQ(&l1da(0), &l2d(0, row));
      if (n > 1) {
        ASSERT_EQ(&l1da(1), &l2d(1, row));
      }
    }

    if (b && n > 13) {
      Kokkos::View<double*, LD, Space> l1db =
          Kokkos::subview(l2d, std::pair<unsigned, unsigned>(2, 13), row);
      ASSERT_EQ(&l1db(0), &l2d(2, row));
      ASSERT_EQ(&l1db(1), &l2d(3, row));
    }

    if (c) {
      Kokkos::View<double*, LD, Space> l1dc =
          Kokkos::subview(l2d, col, Kokkos::ALL);
      ASSERT_EQ(&l1dc(0), &l2d(col, 0));
      if (m > 1) {
        ASSERT_EQ(&l1dc(1), &l2d(col, 1));
      }
    }

    if (d && m > 13) {
      Kokkos::View<double*, LD, Space> l1dd =
          Kokkos::subview(l2d, col, std::pair<unsigned, unsigned>(2, 13));
      ASSERT_EQ(&l1dd(0), &l2d(col, 2));
      ASSERT_EQ(&l1dd(1), &l2d(col, 3));
    }
  }
}

template <class Space>
void test_1d_strided_assignment() {
  test_1d_strided_assignment_impl<Kokkos::LayoutStride, Kokkos::LayoutLeft,
                                  Space>(true, true, true, true, 17, 3);
  test_1d_strided_assignment_impl<Kokkos::LayoutStride, Kokkos::LayoutRight,
                                  Space>(true, true, true, true, 17, 3);

  test_1d_strided_assignment_impl<Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                  Space>(true, true, false, false, 17, 3);
  test_1d_strided_assignment_impl<Kokkos::LayoutRight, Kokkos::LayoutLeft,
                                  Space>(true, true, false, false, 17, 3);
  test_1d_strided_assignment_impl<Kokkos::LayoutLeft, Kokkos::LayoutRight,
                                  Space>(false, false, true, true, 17, 3);
  test_1d_strided_assignment_impl<Kokkos::LayoutRight, Kokkos::LayoutRight,
                                  Space>(false, false, true, true, 17, 3);

  test_1d_strided_assignment_impl<Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                  Space>(true, true, false, false, 17, 1);
  test_1d_strided_assignment_impl<Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                  Space>(true, true, true, true, 1, 17);
  test_1d_strided_assignment_impl<Kokkos::LayoutRight, Kokkos::LayoutLeft,
                                  Space>(true, true, true, true, 1, 17);
  test_1d_strided_assignment_impl<Kokkos::LayoutRight, Kokkos::LayoutLeft,
                                  Space>(true, true, false, false, 17, 1);

  test_1d_strided_assignment_impl<Kokkos::LayoutLeft, Kokkos::LayoutRight,
                                  Space>(true, true, true, true, 17, 1);
  test_1d_strided_assignment_impl<Kokkos::LayoutLeft, Kokkos::LayoutRight,
                                  Space>(false, false, true, true, 1, 17);
  test_1d_strided_assignment_impl<Kokkos::LayoutRight, Kokkos::LayoutRight,
                                  Space>(false, false, true, true, 1, 17);
  test_1d_strided_assignment_impl<Kokkos::LayoutRight, Kokkos::LayoutRight,
                                  Space>(true, true, true, true, 17, 1);
}

template <class NewView, class OrigView, class... Args>
void make_subview(bool use_constructor, NewView& v, OrigView org,
                  Args... args) {
  if (use_constructor) {
    v = NewView(org, args...);
  } else {
    v = Kokkos::subview(org, args...);
  }
}

template <class Space>
void test_left_0(bool constr) {
  using view_static_8_type =
      Kokkos::View<int[2][3][4][5][2][3][4][5], Kokkos::LayoutLeft, Space>;

  if (Kokkos::SpaceAccessibility<Kokkos::HostSpace,
                                 typename Space::memory_space>::accessible) {
    view_static_8_type x_static_8("x_static_left_8");

    ASSERT_TRUE(x_static_8.span_is_contiguous());

    Kokkos::View<int, Kokkos::LayoutLeft, Space> x0;
    make_subview(constr, x0, x_static_8, 0, 0, 0, 0, 0, 0, 0, 0);

    ASSERT_TRUE(x0.span_is_contiguous());
    ASSERT_EQ(x0.span(), 1u);
    ASSERT_EQ(&x0(), &x_static_8(0, 0, 0, 0, 0, 0, 0, 0));

    Kokkos::View<int*, Kokkos::LayoutLeft, Space> x1;
    make_subview(constr, x1, x_static_8, Kokkos::pair<int, int>(0, 2), 1, 2, 3,
                 0, 1, 2, 3);

    ASSERT_TRUE(x1.span_is_contiguous());
    ASSERT_EQ(x1.span(), 2u);
    ASSERT_EQ(&x1(0), &x_static_8(0, 1, 2, 3, 0, 1, 2, 3));
    ASSERT_EQ(&x1(1), &x_static_8(1, 1, 2, 3, 0, 1, 2, 3));

    Kokkos::View<int*, Kokkos::LayoutLeft, Space> x_deg1;
    make_subview(constr, x_deg1, x_static_8, Kokkos::pair<int, int>(0, 0), 1, 2,
                 3, 0, 1, 2, 3);

    ASSERT_TRUE(x_deg1.span_is_contiguous());
    ASSERT_EQ(x_deg1.span(), 0u);
    ASSERT_EQ(x_deg1.data(), &x_static_8(0, 1, 2, 3, 0, 1, 2, 3));

    Kokkos::View<int*, Kokkos::LayoutLeft, Space> x_deg2;
    make_subview(constr, x_deg2, x_static_8, Kokkos::pair<int, int>(2, 2), 2, 3,
                 4, 1, 2, 3, 4);

    ASSERT_TRUE(x_deg2.span_is_contiguous());
    ASSERT_EQ(x_deg2.span(), 0u);
    ASSERT_EQ(x_deg2.data(), x_static_8.data() + x_static_8.span());

    Kokkos::View<int**, Kokkos::LayoutLeft, Space> x2;
    make_subview(constr, x2, x_static_8, Kokkos::pair<int, int>(0, 2), 1, 2, 3,
                 Kokkos::pair<int, int>(0, 2), 1, 2, 3);

    ASSERT_TRUE(!x2.span_is_contiguous());
    ASSERT_EQ(&x2(0, 0), &x_static_8(0, 1, 2, 3, 0, 1, 2, 3));
    ASSERT_EQ(&x2(1, 0), &x_static_8(1, 1, 2, 3, 0, 1, 2, 3));
    ASSERT_EQ(&x2(0, 1), &x_static_8(0, 1, 2, 3, 1, 1, 2, 3));
    ASSERT_EQ(&x2(1, 1), &x_static_8(1, 1, 2, 3, 1, 1, 2, 3));

    // Kokkos::View< int**, Kokkos::LayoutLeft, Space > error_2 =
    Kokkos::View<int**, Kokkos::LayoutStride, Space> sx2;
    make_subview(constr, sx2, x_static_8, 1, Kokkos::pair<int, int>(0, 2), 2, 3,
                 Kokkos::pair<int, int>(0, 2), 1, 2, 3);

    ASSERT_TRUE(!sx2.span_is_contiguous());
    ASSERT_EQ(&sx2(0, 0), &x_static_8(1, 0, 2, 3, 0, 1, 2, 3));
    ASSERT_EQ(&sx2(1, 0), &x_static_8(1, 1, 2, 3, 0, 1, 2, 3));
    ASSERT_EQ(&sx2(0, 1), &x_static_8(1, 0, 2, 3, 1, 1, 2, 3));
    ASSERT_EQ(&sx2(1, 1), &x_static_8(1, 1, 2, 3, 1, 1, 2, 3));

    Kokkos::View<int****, Kokkos::LayoutStride, Space> sx4;
    make_subview(constr, sx4, x_static_8, 0,
                 Kokkos::pair<int, int>(0, 2) /* of [3] */
                 ,
                 1, Kokkos::pair<int, int>(1, 3) /* of [5] */
                 ,
                 1, Kokkos::pair<int, int>(0, 2) /* of [3] */
                 ,
                 2, Kokkos::pair<int, int>(2, 4) /* of [5] */
    );

    ASSERT_TRUE(!sx4.span_is_contiguous());

    for (int i0 = 0; i0 < (int)sx4.extent(0); ++i0)
      for (int i1 = 0; i1 < (int)sx4.extent(1); ++i1)
        for (int i2 = 0; i2 < (int)sx4.extent(2); ++i2)
          for (int i3 = 0; i3 < (int)sx4.extent(3); ++i3) {
            ASSERT_EQ(&sx4(i0, i1, i2, i3),
                      &x_static_8(0, 0 + i0, 1, 1 + i1, 1, 0 + i2, 2, 2 + i3));
          }
  }
}

template <class Space>
void test_left_0() {
  test_left_0<Space>(true);
  test_left_0<Space>(false);
}

template <class Space>
void test_left_1(bool use_constr) {
  using view_type =
      Kokkos::View<int*** * [2][3][4][5], Kokkos::LayoutLeft, Space>;

  if (Kokkos::SpaceAccessibility<Kokkos::HostSpace,
                                 typename Space::memory_space>::accessible) {
    view_type x8("x_left_8", 2, 3, 4, 5);

    ASSERT_TRUE(x8.span_is_contiguous());

    Kokkos::View<int, Kokkos::LayoutLeft, Space> x0;
    make_subview(use_constr, x0, x8, 0, 0, 0, 0, 0, 0, 0, 0);

    ASSERT_TRUE(x0.span_is_contiguous());
    ASSERT_EQ(&x0(), &x8(0, 0, 0, 0, 0, 0, 0, 0));

    Kokkos::View<int*, Kokkos::LayoutLeft, Space> x1;
    make_subview(use_constr, x1, x8, Kokkos::pair<int, int>(0, 2), 1, 2, 3, 0,
                 1, 2, 3);

    ASSERT_TRUE(x1.span_is_contiguous());
    ASSERT_EQ(&x1(0), &x8(0, 1, 2, 3, 0, 1, 2, 3));
    ASSERT_EQ(&x1(1), &x8(1, 1, 2, 3, 0, 1, 2, 3));

    Kokkos::View<int*, Kokkos::LayoutLeft, Space> x1_deg1;
    make_subview(use_constr, x1_deg1, x8, Kokkos::pair<int, int>(0, 0), 1, 2, 3,
                 0, 1, 2, 3);

    ASSERT_TRUE(x1_deg1.span_is_contiguous());
    ASSERT_EQ(0u, x1_deg1.span());
    ASSERT_EQ(x1_deg1.data(), &x8(0, 1, 2, 3, 0, 1, 2, 3));

    Kokkos::View<int*, Kokkos::LayoutLeft, Space> x1_deg2;
    make_subview(use_constr, x1_deg2, x8, Kokkos::pair<int, int>(2, 2), 2, 3, 4,
                 1, 2, 3, 4);

    ASSERT_EQ(0u, x1_deg2.span());
    ASSERT_TRUE(x1_deg2.span_is_contiguous());
    ASSERT_EQ(x1_deg2.data(), x8.data() + x8.span());

    Kokkos::View<int**, Kokkos::LayoutLeft, Space> x2;
    make_subview(use_constr, x2, x8, Kokkos::pair<int, int>(0, 2), 1, 2, 3,
                 Kokkos::pair<int, int>(0, 2), 1, 2, 3);

    ASSERT_TRUE(!x2.span_is_contiguous());
    ASSERT_EQ(&x2(0, 0), &x8(0, 1, 2, 3, 0, 1, 2, 3));
    ASSERT_EQ(&x2(1, 0), &x8(1, 1, 2, 3, 0, 1, 2, 3));
    ASSERT_EQ(&x2(0, 1), &x8(0, 1, 2, 3, 1, 1, 2, 3));
    ASSERT_EQ(&x2(1, 1), &x8(1, 1, 2, 3, 1, 1, 2, 3));

    Kokkos::View<int**, Kokkos::LayoutLeft, Space> x2_deg2;
    make_subview(use_constr, x2_deg2, x8, Kokkos::pair<int, int>(2, 2), 2, 3, 4,
                 1, 2, Kokkos::pair<int, int>(2, 3), 4);
    ASSERT_EQ(0u, x2_deg2.span());

    // Kokkos::View< int**, Kokkos::LayoutLeft, Space > error_2 =
    Kokkos::View<int**, Kokkos::LayoutStride, Space> sx2;
    make_subview(use_constr, sx2, x8, 1, Kokkos::pair<int, int>(0, 2), 2, 3,
                 Kokkos::pair<int, int>(0, 2), 1, 2, 3);

    ASSERT_TRUE(!sx2.span_is_contiguous());
    ASSERT_EQ(&sx2(0, 0), &x8(1, 0, 2, 3, 0, 1, 2, 3));
    ASSERT_EQ(&sx2(1, 0), &x8(1, 1, 2, 3, 0, 1, 2, 3));
    ASSERT_EQ(&sx2(0, 1), &x8(1, 0, 2, 3, 1, 1, 2, 3));
    ASSERT_EQ(&sx2(1, 1), &x8(1, 1, 2, 3, 1, 1, 2, 3));

    Kokkos::View<int**, Kokkos::LayoutStride, Space> sx2_deg;
    make_subview(use_constr, sx2, x8, 1, Kokkos::pair<int, int>(0, 0), 2, 3,
                 Kokkos::pair<int, int>(0, 2), 1, 2, 3);
    ASSERT_EQ(0u, sx2_deg.span());

    Kokkos::View<int****, Kokkos::LayoutStride, Space> sx4;
    make_subview(use_constr, sx4, x8, 0,
                 Kokkos::pair<int, int>(0, 2) /* of [3] */
                 ,
                 1, Kokkos::pair<int, int>(1, 3) /* of [5] */
                 ,
                 1, Kokkos::pair<int, int>(0, 2) /* of [3] */
                 ,
                 2, Kokkos::pair<int, int>(2, 4) /* of [5] */
    );

    ASSERT_TRUE(!sx4.span_is_contiguous());

    for (int i0 = 0; i0 < (int)sx4.extent(0); ++i0)
      for (int i1 = 0; i1 < (int)sx4.extent(1); ++i1)
        for (int i2 = 0; i2 < (int)sx4.extent(2); ++i2)
          for (int i3 = 0; i3 < (int)sx4.extent(3); ++i3) {
            ASSERT_EQ(&sx4(i0, i1, i2, i3),
                      &x8(0, 0 + i0, 1, 1 + i1, 1, 0 + i2, 2, 2 + i3));
          }
  }
}

template <class Space>
void test_left_1() {
  test_left_1<Space>(true);
  test_left_1<Space>(false);
}

template <class Space>
void test_left_2() {
  using view_type = Kokkos::View<int****, Kokkos::LayoutLeft, Space>;

  if (Kokkos::SpaceAccessibility<Kokkos::HostSpace,
                                 typename Space::memory_space>::accessible) {
    view_type x4("x4", 2, 3, 4, 5);

    ASSERT_TRUE(x4.span_is_contiguous());

    Kokkos::View<int, Kokkos::LayoutLeft, Space> x0 =
        Kokkos::subview(x4, 0, 0, 0, 0);

    ASSERT_TRUE(x0.span_is_contiguous());
    ASSERT_EQ(&x0(), &x4(0, 0, 0, 0));

    Kokkos::View<int*, Kokkos::LayoutLeft, Space> x1 =
        Kokkos::subview(x4, Kokkos::pair<int, int>(0, 2), 1, 2, 3);

    ASSERT_TRUE(x1.span_is_contiguous());
    ASSERT_EQ(&x1(0), &x4(0, 1, 2, 3));
    ASSERT_EQ(&x1(1), &x4(1, 1, 2, 3));

    Kokkos::View<int**, Kokkos::LayoutLeft, Space> x2 = Kokkos::subview(
        x4, Kokkos::pair<int, int>(0, 2), 1, Kokkos::pair<int, int>(1, 3), 2);

    ASSERT_TRUE(!x2.span_is_contiguous());
    ASSERT_EQ(&x2(0, 0), &x4(0, 1, 1, 2));
    ASSERT_EQ(&x2(1, 0), &x4(1, 1, 1, 2));
    ASSERT_EQ(&x2(0, 1), &x4(0, 1, 2, 2));
    ASSERT_EQ(&x2(1, 1), &x4(1, 1, 2, 2));

    // Kokkos::View< int**, Kokkos::LayoutLeft, Space > error_2 =
    Kokkos::View<int**, Kokkos::LayoutStride, Space> sx2 = Kokkos::subview(
        x4, 1, Kokkos::pair<int, int>(0, 2), 2, Kokkos::pair<int, int>(1, 4));

    ASSERT_TRUE(!sx2.span_is_contiguous());
    ASSERT_EQ(&sx2(0, 0), &x4(1, 0, 2, 1));
    ASSERT_EQ(&sx2(1, 0), &x4(1, 1, 2, 1));
    ASSERT_EQ(&sx2(0, 1), &x4(1, 0, 2, 2));
    ASSERT_EQ(&sx2(1, 1), &x4(1, 1, 2, 2));
    ASSERT_EQ(&sx2(0, 2), &x4(1, 0, 2, 3));
    ASSERT_EQ(&sx2(1, 2), &x4(1, 1, 2, 3));

    Kokkos::View<int****, Kokkos::LayoutStride, Space> sx4 =
        Kokkos::subview(x4, Kokkos::pair<int, int>(1, 2) /* of [2] */
                        ,
                        Kokkos::pair<int, int>(1, 3) /* of [3] */
                        ,
                        Kokkos::pair<int, int>(0, 4) /* of [4] */
                        ,
                        Kokkos::pair<int, int>(2, 4) /* of [5] */
        );

    ASSERT_TRUE(!sx4.span_is_contiguous());

    for (int i0 = 0; i0 < (int)sx4.extent(0); ++i0)
      for (int i1 = 0; i1 < (int)sx4.extent(1); ++i1)
        for (int i2 = 0; i2 < (int)sx4.extent(2); ++i2)
          for (int i3 = 0; i3 < (int)sx4.extent(3); ++i3) {
            ASSERT_EQ(&sx4(i0, i1, i2, i3),
                      &x4(1 + i0, 1 + i1, 0 + i2, 2 + i3));
          }
  }
}

template <class Space>
void test_left_3() {
  using view_type = Kokkos::View<int**, Kokkos::LayoutLeft, Space>;

  if (Kokkos::SpaceAccessibility<Kokkos::HostSpace,
                                 typename Space::memory_space>::accessible) {
    view_type xm("x4", 10, 5);

    ASSERT_TRUE(xm.span_is_contiguous());

    Kokkos::View<int, Kokkos::LayoutLeft, Space> x0 = Kokkos::subview(xm, 5, 3);

    ASSERT_TRUE(x0.span_is_contiguous());
    ASSERT_EQ(&x0(), &xm(5, 3));

    Kokkos::View<int*, Kokkos::LayoutLeft, Space> x1 =
        Kokkos::subview(xm, Kokkos::ALL, 3);

    ASSERT_TRUE(x1.span_is_contiguous());
    for (int i = 0; i < int(xm.extent(0)); ++i) {
      ASSERT_EQ(&x1(i), &xm(i, 3));
    }

    Kokkos::View<int**, Kokkos::LayoutLeft, Space> x2 =
        Kokkos::subview(xm, Kokkos::pair<int, int>(1, 9), Kokkos::ALL);

    ASSERT_TRUE(!x2.span_is_contiguous());
    for (int j = 0; j < int(x2.extent(1)); ++j)
      for (int i = 0; i < int(x2.extent(0)); ++i) {
        ASSERT_EQ(&x2(i, j), &xm(1 + i, j));
      }

    Kokkos::View<int**, Kokkos::LayoutLeft, Space> x2c =
        Kokkos::subview(xm, Kokkos::ALL, std::pair<int, int>(2, 4));

    ASSERT_TRUE(x2c.span_is_contiguous());
    for (int j = 0; j < int(x2c.extent(1)); ++j)
      for (int i = 0; i < int(x2c.extent(0)); ++i) {
        ASSERT_EQ(&x2c(i, j), &xm(i, 2 + j));
      }

    Kokkos::View<int**, Kokkos::LayoutLeft, Space> x2_n1 =
        Kokkos::subview(xm, std::pair<int, int>(1, 1), Kokkos::ALL);

    ASSERT_EQ(x2_n1.extent(0), 0u);
    ASSERT_EQ(x2_n1.extent(1), xm.extent(1));

    Kokkos::View<int**, Kokkos::LayoutLeft, Space> x2_n2 =
        Kokkos::subview(xm, Kokkos::ALL, std::pair<int, int>(1, 1));

    ASSERT_EQ(x2_n2.extent(0), xm.extent(0));
    ASSERT_EQ(x2_n2.extent(1), 0u);
  }
}

//----------------------------------------------------------------------------

template <class Space>
void test_right_0(bool use_constr) {
  using view_static_8_type =
      Kokkos::View<int[2][3][4][5][2][3][4][5], Kokkos::LayoutRight, Space>;

  if (Kokkos::SpaceAccessibility<Kokkos::HostSpace,
                                 typename Space::memory_space>::accessible) {
    view_static_8_type x_static_8("x_static_right_8");

    Kokkos::View<int, Kokkos::LayoutRight, Space> x0;
    make_subview(use_constr, x0, x_static_8, 0, 0, 0, 0, 0, 0, 0, 0);

    ASSERT_EQ(&x0(), &x_static_8(0, 0, 0, 0, 0, 0, 0, 0));

    Kokkos::View<int*, Kokkos::LayoutRight, Space> x1;
    make_subview(use_constr, x1, x_static_8, 0, 1, 2, 3, 0, 1, 2,
                 Kokkos::pair<int, int>(1, 3));

    ASSERT_EQ(x1.extent(0), 2u);
    ASSERT_EQ(&x1(0), &x_static_8(0, 1, 2, 3, 0, 1, 2, 1));
    ASSERT_EQ(&x1(1), &x_static_8(0, 1, 2, 3, 0, 1, 2, 2));

    Kokkos::View<int**, Kokkos::LayoutRight, Space> x2;
    make_subview(use_constr, x2, x_static_8, 0, 1, 2,
                 Kokkos::pair<int, int>(1, 3), 0, 1, 2,
                 Kokkos::pair<int, int>(1, 3));

    ASSERT_EQ(x2.extent(0), 2u);
    ASSERT_EQ(x2.extent(1), 2u);
    ASSERT_EQ(&x2(0, 0), &x_static_8(0, 1, 2, 1, 0, 1, 2, 1));
    ASSERT_EQ(&x2(1, 0), &x_static_8(0, 1, 2, 2, 0, 1, 2, 1));
    ASSERT_EQ(&x2(0, 1), &x_static_8(0, 1, 2, 1, 0, 1, 2, 2));
    ASSERT_EQ(&x2(1, 1), &x_static_8(0, 1, 2, 2, 0, 1, 2, 2));

    // Kokkos::View< int**, Kokkos::LayoutRight, Space > error_2 =
    Kokkos::View<int**, Kokkos::LayoutStride, Space> sx2;
    make_subview(use_constr, sx2, x_static_8, 1, Kokkos::pair<int, int>(0, 2),
                 2, 3, Kokkos::pair<int, int>(0, 2), 1, 2, 3);

    ASSERT_EQ(sx2.extent(0), 2u);
    ASSERT_EQ(sx2.extent(1), 2u);
    ASSERT_EQ(&sx2(0, 0), &x_static_8(1, 0, 2, 3, 0, 1, 2, 3));
    ASSERT_EQ(&sx2(1, 0), &x_static_8(1, 1, 2, 3, 0, 1, 2, 3));
    ASSERT_EQ(&sx2(0, 1), &x_static_8(1, 0, 2, 3, 1, 1, 2, 3));
    ASSERT_EQ(&sx2(1, 1), &x_static_8(1, 1, 2, 3, 1, 1, 2, 3));

    Kokkos::View<int****, Kokkos::LayoutStride, Space> sx4;
    make_subview(use_constr, sx4, x_static_8, 0,
                 Kokkos::pair<int, int>(0, 2) /* of [3] */
                 ,
                 1, Kokkos::pair<int, int>(1, 3) /* of [5] */
                 ,
                 1, Kokkos::pair<int, int>(0, 2) /* of [3] */
                 ,
                 2, Kokkos::pair<int, int>(2, 4) /* of [5] */
    );

    ASSERT_EQ(sx4.extent(0), 2u);
    ASSERT_EQ(sx4.extent(1), 2u);
    ASSERT_EQ(sx4.extent(2), 2u);
    ASSERT_EQ(sx4.extent(3), 2u);
    for (int i0 = 0; i0 < (int)sx4.extent(0); ++i0)
      for (int i1 = 0; i1 < (int)sx4.extent(1); ++i1)
        for (int i2 = 0; i2 < (int)sx4.extent(2); ++i2)
          for (int i3 = 0; i3 < (int)sx4.extent(3); ++i3) {
            ASSERT_EQ(&sx4(i0, i1, i2, i3),
                      &x_static_8(0, 0 + i0, 1, 1 + i1, 1, 0 + i2, 2, 2 + i3));
          }
  }
}

template <class Space>
void test_right_0() {
  test_right_0<Space>(true);
  test_right_0<Space>(false);
}

template <class Space>
void test_right_1(bool use_constr) {
  using view_type =
      Kokkos::View<int*** * [2][3][4][5], Kokkos::LayoutRight, Space>;

  if (Kokkos::SpaceAccessibility<Kokkos::HostSpace,
                                 typename Space::memory_space>::accessible) {
    view_type x8("x_right_8", 2, 3, 4, 5);

    Kokkos::View<int, Kokkos::LayoutRight, Space> x0;
    make_subview(use_constr, x0, x8, 0, 0, 0, 0, 0, 0, 0, 0);

    ASSERT_EQ(&x0(), &x8(0, 0, 0, 0, 0, 0, 0, 0));

    Kokkos::View<int*, Kokkos::LayoutRight, Space> x1;
    make_subview(use_constr, x1, x8, 0, 1, 2, 3, 0, 1, 2,
                 Kokkos::pair<int, int>(1, 3));

    ASSERT_EQ(&x1(0), &x8(0, 1, 2, 3, 0, 1, 2, 1));
    ASSERT_EQ(&x1(1), &x8(0, 1, 2, 3, 0, 1, 2, 2));

    Kokkos::View<int*, Kokkos::LayoutRight, Space> x1_deg1;
    make_subview(use_constr, x1_deg1, x8, 0, 1, 2, 3, 0, 1, 2,
                 Kokkos::pair<int, int>(3, 3));
    ASSERT_EQ(0u, x1_deg1.span());

    Kokkos::View<int**, Kokkos::LayoutRight, Space> x2;
    make_subview(use_constr, x2, x8, 0, 1, 2, Kokkos::pair<int, int>(1, 3), 0,
                 1, 2, Kokkos::pair<int, int>(1, 3));

    ASSERT_EQ(&x2(0, 0), &x8(0, 1, 2, 1, 0, 1, 2, 1));
    ASSERT_EQ(&x2(1, 0), &x8(0, 1, 2, 2, 0, 1, 2, 1));
    ASSERT_EQ(&x2(0, 1), &x8(0, 1, 2, 1, 0, 1, 2, 2));
    ASSERT_EQ(&x2(1, 1), &x8(0, 1, 2, 2, 0, 1, 2, 2));

    Kokkos::View<int**, Kokkos::LayoutRight, Space> x2_deg2;
    make_subview(use_constr, x2_deg2, x8, 0, 1, 2, Kokkos::pair<int, int>(1, 3),
                 0, 1, 2, Kokkos::pair<int, int>(3, 3));
    ASSERT_EQ(0u, x2_deg2.span());

    // Kokkos::View< int**, Kokkos::LayoutRight, Space > error_2 =
    Kokkos::View<int**, Kokkos::LayoutStride, Space> sx2;
    make_subview(use_constr, sx2, x8, 1, Kokkos::pair<int, int>(0, 2), 2, 3,
                 Kokkos::pair<int, int>(0, 2), 1, 2, 3);

    ASSERT_EQ(&sx2(0, 0), &x8(1, 0, 2, 3, 0, 1, 2, 3));
    ASSERT_EQ(&sx2(1, 0), &x8(1, 1, 2, 3, 0, 1, 2, 3));
    ASSERT_EQ(&sx2(0, 1), &x8(1, 0, 2, 3, 1, 1, 2, 3));
    ASSERT_EQ(&sx2(1, 1), &x8(1, 1, 2, 3, 1, 1, 2, 3));

    Kokkos::View<int**, Kokkos::LayoutStride, Space> sx2_deg;
    make_subview(use_constr, sx2_deg, x8, 1, Kokkos::pair<int, int>(0, 2), 2, 3,
                 1, 1, 2, Kokkos::pair<int, int>(3, 3));
    ASSERT_EQ(0u, sx2_deg.span());

    Kokkos::View<int****, Kokkos::LayoutStride, Space> sx4;
    make_subview(use_constr, sx4, x8, 0,
                 Kokkos::pair<int, int>(0, 2) /* of [3] */
                 ,
                 1, Kokkos::pair<int, int>(1, 3) /* of [5] */
                 ,
                 1, Kokkos::pair<int, int>(0, 2) /* of [3] */
                 ,
                 2, Kokkos::pair<int, int>(2, 4) /* of [5] */
    );

    for (int i0 = 0; i0 < (int)sx4.extent(0); ++i0)
      for (int i1 = 0; i1 < (int)sx4.extent(1); ++i1)
        for (int i2 = 0; i2 < (int)sx4.extent(2); ++i2)
          for (int i3 = 0; i3 < (int)sx4.extent(3); ++i3) {
            ASSERT_EQ(&sx4(i0, i1, i2, i3),
                      &x8(0, 0 + i0, 1, 1 + i1, 1, 0 + i2, 2, 2 + i3));
          }
  }
}

template <class Space>
void test_right_1() {
  test_right_1<Space>(true);
  test_right_1<Space>(false);
}

template <class Space>
void test_right_3() {
  using view_type = Kokkos::View<int**, Kokkos::LayoutRight, Space>;

  if (Kokkos::SpaceAccessibility<Kokkos::HostSpace,
                                 typename Space::memory_space>::accessible) {
    view_type xm("x4", 10, 5);

    ASSERT_TRUE(xm.span_is_contiguous());

    Kokkos::View<int, Kokkos::LayoutRight, Space> x0 =
        Kokkos::subview(xm, 5, 3);

    ASSERT_TRUE(x0.span_is_contiguous());
    ASSERT_EQ(&x0(), &xm(5, 3));

    Kokkos::View<int*, Kokkos::LayoutRight, Space> x1 =
        Kokkos::subview(xm, 3, Kokkos::ALL);

    ASSERT_TRUE(x1.span_is_contiguous());
    for (int i = 0; i < int(xm.extent(1)); ++i) {
      ASSERT_EQ(&x1(i), &xm(3, i));
    }

    Kokkos::View<int**, Kokkos::LayoutRight, Space> x2c =
        Kokkos::subview(xm, Kokkos::pair<int, int>(1, 9), Kokkos::ALL);

    ASSERT_TRUE(x2c.span_is_contiguous());
    for (int j = 0; j < int(x2c.extent(1)); ++j)
      for (int i = 0; i < int(x2c.extent(0)); ++i) {
        ASSERT_EQ(&x2c(i, j), &xm(1 + i, j));
      }

    Kokkos::View<int**, Kokkos::LayoutRight, Space> x2 =
        Kokkos::subview(xm, Kokkos::ALL, std::pair<int, int>(2, 4));

    ASSERT_TRUE(!x2.span_is_contiguous());
    for (int j = 0; j < int(x2.extent(1)); ++j)
      for (int i = 0; i < int(x2.extent(0)); ++i) {
        ASSERT_EQ(&x2(i, j), &xm(i, 2 + j));
      }

    Kokkos::View<int**, Kokkos::LayoutRight, Space> x2_n1 =
        Kokkos::subview(xm, std::pair<int, int>(1, 1), Kokkos::ALL);

    ASSERT_EQ(x2_n1.extent(0), 0u);
    ASSERT_EQ(x2_n1.extent(1), xm.extent(1));

    Kokkos::View<int**, Kokkos::LayoutRight, Space> x2_n2 =
        Kokkos::subview(xm, Kokkos::ALL, std::pair<int, int>(1, 1));

    ASSERT_EQ(x2_n2.extent(0), xm.extent(0));
    ASSERT_EQ(x2_n2.extent(1), 0u);
  }
}

namespace Impl {

constexpr int N0 = 113;
constexpr int N1 = 11;
constexpr int N2 = 17;
constexpr int N3 = 5;
constexpr int N4 = 7;

template <class Layout, class Space>
struct FillView_1D {
  using view_t = Kokkos::View<int*, Layout, Space>;
  view_t a;
  using policy_t = Kokkos::RangePolicy<typename Space::execution_space>;

  FillView_1D(view_t a_) : a(a_) {}

  void run() {
    Kokkos::parallel_for("FillView_1D", policy_t(0, a.extent(0)), *this);
  }
  KOKKOS_INLINE_FUNCTION
  void operator()(int i) const { a(i) = i; }
};

template <class Layout, class Space>
struct FillView_3D {
  using exec_t = typename Space::execution_space;
  using view_t = Kokkos::View<int***, Layout, Space>;
  using rank_t = Kokkos::Rank<
      view_t::rank,
      std::is_same<Layout, Kokkos::LayoutLeft>::value ? Kokkos::Iterate::Left
                                                      : Kokkos::Iterate::Right,
      std::is_same<Layout, Kokkos::LayoutLeft>::value ? Kokkos::Iterate::Left
                                                      : Kokkos::Iterate::Right>;
  using policy_t = Kokkos::MDRangePolicy<exec_t, rank_t>;

  view_t a;

  FillView_3D(view_t a_) : a(a_) {}

  void run() {
    Kokkos::parallel_for(
        "FillView_3D",
        policy_t({0, 0, 0}, {a.extent(0), a.extent(1), a.extent(2)}), *this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(int i0, int i1, int i2) const {
    a(i0, i1, i2) = 1000000 * i0 + 1000 * i1 + i2;
  }
};

template <class Layout, class Space>
struct FillView_4D {
  using exec_t = typename Space::execution_space;
  using view_t = Kokkos::View<int****, Layout, Space>;
  using rank_t = Kokkos::Rank<
      view_t::rank,
      std::is_same<Layout, Kokkos::LayoutLeft>::value ? Kokkos::Iterate::Left
                                                      : Kokkos::Iterate::Right,
      std::is_same<Layout, Kokkos::LayoutLeft>::value ? Kokkos::Iterate::Left
                                                      : Kokkos::Iterate::Right>;
  using policy_t = Kokkos::MDRangePolicy<exec_t, rank_t>;

  view_t a;

  FillView_4D(view_t a_) : a(a_) {}

  void run() {
    Kokkos::parallel_for("FillView_4D",
                         policy_t({0, 0, 0, 0}, {a.extent(0), a.extent(1),
                                                 a.extent(2), a.extent(3)}),
                         *this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(int i0, int i1, int i2, int i3) const {
    a(i0, i1, i2, i3) = 1000000 * i0 + 10000 * i1 + 100 * i2 + i3;
  }
};

template <class Layout, class Space>
struct FillView_5D {
  using exec_t = typename Space::execution_space;
  using view_t = Kokkos::View<int*****, Layout, Space>;
  using rank_t = Kokkos::Rank<
      view_t::rank,
      std::is_same<Layout, Kokkos::LayoutLeft>::value ? Kokkos::Iterate::Left
                                                      : Kokkos::Iterate::Right,
      std::is_same<Layout, Kokkos::LayoutLeft>::value ? Kokkos::Iterate::Left
                                                      : Kokkos::Iterate::Right>;
  using policy_t = Kokkos::MDRangePolicy<exec_t, rank_t>;

  view_t a;

  FillView_5D(view_t a_) : a(a_) {}

  void run() {
    Kokkos::parallel_for(
        "FillView_5D",
        policy_t({0, 0, 0, 0, 0}, {a.extent(0), a.extent(1), a.extent(2),
                                   a.extent(3), a.extent(4)}),
        *this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(int i0, int i1, int i2, int i3, int i4) const {
    a(i0, i1, i2, i3, i4) = 1000000 * i0 + 10000 * i1 + 100 * i2 + 10 * i3 + i4;
  }
};

template <class View, class SubView>
struct CheckSubviewCorrectness_1D_1D {
  using policy_t = Kokkos::RangePolicy<typename View::execution_space>;
  View a;
  SubView b;
  int offset;

  CheckSubviewCorrectness_1D_1D(View a_, SubView b_, int o)
      : a(a_), b(b_), offset(o) {}

  void run() {
    int errors = 0;
    Kokkos::parallel_reduce("CheckSubView_1D_1D", policy_t(0, b.size()), *this,
                            errors);
    ASSERT_EQ(errors, 0);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& i, int& e) const {
    if (a(i + offset) != b(i)) {
      e++;
    }
  }
};

template <class View, class SubView>
struct CheckSubviewCorrectness_1D_2D {
  using policy_t = Kokkos::RangePolicy<typename View::execution_space>;
  View a;
  SubView b;
  int i0;
  int offset;

  CheckSubviewCorrectness_1D_2D(View a_, SubView b_, int i0_, int o)
      : a(a_), b(b_), i0(i0_), offset(o) {}

  void run() {
    int errors = 0;
    Kokkos::parallel_reduce("CheckSubView_1D_2D", policy_t(0, b.size()), *this,
                            errors);
    ASSERT_EQ(errors, 0);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& i1, int& e) const {
    if (a(i0, i1 + offset) != b(i1)) {
      e++;
    }
  }
};

template <class View, class SubView>
struct CheckSubviewCorrectness_2D_3D {
  using policy_t = Kokkos::RangePolicy<typename View::execution_space>;
  using layout   = typename View::array_layout;
  View a;
  SubView b;
  int i0;
  int offset_1;
  int offset_2;

  CheckSubviewCorrectness_2D_3D(View a_, SubView b_, int i0_, int o1, int o2)
      : a(a_), b(b_), i0(i0_), offset_1(o1), offset_2(o2) {}

  void run() {
    int errors = 0;
    Kokkos::parallel_reduce("CheckSubView_2D_3D", policy_t(0, b.size()), *this,
                            errors);
    ASSERT_EQ(errors, 0);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& ii, int& e) const {
    const int i1 = std::is_same<layout, Kokkos::LayoutLeft>::value
                       ? ii % b.extent(0)
                       : ii / b.extent(1);

    const int i2 = std::is_same<layout, Kokkos::LayoutLeft>::value
                       ? ii / b.extent(0)
                       : ii % b.extent(1);

    if (a(i0, i1 + offset_1, i2 + offset_2) != b(i1, i2)) {
      e++;
    }
  }
};

template <class View, class SubView>
struct CheckSubviewCorrectness_3D_3D {
  using policy_t = Kokkos::RangePolicy<typename View::execution_space>;
  using layout   = typename View::array_layout;
  View a;
  SubView b;
  int offset_0;
  int offset_2;

  CheckSubviewCorrectness_3D_3D(View a_, SubView b_, int o0, int o2)
      : a(a_), b(b_), offset_0(o0), offset_2(o2) {}

  void run() {
    int errors = 0;
    Kokkos::parallel_reduce("CheckSubView_3D_3D", policy_t(0, b.size()), *this,
                            errors);
    ASSERT_EQ(errors, 0);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& ii, int& e) const {
    const int i0 = std::is_same<layout, Kokkos::LayoutLeft>::value
                       ? ii % b.extent(0)
                       : ii / (b.extent(1) * b.extent(2));

    const int i1 = std::is_same<layout, Kokkos::LayoutLeft>::value
                       ? (ii / b.extent(0)) % b.extent(1)
                       : (ii / b.extent(2)) % b.extent(1);

    const int i2 = std::is_same<layout, Kokkos::LayoutLeft>::value
                       ? ii / (b.extent(0) * b.extent(1))
                       : ii % b.extent(2);

    if (a(i0 + offset_0, i1, i2 + offset_2) != b(i0, i1, i2)) {
      e++;
    }
  }
};

template <class View, class SubView>
struct CheckSubviewCorrectness_3D_4D {
  using policy_t = Kokkos::RangePolicy<typename View::execution_space>;
  using layout   = typename View::array_layout;
  View a;
  SubView b;
  int index;
  int offset_0, offset_2;

  CheckSubviewCorrectness_3D_4D(View a_, SubView b_, int index_, int o0, int o2)
      : a(a_), b(b_), index(index_), offset_0(o0), offset_2(o2) {}

  void run() {
    int errors = 0;
    Kokkos::parallel_reduce("CheckSubView_3D_4D", policy_t(0, b.size()), *this,
                            errors);
    ASSERT_EQ(errors, 0);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& ii, int& e) const {
    const int i = std::is_same<layout, Kokkos::LayoutLeft>::value
                      ? ii % b.extent(0)
                      : ii / (b.extent(1) * b.extent(2));

    const int j = std::is_same<layout, Kokkos::LayoutLeft>::value
                      ? (ii / b.extent(0)) % b.extent(1)
                      : (ii / b.extent(2)) % b.extent(1);

    const int k = std::is_same<layout, Kokkos::LayoutLeft>::value
                      ? ii / (b.extent(0) * b.extent(1))
                      : ii % b.extent(2);

    int i0, i1, i2, i3;

    if (std::is_same<layout, Kokkos::LayoutLeft>::value) {
      i0 = i + offset_0;
      i1 = j;
      i2 = k + offset_2;
      i3 = index;
    } else {
      i0 = index;
      i1 = i + offset_0;
      i2 = j;
      i3 = k + offset_2;
    }

    if (a(i0, i1, i2, i3) != b(i, j, k)) e++;
  }
};

template <class View, class SubView>
struct CheckSubviewCorrectness_3D_5D {
  using policy_t = Kokkos::RangePolicy<typename View::execution_space>;
  using layout   = typename View::array_layout;
  View a;
  SubView b;
  int i0, i1;
  int offset_2, offset_3, offset_4;

  CheckSubviewCorrectness_3D_5D(View a_, SubView b_, int i0_, int i1_, int o2,
                                int o3, int o4)
      : a(a_),
        b(b_),
        i0(i0_),
        i1(i1_),
        offset_2(o2),
        offset_3(o3),
        offset_4(o4) {}

  void run() {
    int errors = 0;
    Kokkos::parallel_reduce("CheckSubView_3D_5D", policy_t(0, b.size()), *this,
                            errors);
    ASSERT_EQ(errors, 0);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& ii, int& e) const {
    const int i2 = std::is_same<layout, Kokkos::LayoutLeft>::value
                       ? ii % b.extent(0)
                       : ii / (b.extent(1) * b.extent(2));

    const int i3 = std::is_same<layout, Kokkos::LayoutLeft>::value
                       ? (ii / b.extent(0)) % b.extent(1)
                       : (ii / b.extent(2)) % b.extent(1);

    const int i4 = std::is_same<layout, Kokkos::LayoutLeft>::value
                       ? ii / (b.extent(0) * b.extent(1))
                       : ii % b.extent(2);

    if (a(i0, i1, i2 + offset_2, i3 + offset_3, i4 + offset_4) !=
        b(i2, i3, i4)) {
      e++;
    }
  }
};

template <class SubView, class View>
void test_Check1D(SubView a, View b, Kokkos::pair<int, int> range) {
  CheckSubviewCorrectness_1D_1D<View, SubView> check(b, a, range.first);
  check.run();
}

template <class SubView, class View>
void test_Check1D2D(SubView a, View b, int i0, std::pair<int, int> range) {
  CheckSubviewCorrectness_1D_2D<View, SubView> check(b, a, i0, range.first);
  check.run();
}

template <class SubView, class View>
void test_Check2D3D(SubView a, View b, int i0, std::pair<int, int> range1,
                    std::pair<int, int> range2) {
  CheckSubviewCorrectness_2D_3D<View, SubView> check(b, a, i0, range1.first,
                                                     range2.first);
  check.run();
}

template <class SubView, class View>
void test_Check3D5D(SubView a, View b, int i0, int i1,
                    Kokkos::pair<int, int> range2,
                    Kokkos::pair<int, int> range3,
                    Kokkos::pair<int, int> range4) {
  CheckSubviewCorrectness_3D_5D<View, SubView> check(
      b, a, i0, i1, range2.first, range3.first, range4.first);
  check.run();
}

template <class Space, class LayoutSub, class Layout, class LayoutOrg,
          class MemTraits>
void test_1d_assign_impl() {
  {  // Breaks.
    Kokkos::View<int*, LayoutOrg, Space> a_org("A", N0);
    Kokkos::View<int*, LayoutOrg, Space, MemTraits> a(a_org);
    Kokkos::fence();

    Impl::FillView_1D<LayoutOrg, Space> fill(a_org);
    fill.run();

    Kokkos::View<int[N0], Layout, Space, MemTraits> a1(a);
    Kokkos::fence();
    test_Check1D(a1, a, std::pair<int, int>(0, N0));

    Kokkos::View<int[N0], LayoutSub, Space, MemTraits> a2(a1);
    Kokkos::fence();
    test_Check1D(a2, a, std::pair<int, int>(0, N0));
    a1 = a;
    test_Check1D(a1, a, std::pair<int, int>(0, N0));

    // Runtime Fail expected.
    // Kokkos::View< int[N1] > afail1( a );

    // Compile Time Fail expected.
    // Kokkos::View< int[N1] > afail2( a1 );
  }

  {  // Works.
    Kokkos::View<int[N0], LayoutOrg, Space, MemTraits> a("A");
    Kokkos::View<int*, Layout, Space, MemTraits> a1(a);
    Kokkos::fence();
    test_Check1D(a1, a, std::pair<int, int>(0, N0));
    a1 = a;
    Kokkos::fence();
    test_Check1D(a1, a, std::pair<int, int>(0, N0));
  }
}

template <class Space, class Type, class TypeSub, class LayoutSub, class Layout,
          class LayoutOrg, class MemTraits>
void test_2d_subview_3d_impl_type() {
  Kokkos::View<int***, LayoutOrg, Space> a_org("A", N0, N1, N2);
  Kokkos::View<Type, Layout, Space, MemTraits> a(a_org);

  Impl::FillView_3D<LayoutOrg, Space> fill(a_org);
  fill.run();

  Kokkos::View<TypeSub, LayoutSub, Space, MemTraits> a1;
  a1 = Kokkos::subview(a, 3, Kokkos::ALL, Kokkos::ALL);
  Kokkos::fence();
  test_Check2D3D(a1, a, 3, std::pair<int, int>(0, N1),
                 std::pair<int, int>(0, N2));

  Kokkos::View<TypeSub, LayoutSub, Space, MemTraits> a2(a, 3, Kokkos::ALL,
                                                        Kokkos::ALL);
  Kokkos::fence();
  test_Check2D3D(a2, a, 3, std::pair<int, int>(0, N1),
                 std::pair<int, int>(0, N2));
}

template <class Space, class LayoutSub, class Layout, class LayoutOrg,
          class MemTraits>
void test_2d_subview_3d_impl_layout() {
  test_2d_subview_3d_impl_type<Space, int[N0][N1][N2], int[N1][N2], LayoutSub,
                               Layout, LayoutOrg, MemTraits>();
  test_2d_subview_3d_impl_type<Space, int[N0][N1][N2], int * [N2], LayoutSub,
                               Layout, LayoutOrg, MemTraits>();
  test_2d_subview_3d_impl_type<Space, int[N0][N1][N2], int**, LayoutSub, Layout,
                               LayoutOrg, MemTraits>();

  test_2d_subview_3d_impl_type<Space, int * [N1][N2], int[N1][N2], LayoutSub,
                               Layout, LayoutOrg, MemTraits>();
  test_2d_subview_3d_impl_type<Space, int * [N1][N2], int * [N2], LayoutSub,
                               Layout, LayoutOrg, MemTraits>();
  test_2d_subview_3d_impl_type<Space, int * [N1][N2], int**, LayoutSub, Layout,
                               LayoutOrg, MemTraits>();

  test_2d_subview_3d_impl_type<Space, int* * [N2], int[N1][N2], LayoutSub,
                               Layout, LayoutOrg, MemTraits>();
  test_2d_subview_3d_impl_type<Space, int* * [N2], int * [N2], LayoutSub,
                               Layout, LayoutOrg, MemTraits>();
  test_2d_subview_3d_impl_type<Space, int* * [N2], int**, LayoutSub, Layout,
                               LayoutOrg, MemTraits>();

  test_2d_subview_3d_impl_type<Space, int***, int[N1][N2], LayoutSub, Layout,
                               LayoutOrg, MemTraits>();
  test_2d_subview_3d_impl_type<Space, int***, int * [N2], LayoutSub, Layout,
                               LayoutOrg, MemTraits>();
  test_2d_subview_3d_impl_type<Space, int***, int**, LayoutSub, Layout,
                               LayoutOrg, MemTraits>();

  test_2d_subview_3d_impl_type<Space, const int[N0][N1][N2], const int[N1][N2],
                               LayoutSub, Layout, LayoutOrg, MemTraits>();
  test_2d_subview_3d_impl_type<Space, const int[N0][N1][N2], const int * [N2],
                               LayoutSub, Layout, LayoutOrg, MemTraits>();
  test_2d_subview_3d_impl_type<Space, const int[N0][N1][N2], const int**,
                               LayoutSub, Layout, LayoutOrg, MemTraits>();

  test_2d_subview_3d_impl_type<Space, const int * [N1][N2], const int[N1][N2],
                               LayoutSub, Layout, LayoutOrg, MemTraits>();
  test_2d_subview_3d_impl_type<Space, const int * [N1][N2], const int * [N2],
                               LayoutSub, Layout, LayoutOrg, MemTraits>();
  test_2d_subview_3d_impl_type<Space, const int * [N1][N2], const int**,
                               LayoutSub, Layout, LayoutOrg, MemTraits>();

  test_2d_subview_3d_impl_type<Space, const int* * [N2], const int[N1][N2],
                               LayoutSub, Layout, LayoutOrg, MemTraits>();
  test_2d_subview_3d_impl_type<Space, const int* * [N2], const int * [N2],
                               LayoutSub, Layout, LayoutOrg, MemTraits>();
  test_2d_subview_3d_impl_type<Space, const int* * [N2], const int**, LayoutSub,
                               Layout, LayoutOrg, MemTraits>();

  test_2d_subview_3d_impl_type<Space, const int***, const int[N1][N2],
                               LayoutSub, Layout, LayoutOrg, MemTraits>();
  test_2d_subview_3d_impl_type<Space, const int***, const int * [N2], LayoutSub,
                               Layout, LayoutOrg, MemTraits>();
  test_2d_subview_3d_impl_type<Space, const int***, const int**, LayoutSub,
                               Layout, LayoutOrg, MemTraits>();
}

template <class Space, class Type, class TypeSub, class LayoutSub, class Layout,
          class LayoutOrg, class MemTraits>
void test_3d_subview_5d_impl_type() {
  Kokkos::View<int*****, LayoutOrg, Space> a_org("A", N0, N1, N2, N3, N4);
  Kokkos::View<Type, Layout, Space, MemTraits> a(a_org);

  Impl::FillView_5D<LayoutOrg, Space> fill(a_org);
  fill.run();

  Kokkos::View<TypeSub, LayoutSub, Space, MemTraits> a1;
  a1 = Kokkos::subview(a, 3, 5, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
  Kokkos::fence();
  test_Check3D5D(a1, a, 3, 5, std::pair<int, int>(0, N2),
                 std::pair<int, int>(0, N3), std::pair<int, int>(0, N4));

  Kokkos::View<TypeSub, LayoutSub, Space, MemTraits> a2(
      a, 3, 5, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
  Kokkos::fence();
  test_Check3D5D(a2, a, 3, 5, std::pair<int, int>(0, N2),
                 std::pair<int, int>(0, N3), std::pair<int, int>(0, N4));
}

template <class Space, class LayoutSub, class Layout, class LayoutOrg,
          class MemTraits>
void test_3d_subview_5d_impl_layout() {
  test_3d_subview_5d_impl_type<Space, int[N0][N1][N2][N3][N4], int[N2][N3][N4],
                               LayoutSub, Layout, LayoutOrg, MemTraits>();
  test_3d_subview_5d_impl_type<Space, int[N0][N1][N2][N3][N4], int * [N3][N4],
                               LayoutSub, Layout, LayoutOrg, MemTraits>();
  test_3d_subview_5d_impl_type<Space, int[N0][N1][N2][N3][N4], int* * [N4],
                               LayoutSub, Layout, LayoutOrg, MemTraits>();
  test_3d_subview_5d_impl_type<Space, int[N0][N1][N2][N3][N4], int***,
                               LayoutSub, Layout, LayoutOrg, MemTraits>();

  test_3d_subview_5d_impl_type<Space, int * [N1][N2][N3][N4], int[N2][N3][N4],
                               LayoutSub, Layout, LayoutOrg, MemTraits>();
  test_3d_subview_5d_impl_type<Space, int * [N1][N2][N3][N4], int * [N3][N4],
                               LayoutSub, Layout, LayoutOrg, MemTraits>();
  test_3d_subview_5d_impl_type<Space, int * [N1][N2][N3][N4], int* * [N4],
                               LayoutSub, Layout, LayoutOrg, MemTraits>();
  test_3d_subview_5d_impl_type<Space, int * [N1][N2][N3][N4], int***, LayoutSub,
                               Layout, LayoutOrg, MemTraits>();

  test_3d_subview_5d_impl_type<Space, int* * [N2][N3][N4], int[N2][N3][N4],
                               LayoutSub, Layout, LayoutOrg, MemTraits>();
  test_3d_subview_5d_impl_type<Space, int* * [N2][N3][N4], int * [N3][N4],
                               LayoutSub, Layout, LayoutOrg, MemTraits>();
  test_3d_subview_5d_impl_type<Space, int* * [N2][N3][N4], int* * [N4],
                               LayoutSub, Layout, LayoutOrg, MemTraits>();
  test_3d_subview_5d_impl_type<Space, int* * [N2][N3][N4], int***, LayoutSub,
                               Layout, LayoutOrg, MemTraits>();

  test_3d_subview_5d_impl_type<Space, int** * [N3][N4], int[N2][N3][N4],
                               LayoutSub, Layout, LayoutOrg, MemTraits>();
  test_3d_subview_5d_impl_type<Space, int** * [N3][N4], int * [N3][N4],
                               LayoutSub, Layout, LayoutOrg, MemTraits>();
  test_3d_subview_5d_impl_type<Space, int** * [N3][N4], int* * [N4], LayoutSub,
                               Layout, LayoutOrg, MemTraits>();
  test_3d_subview_5d_impl_type<Space, int** * [N3][N4], int***, LayoutSub,
                               Layout, LayoutOrg, MemTraits>();

  test_3d_subview_5d_impl_type<Space, int*** * [N4], int[N2][N3][N4], LayoutSub,
                               Layout, LayoutOrg, MemTraits>();
  test_3d_subview_5d_impl_type<Space, int*** * [N4], int * [N3][N4], LayoutSub,
                               Layout, LayoutOrg, MemTraits>();
  test_3d_subview_5d_impl_type<Space, int*** * [N4], int* * [N4], LayoutSub,
                               Layout, LayoutOrg, MemTraits>();
  test_3d_subview_5d_impl_type<Space, int*** * [N4], int***, LayoutSub, Layout,
                               LayoutOrg, MemTraits>();

  test_3d_subview_5d_impl_type<Space, int*****, int[N2][N3][N4], LayoutSub,
                               Layout, LayoutOrg, MemTraits>();
  test_3d_subview_5d_impl_type<Space, int*****, int * [N3][N4], LayoutSub,
                               Layout, LayoutOrg, MemTraits>();
  test_3d_subview_5d_impl_type<Space, int*****, int* * [N4], LayoutSub, Layout,
                               LayoutOrg, MemTraits>();
  test_3d_subview_5d_impl_type<Space, int*****, int***, LayoutSub, Layout,
                               LayoutOrg, MemTraits>();

  test_3d_subview_5d_impl_type<Space, const int[N0][N1][N2][N3][N4],
                               const int[N2][N3][N4], LayoutSub, Layout,
                               LayoutOrg, MemTraits>();
  test_3d_subview_5d_impl_type<Space, const int[N0][N1][N2][N3][N4],
                               const int * [N3][N4], LayoutSub, Layout,
                               LayoutOrg, MemTraits>();
  test_3d_subview_5d_impl_type<Space, const int[N0][N1][N2][N3][N4],
                               const int* * [N4], LayoutSub, Layout, LayoutOrg,
                               MemTraits>();
  test_3d_subview_5d_impl_type<Space, const int[N0][N1][N2][N3][N4],
                               const int***, LayoutSub, Layout, LayoutOrg,
                               MemTraits>();

  test_3d_subview_5d_impl_type<Space, const int * [N1][N2][N3][N4],
                               const int[N2][N3][N4], LayoutSub, Layout,
                               LayoutOrg, MemTraits>();
  test_3d_subview_5d_impl_type<Space, const int * [N1][N2][N3][N4],
                               const int * [N3][N4], LayoutSub, Layout,
                               LayoutOrg, MemTraits>();
  test_3d_subview_5d_impl_type<Space, const int * [N1][N2][N3][N4],
                               const int* * [N4], LayoutSub, Layout, LayoutOrg,
                               MemTraits>();
  test_3d_subview_5d_impl_type<Space, const int * [N1][N2][N3][N4],
                               const int***, LayoutSub, Layout, LayoutOrg,
                               MemTraits>();

  test_3d_subview_5d_impl_type<Space, const int* * [N2][N3][N4],
                               const int[N2][N3][N4], LayoutSub, Layout,
                               LayoutOrg, MemTraits>();
  test_3d_subview_5d_impl_type<Space, const int* * [N2][N3][N4],
                               const int * [N3][N4], LayoutSub, Layout,
                               LayoutOrg, MemTraits>();
  test_3d_subview_5d_impl_type<Space, const int* * [N2][N3][N4],
                               const int* * [N4], LayoutSub, Layout, LayoutOrg,
                               MemTraits>();
  test_3d_subview_5d_impl_type<Space, const int* * [N2][N3][N4], const int***,
                               LayoutSub, Layout, LayoutOrg, MemTraits>();

  test_3d_subview_5d_impl_type<Space, const int** * [N3][N4],
                               const int[N2][N3][N4], LayoutSub, Layout,
                               LayoutOrg, MemTraits>();
  test_3d_subview_5d_impl_type<Space, const int** * [N3][N4],
                               const int * [N3][N4], LayoutSub, Layout,
                               LayoutOrg, MemTraits>();
  test_3d_subview_5d_impl_type<Space, const int** * [N3][N4], const int* * [N4],
                               LayoutSub, Layout, LayoutOrg, MemTraits>();
  test_3d_subview_5d_impl_type<Space, const int** * [N3][N4], const int***,
                               LayoutSub, Layout, LayoutOrg, MemTraits>();

  test_3d_subview_5d_impl_type<Space, const int*** * [N4],
                               const int[N2][N3][N4], LayoutSub, Layout,
                               LayoutOrg, MemTraits>();
  test_3d_subview_5d_impl_type<Space, const int*** * [N4], const int * [N3][N4],
                               LayoutSub, Layout, LayoutOrg, MemTraits>();
  test_3d_subview_5d_impl_type<Space, const int*** * [N4], const int* * [N4],
                               LayoutSub, Layout, LayoutOrg, MemTraits>();
  test_3d_subview_5d_impl_type<Space, const int*** * [N4], const int***,
                               LayoutSub, Layout, LayoutOrg, MemTraits>();

  test_3d_subview_5d_impl_type<Space, const int*****, const int[N2][N3][N4],
                               LayoutSub, Layout, LayoutOrg, MemTraits>();
  test_3d_subview_5d_impl_type<Space, const int*****, const int * [N3][N4],
                               LayoutSub, Layout, LayoutOrg, MemTraits>();
  test_3d_subview_5d_impl_type<Space, const int*****, const int* * [N4],
                               LayoutSub, Layout, LayoutOrg, MemTraits>();
  test_3d_subview_5d_impl_type<Space, const int*****, const int***, LayoutSub,
                               Layout, LayoutOrg, MemTraits>();
}

inline void test_subview_legal_args_right() {
  ASSERT_EQ(
      0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
             Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0, Kokkos::ALL_t,
             Kokkos::ALL_t, Kokkos::pair<int, int>, int, int>::value));
  ASSERT_EQ(
      0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
             Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0, Kokkos::ALL_t,
             Kokkos::ALL_t, Kokkos::ALL_t, int, int>::value));
  ASSERT_EQ(
      0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
             Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0, Kokkos::ALL_t,
             Kokkos::pair<int, int>, Kokkos::pair<int, int>, int, int>::value));
  ASSERT_EQ(
      0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
             Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0, Kokkos::ALL_t,
             Kokkos::pair<int, int>, Kokkos::ALL_t, int, int>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0,
                   Kokkos::pair<int, int>, Kokkos::ALL_t,
                   Kokkos::pair<int, int>, int, int>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0,
                   Kokkos::pair<int, int>, Kokkos::ALL_t, Kokkos::ALL_t, int,
                   int>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0,
                   Kokkos::pair<int, int>, Kokkos::pair<int, int>,
                   Kokkos::pair<int, int>, int, int>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0,
                   Kokkos::pair<int, int>, Kokkos::pair<int, int>,
                   Kokkos::ALL_t, int, int>::value));

  ASSERT_EQ(
      0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
             Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0, Kokkos::ALL_t,
             int, Kokkos::ALL_t, Kokkos::pair<int, int>, int>::value));
  ASSERT_EQ(
      0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
             Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0, Kokkos::ALL_t,
             int, Kokkos::ALL_t, Kokkos::ALL_t, int>::value));
  ASSERT_EQ(
      0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
             Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0, Kokkos::ALL_t,
             int, Kokkos::pair<int, int>, Kokkos::pair<int, int>, int>::value));
  ASSERT_EQ(
      0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
             Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0, Kokkos::ALL_t,
             int, Kokkos::pair<int, int>, Kokkos::ALL_t, int>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0,
                   Kokkos::pair<int, int>, int, Kokkos::ALL_t,
                   Kokkos::pair<int, int>, int>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0,
                   Kokkos::pair<int, int>, int, Kokkos::ALL_t, Kokkos::ALL_t,
                   int>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0,
                   Kokkos::pair<int, int>, int, Kokkos::pair<int, int>,
                   Kokkos::pair<int, int>, int>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0,
                   Kokkos::pair<int, int>, int, Kokkos::ALL_t,
                   Kokkos::pair<int, int>, int>::value));

  ASSERT_EQ(
      0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
             Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0, Kokkos::ALL_t,
             Kokkos::ALL_t, int, Kokkos::pair<int, int>, int>::value));
  ASSERT_EQ(
      0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
             Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0, Kokkos::ALL_t,
             Kokkos::ALL_t, int, Kokkos::ALL_t, int>::value));
  ASSERT_EQ(
      0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
             Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0, Kokkos::ALL_t,
             Kokkos::pair<int, int>, int, Kokkos::pair<int, int>, int>::value));
  ASSERT_EQ(
      0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
             Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0, Kokkos::ALL_t,
             Kokkos::pair<int, int>, int, Kokkos::ALL_t, int>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0,
                   Kokkos::pair<int, int>, Kokkos::ALL_t, int,
                   Kokkos::pair<int, int>, int>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0,
                   Kokkos::pair<int, int>, Kokkos::ALL_t, int, Kokkos::ALL_t,
                   int>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0,
                   Kokkos::pair<int, int>, Kokkos::pair<int, int>, int,
                   Kokkos::pair<int, int>, int>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0,
                   Kokkos::pair<int, int>, Kokkos::ALL_t, int,
                   Kokkos::pair<int, int>, int>::value));

  ASSERT_EQ(
      0,
      (Kokkos::Impl::SubviewLegalArgsCompileTime<
          Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0, int, Kokkos::ALL_t,
          Kokkos::ALL_t, Kokkos::pair<int, int>, int>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0, int,
                   Kokkos::ALL_t, Kokkos::ALL_t, Kokkos::ALL_t, int>::value));
  ASSERT_EQ(
      0,
      (Kokkos::Impl::SubviewLegalArgsCompileTime<
          Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0, int, Kokkos::ALL_t,
          Kokkos::pair<int, int>, Kokkos::pair<int, int>, int>::value));
  ASSERT_EQ(
      0,
      (Kokkos::Impl::SubviewLegalArgsCompileTime<
          Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0, int, Kokkos::ALL_t,
          Kokkos::pair<int, int>, Kokkos::ALL_t, int>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0, int,
                   Kokkos::pair<int, int>, Kokkos::ALL_t,
                   Kokkos::pair<int, int>, int>::value));
  ASSERT_EQ(
      0,
      (Kokkos::Impl::SubviewLegalArgsCompileTime<
          Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0, int,
          Kokkos::pair<int, int>, Kokkos::ALL_t, Kokkos::ALL_t, int>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0, int,
                   Kokkos::pair<int, int>, Kokkos::pair<int, int>,
                   Kokkos::pair<int, int>, int>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0, int,
                   Kokkos::pair<int, int>, Kokkos::pair<int, int>,
                   Kokkos::ALL_t, int>::value));

  ASSERT_EQ(
      0,
      (Kokkos::Impl::SubviewLegalArgsCompileTime<
          Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0, int, Kokkos::ALL_t,
          Kokkos::ALL_t, int, Kokkos::pair<int, int>>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0, int,
                   Kokkos::ALL_t, Kokkos::ALL_t, int, Kokkos::ALL_t>::value));
  ASSERT_EQ(
      0,
      (Kokkos::Impl::SubviewLegalArgsCompileTime<
          Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0, int, Kokkos::ALL_t,
          Kokkos::pair<int, int>, int, Kokkos::pair<int, int>>::value));
  ASSERT_EQ(
      0,
      (Kokkos::Impl::SubviewLegalArgsCompileTime<
          Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0, int, Kokkos::ALL_t,
          Kokkos::pair<int, int>, int, Kokkos::ALL_t>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0, int,
                   Kokkos::pair<int, int>, Kokkos::ALL_t, int,
                   Kokkos::pair<int, int>>::value));
  ASSERT_EQ(
      0,
      (Kokkos::Impl::SubviewLegalArgsCompileTime<
          Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0, int,
          Kokkos::pair<int, int>, Kokkos::ALL_t, int, Kokkos::ALL_t>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0, int,
                   Kokkos::pair<int, int>, Kokkos::pair<int, int>, int,
                   Kokkos::pair<int, int>>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0, int,
                   Kokkos::pair<int, int>, Kokkos::pair<int, int>, int,
                   Kokkos::ALL_t>::value));

  ASSERT_EQ(0,
            (Kokkos::Impl::SubviewLegalArgsCompileTime<
                Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0, int, int,
                Kokkos::ALL_t, Kokkos::ALL_t, Kokkos::pair<int, int>>::value));
  ASSERT_EQ(1, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0, int, int,
                   Kokkos::ALL_t, Kokkos::ALL_t, Kokkos::ALL_t>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0, int, int,
                   Kokkos::ALL_t, Kokkos::pair<int, int>,
                   Kokkos::pair<int, int>>::value));
  ASSERT_EQ(0,
            (Kokkos::Impl::SubviewLegalArgsCompileTime<
                Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0, int, int,
                Kokkos::ALL_t, Kokkos::pair<int, int>, Kokkos::ALL_t>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0, int, int,
                   Kokkos::pair<int, int>, Kokkos::ALL_t,
                   Kokkos::pair<int, int>>::value));
  ASSERT_EQ(1,
            (Kokkos::Impl::SubviewLegalArgsCompileTime<
                Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0, int, int,
                Kokkos::pair<int, int>, Kokkos::ALL_t, Kokkos::ALL_t>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0, int, int,
                   Kokkos::pair<int, int>, Kokkos::pair<int, int>,
                   Kokkos::pair<int, int>>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 5, 0, int, int,
                   Kokkos::pair<int, int>, Kokkos::pair<int, int>,
                   Kokkos::ALL_t>::value));

  ASSERT_EQ(1, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 3, 0,
                   Kokkos::ALL_t, Kokkos::ALL_t, Kokkos::ALL_t>::value));
  ASSERT_EQ(
      0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
             Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 3, 0, Kokkos::ALL_t,
             Kokkos::ALL_t, Kokkos::pair<int, int>>::value));
  ASSERT_EQ(1,
            (Kokkos::Impl::SubviewLegalArgsCompileTime<
                Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 3, 0,
                Kokkos::pair<int, int>, Kokkos::ALL_t, Kokkos::ALL_t>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 3, 0,
                   Kokkos::pair<int, int>, Kokkos::ALL_t,
                   Kokkos::pair<int, int>>::value));
  ASSERT_EQ(
      0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
             Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 3, 0, Kokkos::ALL_t,
             Kokkos::pair<int, int>, Kokkos::ALL_t>::value));
  ASSERT_EQ(
      0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
             Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 3, 0, Kokkos::ALL_t,
             Kokkos::pair<int, int>, Kokkos::pair<int, int>>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 3, 0,
                   Kokkos::pair<int, int>, Kokkos::pair<int, int>,
                   Kokkos::ALL_t>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutRight, Kokkos::LayoutRight, 3, 3, 0,
                   Kokkos::pair<int, int>, Kokkos::pair<int, int>,
                   Kokkos::pair<int, int>>::value));
}

inline void test_subview_legal_args_left() {
  ASSERT_EQ(1,
            (Kokkos::Impl::SubviewLegalArgsCompileTime<
                Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0, Kokkos::ALL_t,
                Kokkos::ALL_t, Kokkos::pair<int, int>, int, int>::value));
  ASSERT_EQ(1,
            (Kokkos::Impl::SubviewLegalArgsCompileTime<
                Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0, Kokkos::ALL_t,
                Kokkos::ALL_t, Kokkos::ALL_t, int, int>::value));
  ASSERT_EQ(
      0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
             Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0, Kokkos::ALL_t,
             Kokkos::pair<int, int>, Kokkos::pair<int, int>, int, int>::value));
  ASSERT_EQ(0,
            (Kokkos::Impl::SubviewLegalArgsCompileTime<
                Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0, Kokkos::ALL_t,
                Kokkos::pair<int, int>, Kokkos::ALL_t, int, int>::value));
  ASSERT_EQ(1, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0,
                   Kokkos::pair<int, int>, Kokkos::ALL_t,
                   Kokkos::pair<int, int>, int, int>::value));
  ASSERT_EQ(1, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0,
                   Kokkos::pair<int, int>, Kokkos::ALL_t, Kokkos::ALL_t, int,
                   int>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0,
                   Kokkos::pair<int, int>, Kokkos::pair<int, int>,
                   Kokkos::pair<int, int>, int, int>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0,
                   Kokkos::pair<int, int>, Kokkos::pair<int, int>,
                   Kokkos::ALL_t, int, int>::value));

  ASSERT_EQ(0,
            (Kokkos::Impl::SubviewLegalArgsCompileTime<
                Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0, Kokkos::ALL_t,
                int, Kokkos::ALL_t, Kokkos::pair<int, int>, int>::value));
  ASSERT_EQ(0,
            (Kokkos::Impl::SubviewLegalArgsCompileTime<
                Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0, Kokkos::ALL_t,
                int, Kokkos::ALL_t, Kokkos::ALL_t, int>::value));
  ASSERT_EQ(
      0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
             Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0, Kokkos::ALL_t,
             int, Kokkos::pair<int, int>, Kokkos::pair<int, int>, int>::value));
  ASSERT_EQ(0,
            (Kokkos::Impl::SubviewLegalArgsCompileTime<
                Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0, Kokkos::ALL_t,
                int, Kokkos::pair<int, int>, Kokkos::ALL_t, int>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0,
                   Kokkos::pair<int, int>, int, Kokkos::ALL_t,
                   Kokkos::pair<int, int>, int>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0,
                   Kokkos::pair<int, int>, int, Kokkos::ALL_t, Kokkos::ALL_t,
                   int>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0,
                   Kokkos::pair<int, int>, int, Kokkos::pair<int, int>,
                   Kokkos::pair<int, int>, int>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0,
                   Kokkos::pair<int, int>, int, Kokkos::ALL_t,
                   Kokkos::pair<int, int>, int>::value));

  ASSERT_EQ(0,
            (Kokkos::Impl::SubviewLegalArgsCompileTime<
                Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0, Kokkos::ALL_t,
                Kokkos::ALL_t, int, Kokkos::pair<int, int>, int>::value));
  ASSERT_EQ(0,
            (Kokkos::Impl::SubviewLegalArgsCompileTime<
                Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0, Kokkos::ALL_t,
                Kokkos::ALL_t, int, Kokkos::ALL_t, int>::value));
  ASSERT_EQ(
      0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
             Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0, Kokkos::ALL_t,
             Kokkos::pair<int, int>, int, Kokkos::pair<int, int>, int>::value));
  ASSERT_EQ(0,
            (Kokkos::Impl::SubviewLegalArgsCompileTime<
                Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0, Kokkos::ALL_t,
                Kokkos::pair<int, int>, int, Kokkos::ALL_t, int>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0,
                   Kokkos::pair<int, int>, Kokkos::ALL_t, int,
                   Kokkos::pair<int, int>, int>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0,
                   Kokkos::pair<int, int>, Kokkos::ALL_t, int, Kokkos::ALL_t,
                   int>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0,
                   Kokkos::pair<int, int>, Kokkos::pair<int, int>, int,
                   Kokkos::pair<int, int>, int>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0,
                   Kokkos::pair<int, int>, Kokkos::ALL_t, int,
                   Kokkos::pair<int, int>, int>::value));

  ASSERT_EQ(
      0,
      (Kokkos::Impl::SubviewLegalArgsCompileTime<
          Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0, int, Kokkos::ALL_t,
          Kokkos::ALL_t, Kokkos::pair<int, int>, int>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0, int,
                   Kokkos::ALL_t, Kokkos::ALL_t, Kokkos::ALL_t, int>::value));
  ASSERT_EQ(
      0,
      (Kokkos::Impl::SubviewLegalArgsCompileTime<
          Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0, int, Kokkos::ALL_t,
          Kokkos::pair<int, int>, Kokkos::pair<int, int>, int>::value));
  ASSERT_EQ(
      0,
      (Kokkos::Impl::SubviewLegalArgsCompileTime<
          Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0, int, Kokkos::ALL_t,
          Kokkos::pair<int, int>, Kokkos::ALL_t, int>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0, int,
                   Kokkos::pair<int, int>, Kokkos::ALL_t,
                   Kokkos::pair<int, int>, int>::value));
  ASSERT_EQ(
      0,
      (Kokkos::Impl::SubviewLegalArgsCompileTime<
          Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0, int,
          Kokkos::pair<int, int>, Kokkos::ALL_t, Kokkos::ALL_t, int>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0, int,
                   Kokkos::pair<int, int>, Kokkos::pair<int, int>,
                   Kokkos::pair<int, int>, int>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0, int,
                   Kokkos::pair<int, int>, Kokkos::pair<int, int>,
                   Kokkos::ALL_t, int>::value));

  ASSERT_EQ(
      0,
      (Kokkos::Impl::SubviewLegalArgsCompileTime<
          Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0, int, Kokkos::ALL_t,
          Kokkos::ALL_t, int, Kokkos::pair<int, int>>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0, int,
                   Kokkos::ALL_t, Kokkos::ALL_t, int, Kokkos::ALL_t>::value));
  ASSERT_EQ(
      0,
      (Kokkos::Impl::SubviewLegalArgsCompileTime<
          Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0, int, Kokkos::ALL_t,
          Kokkos::pair<int, int>, int, Kokkos::pair<int, int>>::value));
  ASSERT_EQ(
      0,
      (Kokkos::Impl::SubviewLegalArgsCompileTime<
          Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0, int, Kokkos::ALL_t,
          Kokkos::pair<int, int>, int, Kokkos::ALL_t>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0, int,
                   Kokkos::pair<int, int>, Kokkos::ALL_t, int,
                   Kokkos::pair<int, int>>::value));
  ASSERT_EQ(
      0,
      (Kokkos::Impl::SubviewLegalArgsCompileTime<
          Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0, int,
          Kokkos::pair<int, int>, Kokkos::ALL_t, int, Kokkos::ALL_t>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0, int,
                   Kokkos::pair<int, int>, Kokkos::pair<int, int>, int,
                   Kokkos::pair<int, int>>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0, int,
                   Kokkos::pair<int, int>, Kokkos::pair<int, int>, int,
                   Kokkos::ALL_t>::value));

  ASSERT_EQ(0,
            (Kokkos::Impl::SubviewLegalArgsCompileTime<
                Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0, int, int,
                Kokkos::ALL_t, Kokkos::ALL_t, Kokkos::pair<int, int>>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0, int, int,
                   Kokkos::ALL_t, Kokkos::ALL_t, Kokkos::ALL_t>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0, int, int,
                   Kokkos::ALL_t, Kokkos::pair<int, int>,
                   Kokkos::pair<int, int>>::value));
  ASSERT_EQ(0,
            (Kokkos::Impl::SubviewLegalArgsCompileTime<
                Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0, int, int,
                Kokkos::ALL_t, Kokkos::pair<int, int>, Kokkos::ALL_t>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0, int, int,
                   Kokkos::pair<int, int>, Kokkos::ALL_t,
                   Kokkos::pair<int, int>>::value));
  ASSERT_EQ(0,
            (Kokkos::Impl::SubviewLegalArgsCompileTime<
                Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0, int, int,
                Kokkos::pair<int, int>, Kokkos::ALL_t, Kokkos::ALL_t>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0, int, int,
                   Kokkos::pair<int, int>, Kokkos::pair<int, int>,
                   Kokkos::pair<int, int>>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 5, 0, int, int,
                   Kokkos::pair<int, int>, Kokkos::pair<int, int>,
                   Kokkos::ALL_t>::value));

  ASSERT_EQ(1,
            (Kokkos::Impl::SubviewLegalArgsCompileTime<
                Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 3, 0, Kokkos::ALL_t,
                Kokkos::ALL_t, Kokkos::pair<int, int>>::value));
  ASSERT_EQ(1, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 3, 0,
                   Kokkos::ALL_t, Kokkos::ALL_t, Kokkos::ALL_t>::value));
  ASSERT_EQ(1, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 3, 0,
                   Kokkos::pair<int, int>, Kokkos::ALL_t,
                   Kokkos::pair<int, int>>::value));
  ASSERT_EQ(1,
            (Kokkos::Impl::SubviewLegalArgsCompileTime<
                Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 3, 0,
                Kokkos::pair<int, int>, Kokkos::ALL_t, Kokkos::ALL_t>::value));
  ASSERT_EQ(0,
            (Kokkos::Impl::SubviewLegalArgsCompileTime<
                Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 3, 0, Kokkos::ALL_t,
                Kokkos::pair<int, int>, Kokkos::ALL_t>::value));
  ASSERT_EQ(0,
            (Kokkos::Impl::SubviewLegalArgsCompileTime<
                Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 3, 0, Kokkos::ALL_t,
                Kokkos::pair<int, int>, Kokkos::pair<int, int>>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 3, 0,
                   Kokkos::pair<int, int>, Kokkos::pair<int, int>,
                   Kokkos::ALL_t>::value));
  ASSERT_EQ(0, (Kokkos::Impl::SubviewLegalArgsCompileTime<
                   Kokkos::LayoutLeft, Kokkos::LayoutLeft, 3, 3, 0,
                   Kokkos::pair<int, int>, Kokkos::pair<int, int>,
                   Kokkos::pair<int, int>>::value));
}

}  // namespace Impl

template <class Space, class MemTraits = void>
void test_1d_assign() {
  Impl::test_1d_assign_impl<Space, Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                            Kokkos::LayoutLeft, MemTraits>();
  // Impl::test_1d_assign_impl< Space, Kokkos::LayoutRight, Kokkos::LayoutLeft,
  // Kokkos::LayoutLeft >();
  Impl::test_1d_assign_impl<Space, Kokkos::LayoutStride, Kokkos::LayoutLeft,
                            Kokkos::LayoutLeft, MemTraits>();
  // Impl::test_1d_assign_impl< Space, Kokkos::LayoutLeft, Kokkos::LayoutRight,
  // Kokkos::LayoutLeft >();
  Impl::test_1d_assign_impl<Space, Kokkos::LayoutRight, Kokkos::LayoutRight,
                            Kokkos::LayoutRight, MemTraits>();
  Impl::test_1d_assign_impl<Space, Kokkos::LayoutStride, Kokkos::LayoutRight,
                            Kokkos::LayoutRight, MemTraits>();
  // Impl::test_1d_assign_impl< Space, Kokkos::LayoutLeft, Kokkos::LayoutStride,
  // Kokkos::LayoutLeft >(); Impl::test_1d_assign_impl< Space,
  // Kokkos::LayoutRight, Kokkos::LayoutStride, Kokkos::LayoutLeft >();
  Impl::test_1d_assign_impl<Space, Kokkos::LayoutStride, Kokkos::LayoutStride,
                            Kokkos::LayoutLeft, MemTraits>();
}

template <class Space, class MemTraits = void>
void test_2d_subview_3d() {
  Impl::test_2d_subview_3d_impl_layout<Space, Kokkos::LayoutRight,
                                       Kokkos::LayoutRight, Kokkos::LayoutRight,
                                       MemTraits>();
  Impl::test_2d_subview_3d_impl_layout<Space, Kokkos::LayoutStride,
                                       Kokkos::LayoutRight, Kokkos::LayoutRight,
                                       MemTraits>();
  Impl::test_2d_subview_3d_impl_layout<Space, Kokkos::LayoutStride,
                                       Kokkos::LayoutStride,
                                       Kokkos::LayoutRight, MemTraits>();
  Impl::test_2d_subview_3d_impl_layout<Space, Kokkos::LayoutStride,
                                       Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                       MemTraits>();
  Impl::test_2d_subview_3d_impl_layout<Space, Kokkos::LayoutStride,
                                       Kokkos::LayoutStride, Kokkos::LayoutLeft,
                                       MemTraits>();
}

template <class Space, class MemTraits = void>
void test_3d_subview_5d_right() {
  Impl::test_3d_subview_5d_impl_layout<Space, Kokkos::LayoutStride,
                                       Kokkos::LayoutRight, Kokkos::LayoutRight,
                                       MemTraits>();
  Impl::test_3d_subview_5d_impl_layout<Space, Kokkos::LayoutStride,
                                       Kokkos::LayoutStride,
                                       Kokkos::LayoutRight, MemTraits>();
}

template <class Space, class MemTraits = void>
void test_3d_subview_5d_left() {
  Impl::test_3d_subview_5d_impl_layout<Space, Kokkos::LayoutStride,
                                       Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                       MemTraits>();
  Impl::test_3d_subview_5d_impl_layout<Space, Kokkos::LayoutStride,
                                       Kokkos::LayoutStride, Kokkos::LayoutLeft,
                                       MemTraits>();
}

template <class Space, class MemTraits = void>
void test_layoutleft_to_layoutleft() {
  Impl::test_subview_legal_args_left();

  using view3D_t = Kokkos::View<int***, Kokkos::LayoutLeft, Space>;
  using view4D_t = Kokkos::View<int****, Kokkos::LayoutLeft, Space>;
  {
    view3D_t a("A", 100, 4, 3);
    view3D_t b(a, Kokkos::pair<int, int>(16, 32), Kokkos::ALL, Kokkos::ALL);

    Impl::FillView_3D<Kokkos::LayoutLeft, Space> fill(a);
    fill.run();

    Impl::CheckSubviewCorrectness_3D_3D<view3D_t, view3D_t> check(a, b, 16, 0);
    check.run();
  }

  {
    view3D_t a("A", 100, 4, 5);
    view3D_t b(a, Kokkos::pair<int, int>(16, 32), Kokkos::ALL,
               Kokkos::pair<int, int>(1, 3));

    Impl::FillView_3D<Kokkos::LayoutLeft, Space> fill(a);
    fill.run();

    Impl::CheckSubviewCorrectness_3D_3D<view3D_t, view3D_t> check(a, b, 16, 1);
    check.run();
  }

  {
    view4D_t a("A", 100, 4, 5, 3);
    view3D_t b(a, Kokkos::pair<int, int>(16, 32), Kokkos::ALL,
               Kokkos::pair<int, int>(1, 3), 1);

    Impl::FillView_4D<Kokkos::LayoutLeft, Space> fill(a);
    fill.run();

    Impl::CheckSubviewCorrectness_3D_4D<view4D_t, view3D_t> check(a, b, 1, 16,
                                                                  1);
    check.run();
  }
}

template <class Space, class MemTraits = void>
void test_layoutright_to_layoutright() {
  Impl::test_subview_legal_args_right();

  using view3D_t = Kokkos::View<int***, Kokkos::LayoutRight, Space>;
  using view4D_t = Kokkos::View<int****, Kokkos::LayoutRight, Space>;
  {
    view3D_t a("A", 100, 4, 3);
    view3D_t b(a, Kokkos::pair<int, int>(16, 32), Kokkos::ALL, Kokkos::ALL);

    Impl::FillView_3D<Kokkos::LayoutRight, Space> fill(a);
    fill.run();

    Impl::CheckSubviewCorrectness_3D_3D<view3D_t, view3D_t> check(a, b, 16, 0);
    check.run();
  }
  {
    view4D_t a("A", 3, 4, 5, 100);
    view3D_t b(a, 1, Kokkos::pair<int, int>(1, 3), Kokkos::ALL, Kokkos::ALL);

    Impl::FillView_4D<Kokkos::LayoutRight, Space> fill(a);
    fill.run();

    Impl::CheckSubviewCorrectness_3D_4D<view4D_t, view3D_t> check(a, b, 1, 1,
                                                                  0);
    check.run();
  }
}
//----------------------------------------------------------------------------

template <class Space>
struct TestUnmanagedSubviewReset {
  Kokkos::View<int****, Space> a;

  KOKKOS_INLINE_FUNCTION
  void operator()(int) const noexcept {
    auto sub_a = Kokkos::subview(a, 0, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);

    for (int i = 0; i < int(a.extent(0)); ++i) {
      sub_a.assign_data(&a(i, 0, 0, 0));
      if (&sub_a(1, 1, 1) != &a(i, 1, 1, 1)) {
        Kokkos::abort("TestUnmanagedSubviewReset");
      }
    }
  }

  TestUnmanagedSubviewReset() : a(Kokkos::view_alloc(), 20, 10, 5, 2) {}
};

template <class Space>
void test_unmanaged_subview_reset() {
  Kokkos::parallel_for(
      Kokkos::RangePolicy<typename Space::execution_space>(0, 1),
      TestUnmanagedSubviewReset<Space>());
}

//----------------------------------------------------------------------------

template <std::underlying_type_t<Kokkos::MemoryTraitsFlags> MTF>
struct TestSubviewMemoryTraitsConstruction {
  void operator()() const noexcept {
    using memory_traits_type = Kokkos::MemoryTraits<MTF>;
    using view_type =
        Kokkos::View<double*, Kokkos::HostSpace, memory_traits_type>;
    using size_type = typename view_type::size_type;

    // Create a managed View first and then apply the desired memory traits to
    // an unmanaged version of it since a managed View can't use the Unmanaged
    // trait.
    Kokkos::View<double*, Kokkos::HostSpace> v_original("v", 7);
    view_type v(v_original.data(), v_original.size());
    for (size_type i = 0; i != v.size(); ++i) v[i] = static_cast<double>(i);

    std::pair<int, int> range(3, 5);
    auto sv = Kokkos::subview(v, range);

    // check that the subview memory traits are the same as the original view
    // (with the Aligned trait stripped).
    using view_memory_traits    = typename decltype(v)::memory_traits;
    using subview_memory_traits = typename decltype(sv)::memory_traits;
    static_assert(view_memory_traits::impl_value ==
                  memory_traits_type::impl_value);
    if constexpr (memory_traits_type::is_aligned)
      static_assert(subview_memory_traits::impl_value + Kokkos::Aligned ==
                    memory_traits_type::impl_value);
    else
      static_assert(subview_memory_traits::impl_value ==
                    memory_traits_type::impl_value);

    ASSERT_EQ(2u, sv.size());
    EXPECT_EQ(3., sv[0]);
    EXPECT_EQ(4., sv[1]);
  }
};

inline void test_subview_memory_traits_construction() {
  // Test all combinations of MemoryTraits:
  // Unmanaged (1)
  // RandomAccess (2)
  // Atomic (4)
  // Restricted (8)
  // Aligned (16)
  TestSubviewMemoryTraitsConstruction<0>()();
  TestSubviewMemoryTraitsConstruction<1>()();
  TestSubviewMemoryTraitsConstruction<2>()();
  TestSubviewMemoryTraitsConstruction<3>()();
  TestSubviewMemoryTraitsConstruction<4>()();
  TestSubviewMemoryTraitsConstruction<5>()();
  TestSubviewMemoryTraitsConstruction<6>()();
  TestSubviewMemoryTraitsConstruction<7>()();
  TestSubviewMemoryTraitsConstruction<8>()();
  TestSubviewMemoryTraitsConstruction<9>()();
  TestSubviewMemoryTraitsConstruction<10>()();
  TestSubviewMemoryTraitsConstruction<11>()();
  TestSubviewMemoryTraitsConstruction<12>()();
  TestSubviewMemoryTraitsConstruction<13>()();
  TestSubviewMemoryTraitsConstruction<14>()();
  TestSubviewMemoryTraitsConstruction<15>()();
  TestSubviewMemoryTraitsConstruction<16>()();
  TestSubviewMemoryTraitsConstruction<17>()();
  TestSubviewMemoryTraitsConstruction<18>()();
  TestSubviewMemoryTraitsConstruction<19>()();
  TestSubviewMemoryTraitsConstruction<20>()();
  TestSubviewMemoryTraitsConstruction<21>()();
  TestSubviewMemoryTraitsConstruction<22>()();
  TestSubviewMemoryTraitsConstruction<23>()();
  TestSubviewMemoryTraitsConstruction<24>()();
  TestSubviewMemoryTraitsConstruction<25>()();
  TestSubviewMemoryTraitsConstruction<26>()();
  TestSubviewMemoryTraitsConstruction<27>()();
  TestSubviewMemoryTraitsConstruction<28>()();
  TestSubviewMemoryTraitsConstruction<29>()();
  TestSubviewMemoryTraitsConstruction<30>()();
  TestSubviewMemoryTraitsConstruction<31>()();
}

//----------------------------------------------------------------------------

template <class T>
struct get_view_type;

template <class T, class... Args>
struct get_view_type<Kokkos::View<T, Args...>> {
  using type = T;
};

template <class T>
struct
    ___________________________________TYPE_DISPLAY________________________________________;
#define TYPE_DISPLAY(...)                                                                           \
  typename ___________________________________TYPE_DISPLAY________________________________________< \
      __VA_ARGS__>::type notdefined;

template <class Space, class Layout>
struct TestSubviewStaticSizes {
  Kokkos::View<int * [10][5][2], Layout, Space> a;
  Kokkos::View<int[6][7][8], Layout, Space> b;

  KOKKOS_INLINE_FUNCTION
  int operator()() const noexcept {
    /* Doesn't actually do anything; just static assertions */

    auto sub_a = Kokkos::subview(a, 0, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
    typename static_expect_same<
        /* expected */ int[10][5][2],
        /*  actual  */ typename get_view_type<decltype(sub_a)>::type>::type
        test_1 = 0;

    auto sub_a_2 = Kokkos::subview(a, 0, 0, Kokkos::ALL, Kokkos::ALL);
    typename static_expect_same<
        /* expected */ int[5][2],
        /*  actual  */ typename get_view_type<decltype(sub_a_2)>::type>::type
        test_2 = 0;

    auto sub_a_3 = Kokkos::subview(a, 0, 0, Kokkos::ALL, 0);
    typename static_expect_same<
        /* expected */ int[5],
        /*  actual  */ typename get_view_type<decltype(sub_a_3)>::type>::type
        test_3 = 0;

    auto sub_a_4 = Kokkos::subview(a, Kokkos::ALL, 0, Kokkos::ALL, Kokkos::ALL);
    typename static_expect_same<
        /* expected */ int * [5][2],
        /*  actual  */ typename get_view_type<decltype(sub_a_4)>::type>::type
        test_4 = 0;

    // TODO we'll need to update this test once we allow interleaving of static
    // and dynamic
    auto sub_a_5 = Kokkos::subview(a, Kokkos::ALL, 0, Kokkos::ALL,
                                   Kokkos::make_pair(0, 1));
    typename static_expect_same<
        /* expected */ int***,
        /*  actual  */ typename get_view_type<decltype(sub_a_5)>::type>::type
        test_5 = 0;

    auto sub_a_sub = Kokkos::subview(sub_a_5, 0, Kokkos::ALL, 0);
    typename static_expect_same<
        /* expected */ int*,
        /*  actual  */ typename get_view_type<decltype(sub_a_sub)>::type>::type
        test_sub = 0;

    auto sub_a_7 = Kokkos::subview(a, Kokkos::ALL, 0, Kokkos::make_pair(0, 1),
                                   Kokkos::ALL);
    typename static_expect_same<
        /* expected */ int* * [2],
        /*  actual  */ typename get_view_type<decltype(sub_a_7)>::type>::type
        test_7 = 0;

    auto sub_a_8 =
        Kokkos::subview(a, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
    typename static_expect_same<
        /* expected */ int * [10][5][2],
        /*  actual  */ typename get_view_type<decltype(sub_a_8)>::type>::type
        test_8 = 0;

    auto sub_b = Kokkos::subview(b, Kokkos::ALL, Kokkos::ALL, Kokkos::ALL);
    typename static_expect_same<
        /* expected */ int[6][7][8],
        /*  actual  */ typename get_view_type<decltype(sub_b)>::type>::type
        test_9 = 0;

    auto sub_b_2 = Kokkos::subview(b, 0, Kokkos::ALL, Kokkos::ALL);
    typename static_expect_same<
        /* expected */ int[7][8],
        /*  actual  */ typename get_view_type<decltype(sub_b_2)>::type>::type
        test_10 = 0;

    auto sub_b_3 =
        Kokkos::subview(b, Kokkos::make_pair(2, 3), Kokkos::ALL, Kokkos::ALL);
    typename static_expect_same<
        /* expected */ int * [7][8],
        /*  actual  */ typename get_view_type<decltype(sub_b_3)>::type>::type
        test_11 = 0;

    return test_1 + test_2 + test_3 + test_4 + test_5 + test_sub + test_7 +
           test_8 + test_9 + test_10 + test_11;
  }

  TestSubviewStaticSizes() : a(Kokkos::view_alloc("a"), 20), b("b") {}
};

template <class Space>
struct TestExtentsStaticTests {
  using test1 = typename static_expect_same<
      /* expected */
      Kokkos::Experimental::Extents<Kokkos::Experimental::dynamic_extent,
                                    Kokkos::Experimental::dynamic_extent, 1, 2,
                                    3>,
      /* actual */
      typename Kokkos::Impl::ParseViewExtents<double* * [1][2][3]>::type>::type;

  using test2 = typename static_expect_same<
      /* expected */
      Kokkos::Experimental::Extents<1, 2, 3>,
      /* actual */
      typename Kokkos::Impl::ParseViewExtents<double[1][2][3]>::type>::type;

  using test3 = typename static_expect_same<
      /* expected */
      Kokkos::Experimental::Extents<3>,
      /* actual */
      typename Kokkos::Impl::ParseViewExtents<double[3]>::type>::type;

  using test4 = typename static_expect_same<
      /* expected */
      Kokkos::Experimental::Extents<>,
      /* actual */
      typename Kokkos::Impl::ParseViewExtents<double>::type>::type;
};

}  // namespace TestViewSubview

#endif
