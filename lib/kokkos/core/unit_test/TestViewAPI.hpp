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
#include <sstream>
#include <iostream>

namespace Test {

template <class T, class... P>
size_t allocation_count(const Kokkos::View<T, P...> &view) {
  const size_t card  = view.size();
  const size_t alloc = view.span();

  const int memory_span = Kokkos::View<int *>::required_allocation_size(100);
  const int memory_span_layout =
      Kokkos::View<int *, Kokkos::LayoutRight>::required_allocation_size(
          Kokkos::LayoutRight(100));

  return ((card <= alloc) && (memory_span == 400) &&
          (memory_span_layout == 400))
             ? alloc
             : 0;
}

/*--------------------------------------------------------------------------*/

template <typename T, class DeviceType>
struct TestViewOperator {
  using execution_space = typename DeviceType::execution_space;

  enum { N = 1000 };
  enum { D = 3 };

  using view_type = Kokkos::View<T * [D], execution_space>;

  const view_type v1;
  const view_type v2;

  TestViewOperator() : v1("v1", N), v2("v2", N) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const unsigned i) const {
    const unsigned X = 0;
    const unsigned Y = 1;
    const unsigned Z = 2;

    v2(i, X) = v1(i, X);
    v2(i, Y) = v1(i, Y);
    v2(i, Z) = v1(i, Z);
  }
};

/*--------------------------------------------------------------------------*/

template <class DataType, class DeviceType,
          unsigned Rank = Kokkos::ViewTraits<DataType>::rank>
struct TestViewOperator_LeftAndRight;

template <class DataType, class DeviceType>
struct TestViewOperator_LeftAndRight<DataType, DeviceType, 8> {
  using execution_space = typename DeviceType::execution_space;
  using memory_space    = typename DeviceType::memory_space;
  using size_type       = typename execution_space::size_type;

  using value_type = int;

  KOKKOS_INLINE_FUNCTION
  static void join(value_type &update, const value_type &input) {
    update |= input;
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type &update) { update = 0; }

  using left_view = Kokkos::View<DataType, Kokkos::LayoutLeft, execution_space>;
  using right_view =
      Kokkos::View<DataType, Kokkos::LayoutRight, execution_space>;
  using stride_view =
      Kokkos::View<DataType, Kokkos::LayoutStride, execution_space>;

  left_view left;
  right_view right;
  stride_view left_stride;
  stride_view right_stride;
  long left_alloc;
  long right_alloc;

  TestViewOperator_LeftAndRight()
      : left("left"),
        right("right"),
        left_stride(left),
        right_stride(right),
        left_alloc(allocation_count(left)),
        right_alloc(allocation_count(right)) {}

  void testit() {
    int error_flag = 0;

    Kokkos::parallel_reduce(1, *this, error_flag);

    ASSERT_EQ(error_flag, 0);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_type, value_type &update) const {
    long offset = -1;

    for (unsigned i7 = 0; i7 < unsigned(left.extent(7)); ++i7)
      for (unsigned i6 = 0; i6 < unsigned(left.extent(6)); ++i6)
        for (unsigned i5 = 0; i5 < unsigned(left.extent(5)); ++i5)
          for (unsigned i4 = 0; i4 < unsigned(left.extent(4)); ++i4)
            for (unsigned i3 = 0; i3 < unsigned(left.extent(3)); ++i3)
              for (unsigned i2 = 0; i2 < unsigned(left.extent(2)); ++i2)
                for (unsigned i1 = 0; i1 < unsigned(left.extent(1)); ++i1)
                  for (unsigned i0 = 0; i0 < unsigned(left.extent(0)); ++i0) {
                    const long j = &left(i0, i1, i2, i3, i4, i5, i6, i7) -
                                   &left(0, 0, 0, 0, 0, 0, 0, 0);
                    if (j <= offset || left_alloc <= j) {
                      update |= 1;
                    }
                    offset = j;

                    if (&left(i0, i1, i2, i3, i4, i5, i6, i7) !=
                        &left_stride(i0, i1, i2, i3, i4, i5, i6, i7)) {
                      update |= 4;
                    }
                  }

    offset = -1;

    for (unsigned i0 = 0; i0 < unsigned(right.extent(0)); ++i0)
      for (unsigned i1 = 0; i1 < unsigned(right.extent(1)); ++i1)
        for (unsigned i2 = 0; i2 < unsigned(right.extent(2)); ++i2)
          for (unsigned i3 = 0; i3 < unsigned(right.extent(3)); ++i3)
            for (unsigned i4 = 0; i4 < unsigned(right.extent(4)); ++i4)
              for (unsigned i5 = 0; i5 < unsigned(right.extent(5)); ++i5)
                for (unsigned i6 = 0; i6 < unsigned(right.extent(6)); ++i6)
                  for (unsigned i7 = 0; i7 < unsigned(right.extent(7)); ++i7) {
                    const long j = &right(i0, i1, i2, i3, i4, i5, i6, i7) -
                                   &right(0, 0, 0, 0, 0, 0, 0, 0);
                    if (j <= offset || right_alloc <= j) {
                      update |= 2;
                    }
                    offset = j;

                    if (&right(i0, i1, i2, i3, i4, i5, i6, i7) !=
                        &right_stride(i0, i1, i2, i3, i4, i5, i6, i7)) {
                      update |= 8;
                    }
                  }
  }
};

template <class DataType, class DeviceType>
struct TestViewOperator_LeftAndRight<DataType, DeviceType, 7> {
  using execution_space = typename DeviceType::execution_space;
  using memory_space    = typename DeviceType::memory_space;
  using size_type       = typename execution_space::size_type;

  using value_type = int;

  KOKKOS_INLINE_FUNCTION
  static void join(value_type &update, const value_type &input) {
    update |= input;
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type &update) { update = 0; }

  using left_view = Kokkos::View<DataType, Kokkos::LayoutLeft, execution_space>;
  using right_view =
      Kokkos::View<DataType, Kokkos::LayoutRight, execution_space>;

  left_view left;
  right_view right;
  long left_alloc;
  long right_alloc;

  TestViewOperator_LeftAndRight()
      : left("left"),
        right("right"),
        left_alloc(allocation_count(left)),
        right_alloc(allocation_count(right)) {}

  void testit() {
    int error_flag = 0;

    Kokkos::parallel_reduce(1, *this, error_flag);

    ASSERT_EQ(error_flag, 0);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_type, value_type &update) const {
    long offset = -1;

    for (unsigned i6 = 0; i6 < unsigned(left.extent(6)); ++i6)
      for (unsigned i5 = 0; i5 < unsigned(left.extent(5)); ++i5)
        for (unsigned i4 = 0; i4 < unsigned(left.extent(4)); ++i4)
          for (unsigned i3 = 0; i3 < unsigned(left.extent(3)); ++i3)
            for (unsigned i2 = 0; i2 < unsigned(left.extent(2)); ++i2)
              for (unsigned i1 = 0; i1 < unsigned(left.extent(1)); ++i1)
                for (unsigned i0 = 0; i0 < unsigned(left.extent(0)); ++i0) {
                  const long j = &left(i0, i1, i2, i3, i4, i5, i6) -
                                 &left(0, 0, 0, 0, 0, 0, 0);
                  if (j <= offset || left_alloc <= j) {
                    update |= 1;
                  }
                  offset = j;
                }

    offset = -1;

    for (unsigned i0 = 0; i0 < unsigned(right.extent(0)); ++i0)
      for (unsigned i1 = 0; i1 < unsigned(right.extent(1)); ++i1)
        for (unsigned i2 = 0; i2 < unsigned(right.extent(2)); ++i2)
          for (unsigned i3 = 0; i3 < unsigned(right.extent(3)); ++i3)
            for (unsigned i4 = 0; i4 < unsigned(right.extent(4)); ++i4)
              for (unsigned i5 = 0; i5 < unsigned(right.extent(5)); ++i5)
                for (unsigned i6 = 0; i6 < unsigned(right.extent(6)); ++i6) {
                  const long j = &right(i0, i1, i2, i3, i4, i5, i6) -
                                 &right(0, 0, 0, 0, 0, 0, 0);
                  if (j <= offset || right_alloc <= j) {
                    update |= 2;
                  }
                  offset = j;
                }
  }
};

template <class DataType, class DeviceType>
struct TestViewOperator_LeftAndRight<DataType, DeviceType, 6> {
  using execution_space = typename DeviceType::execution_space;
  using memory_space    = typename DeviceType::memory_space;
  using size_type       = typename execution_space::size_type;

  using value_type = int;

  KOKKOS_INLINE_FUNCTION
  static void join(value_type &update, const value_type &input) {
    update |= input;
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type &update) { update = 0; }

  using left_view = Kokkos::View<DataType, Kokkos::LayoutLeft, execution_space>;
  using right_view =
      Kokkos::View<DataType, Kokkos::LayoutRight, execution_space>;

  left_view left;
  right_view right;
  long left_alloc;
  long right_alloc;

  TestViewOperator_LeftAndRight()
      : left("left"),
        right("right"),
        left_alloc(allocation_count(left)),
        right_alloc(allocation_count(right)) {}

  void testit() {
    int error_flag = 0;

    Kokkos::parallel_reduce(1, *this, error_flag);

    ASSERT_EQ(error_flag, 0);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_type, value_type &update) const {
    long offset = -1;

    for (unsigned i5 = 0; i5 < unsigned(left.extent(5)); ++i5)
      for (unsigned i4 = 0; i4 < unsigned(left.extent(4)); ++i4)
        for (unsigned i3 = 0; i3 < unsigned(left.extent(3)); ++i3)
          for (unsigned i2 = 0; i2 < unsigned(left.extent(2)); ++i2)
            for (unsigned i1 = 0; i1 < unsigned(left.extent(1)); ++i1)
              for (unsigned i0 = 0; i0 < unsigned(left.extent(0)); ++i0) {
                const long j =
                    &left(i0, i1, i2, i3, i4, i5) - &left(0, 0, 0, 0, 0, 0);
                if (j <= offset || left_alloc <= j) {
                  update |= 1;
                }
                offset = j;
              }

    offset = -1;

    for (unsigned i0 = 0; i0 < unsigned(right.extent(0)); ++i0)
      for (unsigned i1 = 0; i1 < unsigned(right.extent(1)); ++i1)
        for (unsigned i2 = 0; i2 < unsigned(right.extent(2)); ++i2)
          for (unsigned i3 = 0; i3 < unsigned(right.extent(3)); ++i3)
            for (unsigned i4 = 0; i4 < unsigned(right.extent(4)); ++i4)
              for (unsigned i5 = 0; i5 < unsigned(right.extent(5)); ++i5) {
                const long j =
                    &right(i0, i1, i2, i3, i4, i5) - &right(0, 0, 0, 0, 0, 0);
                if (j <= offset || right_alloc <= j) {
                  update |= 2;
                }
                offset = j;
              }
  }
};

template <class DataType, class DeviceType>
struct TestViewOperator_LeftAndRight<DataType, DeviceType, 5> {
  using execution_space = typename DeviceType::execution_space;
  using memory_space    = typename DeviceType::memory_space;
  using size_type       = typename execution_space::size_type;

  using value_type = int;

  KOKKOS_INLINE_FUNCTION
  static void join(value_type &update, const value_type &input) {
    update |= input;
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type &update) { update = 0; }

  using left_view = Kokkos::View<DataType, Kokkos::LayoutLeft, execution_space>;
  using right_view =
      Kokkos::View<DataType, Kokkos::LayoutRight, execution_space>;
  using stride_view =
      Kokkos::View<DataType, Kokkos::LayoutStride, execution_space>;

  left_view left;
  right_view right;
  stride_view left_stride;
  stride_view right_stride;
  long left_alloc;
  long right_alloc;

  TestViewOperator_LeftAndRight()
      : left("left"),
        right("right"),
        left_stride(left),
        right_stride(right),
        left_alloc(allocation_count(left)),
        right_alloc(allocation_count(right)) {}

  void testit() {
    int error_flag = 0;

    Kokkos::parallel_reduce(1, *this, error_flag);

    ASSERT_EQ(error_flag, 0);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_type, value_type &update) const {
    long offset = -1;

    for (unsigned i4 = 0; i4 < unsigned(left.extent(4)); ++i4)
      for (unsigned i3 = 0; i3 < unsigned(left.extent(3)); ++i3)
        for (unsigned i2 = 0; i2 < unsigned(left.extent(2)); ++i2)
          for (unsigned i1 = 0; i1 < unsigned(left.extent(1)); ++i1)
            for (unsigned i0 = 0; i0 < unsigned(left.extent(0)); ++i0) {
              const long j = &left(i0, i1, i2, i3, i4) - &left(0, 0, 0, 0, 0);
              if (j <= offset || left_alloc <= j) {
                update |= 1;
              }
              offset = j;

              if (&left(i0, i1, i2, i3, i4) !=
                  &left_stride(i0, i1, i2, i3, i4)) {
                update |= 4;
              }
            }

    offset = -1;

    for (unsigned i0 = 0; i0 < unsigned(right.extent(0)); ++i0)
      for (unsigned i1 = 0; i1 < unsigned(right.extent(1)); ++i1)
        for (unsigned i2 = 0; i2 < unsigned(right.extent(2)); ++i2)
          for (unsigned i3 = 0; i3 < unsigned(right.extent(3)); ++i3)
            for (unsigned i4 = 0; i4 < unsigned(right.extent(4)); ++i4) {
              const long j = &right(i0, i1, i2, i3, i4) - &right(0, 0, 0, 0, 0);
              if (j <= offset || right_alloc <= j) {
                update |= 2;
              }
              offset = j;

              if (&right(i0, i1, i2, i3, i4) !=
                  &right_stride(i0, i1, i2, i3, i4)) {
                update |= 8;
              }
            }
  }
};

template <class DataType, class DeviceType>
struct TestViewOperator_LeftAndRight<DataType, DeviceType, 4> {
  using execution_space = typename DeviceType::execution_space;
  using memory_space    = typename DeviceType::memory_space;
  using size_type       = typename execution_space::size_type;

  using value_type = int;

  KOKKOS_INLINE_FUNCTION
  static void join(value_type &update, const value_type &input) {
    update |= input;
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type &update) { update = 0; }

  using left_view = Kokkos::View<DataType, Kokkos::LayoutLeft, execution_space>;
  using right_view =
      Kokkos::View<DataType, Kokkos::LayoutRight, execution_space>;

  left_view left;
  right_view right;
  long left_alloc;
  long right_alloc;

  TestViewOperator_LeftAndRight()
      : left("left"),
        right("right"),
        left_alloc(allocation_count(left)),
        right_alloc(allocation_count(right)) {}

  void testit() {
    int error_flag = 0;

    Kokkos::parallel_reduce(1, *this, error_flag);

    ASSERT_EQ(error_flag, 0);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_type, value_type &update) const {
    long offset = -1;

    for (unsigned i3 = 0; i3 < unsigned(left.extent(3)); ++i3)
      for (unsigned i2 = 0; i2 < unsigned(left.extent(2)); ++i2)
        for (unsigned i1 = 0; i1 < unsigned(left.extent(1)); ++i1)
          for (unsigned i0 = 0; i0 < unsigned(left.extent(0)); ++i0) {
            const long j = &left(i0, i1, i2, i3) - &left(0, 0, 0, 0);
            if (j <= offset || left_alloc <= j) {
              update |= 1;
            }
            offset = j;
          }

    offset = -1;

    for (unsigned i0 = 0; i0 < unsigned(right.extent(0)); ++i0)
      for (unsigned i1 = 0; i1 < unsigned(right.extent(1)); ++i1)
        for (unsigned i2 = 0; i2 < unsigned(right.extent(2)); ++i2)
          for (unsigned i3 = 0; i3 < unsigned(right.extent(3)); ++i3) {
            const long j = &right(i0, i1, i2, i3) - &right(0, 0, 0, 0);
            if (j <= offset || right_alloc <= j) {
              update |= 2;
            }
            offset = j;
          }
  }
};

template <class DataType, class DeviceType>
struct TestViewOperator_LeftAndRight<DataType, DeviceType, 3> {
  using execution_space = typename DeviceType::execution_space;
  using memory_space    = typename DeviceType::memory_space;
  using size_type       = typename execution_space::size_type;

  using value_type = int;

  KOKKOS_INLINE_FUNCTION
  static void join(value_type &update, const value_type &input) {
    update |= input;
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type &update) { update = 0; }

  using left_view = Kokkos::View<DataType, Kokkos::LayoutLeft, execution_space>;
  using right_view =
      Kokkos::View<DataType, Kokkos::LayoutRight, execution_space>;
  using stride_view =
      Kokkos::View<DataType, Kokkos::LayoutStride, execution_space>;

  left_view left;
  right_view right;
  stride_view left_stride;
  stride_view right_stride;
  long left_alloc;
  long right_alloc;

  TestViewOperator_LeftAndRight()
      : left(std::string("left")),
        right(std::string("right")),
        left_stride(left),
        right_stride(right),
        left_alloc(allocation_count(left)),
        right_alloc(allocation_count(right)) {}

  void testit() {
    int error_flag = 0;

    Kokkos::parallel_reduce(1, *this, error_flag);

    ASSERT_EQ(error_flag, 0);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_type, value_type &update) const {
    long offset = -1;

    for (unsigned i2 = 0; i2 < unsigned(left.extent(2)); ++i2)
      for (unsigned i1 = 0; i1 < unsigned(left.extent(1)); ++i1)
        for (unsigned i0 = 0; i0 < unsigned(left.extent(0)); ++i0) {
          const long j = &left(i0, i1, i2) - &left(0, 0, 0);
          if (j <= offset || left_alloc <= j) {
            update |= 1;
          }
          offset = j;

          if (&left(i0, i1, i2) != &left_stride(i0, i1, i2)) {
            update |= 4;
          }
        }

    offset = -1;

    for (unsigned i0 = 0; i0 < unsigned(right.extent(0)); ++i0)
      for (unsigned i1 = 0; i1 < unsigned(right.extent(1)); ++i1)
        for (unsigned i2 = 0; i2 < unsigned(right.extent(2)); ++i2) {
          const long j = &right(i0, i1, i2) - &right(0, 0, 0);
          if (j <= offset || right_alloc <= j) {
            update |= 2;
          }
          offset = j;

          if (&right(i0, i1, i2) != &right_stride(i0, i1, i2)) {
            update |= 8;
          }
        }

    for (unsigned i0 = 0; i0 < unsigned(left.extent(0)); ++i0)
      for (unsigned i1 = 0; i1 < unsigned(left.extent(1)); ++i1)
        for (unsigned i2 = 0; i2 < unsigned(left.extent(2)); ++i2) {
          if (&left(i0, i1, i2) != &left.access(i0, i1, i2, 0, 0, 0, 0, 0)) {
            update |= 3;
          }
          if (&right(i0, i1, i2) != &right.access(i0, i1, i2, 0, 0, 0, 0, 0)) {
            update |= 3;
          }
        }
  }
};

template <class DataType, class DeviceType>
struct TestViewOperator_LeftAndRight<DataType, DeviceType, 2> {
  using execution_space = typename DeviceType::execution_space;
  using memory_space    = typename DeviceType::memory_space;
  using size_type       = typename execution_space::size_type;

  using value_type = int;

  KOKKOS_INLINE_FUNCTION
  static void join(value_type &update, const value_type &input) {
    update |= input;
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type &update) { update = 0; }

  using left_view = Kokkos::View<DataType, Kokkos::LayoutLeft, execution_space>;
  using right_view =
      Kokkos::View<DataType, Kokkos::LayoutRight, execution_space>;

  left_view left;
  right_view right;
  long left_alloc;
  long right_alloc;

  TestViewOperator_LeftAndRight()
      : left("left"),
        right("right"),
        left_alloc(allocation_count(left)),
        right_alloc(allocation_count(right)) {}

  void testit() {
    int error_flag = 0;

    Kokkos::parallel_reduce(1, *this, error_flag);

    ASSERT_EQ(error_flag, 0);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_type, value_type &update) const {
    long offset = -1;

    for (unsigned i1 = 0; i1 < unsigned(left.extent(1)); ++i1)
      for (unsigned i0 = 0; i0 < unsigned(left.extent(0)); ++i0) {
        const long j = &left(i0, i1) - &left(0, 0);
        if (j <= offset || left_alloc <= j) {
          update |= 1;
        }
        offset = j;
      }

    offset = -1;

    for (unsigned i0 = 0; i0 < unsigned(right.extent(0)); ++i0)
      for (unsigned i1 = 0; i1 < unsigned(right.extent(1)); ++i1) {
        const long j = &right(i0, i1) - &right(0, 0);
        if (j <= offset || right_alloc <= j) {
          update |= 2;
        }
        offset = j;
      }

    for (unsigned i0 = 0; i0 < unsigned(left.extent(0)); ++i0)
      for (unsigned i1 = 0; i1 < unsigned(left.extent(1)); ++i1) {
        if (&left(i0, i1) != &left.access(i0, i1, 0, 0, 0, 0, 0, 0)) {
          update |= 3;
        }
        if (&right(i0, i1) != &right.access(i0, i1, 0, 0, 0, 0, 0, 0)) {
          update |= 3;
        }
      }
  }
};

template <class DataType, class DeviceType>
struct TestViewOperator_LeftAndRight<DataType, DeviceType, 1> {
  using execution_space = typename DeviceType::execution_space;
  using memory_space    = typename DeviceType::memory_space;
  using size_type       = typename execution_space::size_type;

  using value_type = int;

  KOKKOS_INLINE_FUNCTION
  static void join(value_type &update, const value_type &input) {
    update |= input;
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type &update) { update = 0; }

  using left_view = Kokkos::View<DataType, Kokkos::LayoutLeft, execution_space>;
  using right_view =
      Kokkos::View<DataType, Kokkos::LayoutRight, execution_space>;
  using stride_view =
      Kokkos::View<DataType, Kokkos::LayoutStride, execution_space>;

  left_view left;
  right_view right;
  stride_view left_stride;
  stride_view right_stride;
  long left_alloc;
  long right_alloc;

  TestViewOperator_LeftAndRight()
      : left("left"),
        right("right"),
        left_stride(left),
        right_stride(right),
        left_alloc(allocation_count(left)),
        right_alloc(allocation_count(right)) {}

  void testit() {
    TestViewOperator_LeftAndRight driver;

    int error_flag = 0;

    Kokkos::parallel_reduce(1, *this, error_flag);

    ASSERT_EQ(error_flag, 0);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_type, value_type &update) const {
    for (unsigned i0 = 0; i0 < unsigned(left.extent(0)); ++i0) {
      if (&left(i0) != &left.access(i0, 0, 0, 0, 0, 0, 0, 0)) {
        update |= 3;
      }
      if (&right(i0) != &right.access(i0, 0, 0, 0, 0, 0, 0, 0)) {
        update |= 3;
      }
      if (&left(i0) != &left_stride(i0)) {
        update |= 4;
      }
      if (&right(i0) != &right_stride(i0)) {
        update |= 8;
      }
    }
  }
};

template <class Layout, class DeviceType>
struct TestViewMirror {
  template <class MemoryTraits>
  void static test_mirror() {
    Kokkos::View<double *, Layout, Kokkos::HostSpace> a_org("A", 1000);
    Kokkos::View<double *, Layout, Kokkos::HostSpace, MemoryTraits> a_h = a_org;
    auto a_h2 = Kokkos::create_mirror(Kokkos::HostSpace(), a_h);
    auto a_d  = Kokkos::create_mirror(DeviceType(), a_h);

    int equal_ptr_h_h2 = (a_h.data() == a_h2.data()) ? 1 : 0;
    int equal_ptr_h_d  = (a_h.data() == a_d.data()) ? 1 : 0;
    int equal_ptr_h2_d = (a_h2.data() == a_d.data()) ? 1 : 0;

    ASSERT_EQ(equal_ptr_h_h2, 0);
    ASSERT_EQ(equal_ptr_h_d, 0);
    ASSERT_EQ(equal_ptr_h2_d, 0);

    ASSERT_EQ(a_h.extent(0), a_h2.extent(0));
    ASSERT_EQ(a_h.extent(0), a_d.extent(0));
  }

  template <class MemoryTraits>
  void static test_mirror_view() {
    Kokkos::View<double *, Layout, Kokkos::HostSpace> a_org("A", 1000);
    Kokkos::View<double *, Layout, Kokkos::HostSpace, MemoryTraits> a_h = a_org;
    auto a_h2 = Kokkos::create_mirror_view(Kokkos::HostSpace(), a_h);
    auto a_d  = Kokkos::create_mirror_view(DeviceType(), a_h);

    int equal_ptr_h_h2 = a_h.data() == a_h2.data() ? 1 : 0;
    int equal_ptr_h_d  = a_h.data() == a_d.data() ? 1 : 0;
    int equal_ptr_h2_d = a_h2.data() == a_d.data() ? 1 : 0;

    int is_same_memspace =
        std::is_same<Kokkos::HostSpace,
                     typename DeviceType::memory_space>::value
            ? 1
            : 0;
    ASSERT_EQ(equal_ptr_h_h2, 1);
    ASSERT_EQ(equal_ptr_h_d, is_same_memspace);
    ASSERT_EQ(equal_ptr_h2_d, is_same_memspace);

    ASSERT_EQ(a_h.extent(0), a_h2.extent(0));
    ASSERT_EQ(a_h.extent(0), a_d.extent(0));
  }

  template <class MemoryTraits>
  void static test_mirror_copy() {
    Kokkos::View<double *, Layout, Kokkos::HostSpace> a_org("A", 10);
    a_org(5)                                                            = 42.0;
    Kokkos::View<double *, Layout, Kokkos::HostSpace, MemoryTraits> a_h = a_org;
    auto a_h2 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), a_h);
    auto a_d  = Kokkos::create_mirror_view_and_copy(DeviceType(), a_h);
    auto a_h3 = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), a_d);

    int equal_ptr_h_h2 = a_h.data() == a_h2.data() ? 1 : 0;
    int equal_ptr_h_d  = a_h.data() == a_d.data() ? 1 : 0;
    int equal_ptr_h2_d = a_h2.data() == a_d.data() ? 1 : 0;
    int equal_ptr_h3_d = a_h3.data() == a_d.data() ? 1 : 0;

    int is_same_memspace =
        std::is_same<Kokkos::HostSpace,
                     typename DeviceType::memory_space>::value
            ? 1
            : 0;
    ASSERT_EQ(equal_ptr_h_h2, 1);
    ASSERT_EQ(equal_ptr_h_d, is_same_memspace);
    ASSERT_EQ(equal_ptr_h2_d, is_same_memspace);
    ASSERT_EQ(equal_ptr_h3_d, is_same_memspace);

    ASSERT_EQ(a_h.extent(0), a_h3.extent(0));
    ASSERT_EQ(a_h.extent(0), a_h2.extent(0));
    ASSERT_EQ(a_h.extent(0), a_d.extent(0));
    ASSERT_EQ(a_org(5), a_h3(5));
  }

  template <typename View>
  static typename View::const_type view_const_cast(View const &v) {
    return v;
  }

  static void test_allocated() {
    using ExecutionSpace = typename DeviceType::execution_space;
    using dynamic_view   = Kokkos::View<int *, ExecutionSpace>;
    using static_view    = Kokkos::View<int[5], ExecutionSpace>;
    using unmanaged_view =
        Kokkos::View<int *, ExecutionSpace,
                     Kokkos::MemoryTraits<Kokkos::Unmanaged> >;
    int const N = 100;

    dynamic_view d1;
    static_view s1;
    unmanaged_view u1;
    ASSERT_FALSE(d1.is_allocated());
    ASSERT_FALSE(s1.is_allocated());
    ASSERT_FALSE(u1.is_allocated());

    d1 = dynamic_view("d1", N);
    dynamic_view d2(d1);
    dynamic_view d3("d3", N);
    ASSERT_TRUE(d1.is_allocated());
    ASSERT_TRUE(d2.is_allocated());
    ASSERT_TRUE(d3.is_allocated());

    s1 = static_view("s1");
    static_view s2(s1);
    static_view s3("s3");
    ASSERT_TRUE(s1.is_allocated());
    ASSERT_TRUE(s2.is_allocated());
    ASSERT_TRUE(s3.is_allocated());

    u1 = unmanaged_view(d1.data(), N);
    unmanaged_view u2(u1);
    unmanaged_view u3(d1.data(), N);
    ASSERT_TRUE(u1.is_allocated());
    ASSERT_TRUE(u2.is_allocated());
    ASSERT_TRUE(u3.is_allocated());
  }

  static void test_mirror_copy_const_data_type() {
    using ExecutionSpace = typename DeviceType::execution_space;
    int const N          = 100;
    Kokkos::View<int *, ExecutionSpace> v("v", N);
    Kokkos::deep_copy(v, 255);
    auto v_m1 = Kokkos::create_mirror_view_and_copy(
        Kokkos::DefaultHostExecutionSpace(), view_const_cast(v));
    auto v_m2 = Kokkos::create_mirror_view_and_copy(ExecutionSpace(),
                                                    view_const_cast(v));
  }

  template <class MemoryTraits, class Space>
  struct CopyUnInit {
    using mirror_view_type = typename Kokkos::Impl::MirrorViewType<
        Space, double *, Layout, Kokkos::HostSpace, MemoryTraits>::view_type;

    mirror_view_type a_d;

    KOKKOS_INLINE_FUNCTION
    CopyUnInit(mirror_view_type &a_d_) : a_d(a_d_) {}

    KOKKOS_INLINE_FUNCTION
    void operator()(const typename Space::size_type i) const {
      a_d(i) = (double)(10 - i);
    }
  };

  template <class MemoryTraits>
  void static test_mirror_no_initialize() {
    Kokkos::View<double *, Layout, Kokkos::HostSpace> a_org("A", 10);
    Kokkos::View<double *, Layout, Kokkos::HostSpace, MemoryTraits> a_h = a_org;

    for (int i = 0; i < 10; i++) {
      a_h(i) = (double)i;
    }
    auto a_d = Kokkos::create_mirror_view(Kokkos::WithoutInitializing,
                                          DeviceType(), a_h);

    int equal_ptr_h_d = (a_h.data() == a_d.data()) ? 1 : 0;
    constexpr int is_same_memspace =
        std::is_same<Kokkos::HostSpace,
                     typename DeviceType::memory_space>::value
            ? 1
            : 0;

    ASSERT_EQ(equal_ptr_h_d, is_same_memspace);

    Kokkos::parallel_for(
        Kokkos::RangePolicy<typename DeviceType::execution_space>(0, int(10)),
        CopyUnInit<MemoryTraits, DeviceType>(a_d));

    Kokkos::deep_copy(a_h, a_d);

    for (int i = 0; i < 10; i++) {
      ASSERT_EQ(a_h(i), (double)(10 - i));
    }
  }

  void static testit() {
    test_mirror<Kokkos::MemoryTraits<0> >();
    test_mirror<Kokkos::MemoryTraits<Kokkos::Unmanaged> >();
    test_mirror_view<Kokkos::MemoryTraits<0> >();
    test_mirror_view<Kokkos::MemoryTraits<Kokkos::Unmanaged> >();
    test_mirror_copy<Kokkos::MemoryTraits<0> >();
    test_mirror_copy<Kokkos::MemoryTraits<Kokkos::Unmanaged> >();
    test_mirror_copy_const_data_type();
    test_allocated();
    test_mirror_no_initialize<Kokkos::MemoryTraits<0> >();
    test_mirror_no_initialize<Kokkos::MemoryTraits<Kokkos::Unmanaged> >();
  }
};

/*--------------------------------------------------------------------------*/

template <typename T, class DeviceType>
class TestViewAPI {
 public:
  using device = DeviceType;

  enum { N0 = 1000, N1 = 3, N2 = 5, N3 = 7 };

  using dView0       = Kokkos::View<T, device>;
  using dView1       = Kokkos::View<T *, device>;
  using dView2       = Kokkos::View<T * [N1], device>;
  using dView3       = Kokkos::View<T * [N1][N2], device>;
  using dView4       = Kokkos::View<T * [N1][N2][N3], device>;
  using const_dView4 = Kokkos::View<const T * [N1][N2][N3], device>;
  using dView4_unmanaged =
      Kokkos::View<T ****, device, Kokkos::MemoryUnmanaged>;
  using host = typename dView0::host_mirror_space;

  static void run_test_view_operator_a() {
    {
      TestViewOperator<T, device> f;
      Kokkos::parallel_for(int(N0), f);
      Kokkos::fence();
    }
#ifndef KOKKOS_ENABLE_OPENMPTARGET
    TestViewOperator_LeftAndRight<int[2][3][4][2][3][4], device> f6;
    f6.testit();
    TestViewOperator_LeftAndRight<int[2][3][4][2][3], device> f5;
    f5.testit();
    TestViewOperator_LeftAndRight<int[2][3][4][2], device> f4;
    f4.testit();
    TestViewOperator_LeftAndRight<int[2][3][4], device> f3;
    f3.testit();
    TestViewOperator_LeftAndRight<int[2][3], device> f2;
    f2.testit();
    TestViewOperator_LeftAndRight<int[2], device> f1;
    f1.testit();
#endif
  }

  static void run_test_view_operator_b() {
#ifndef KOKKOS_ENABLE_OPENMPTARGET
    TestViewOperator_LeftAndRight<int[2][3][4][2][3][4][2], device> f7;
    f7.testit();
#endif
  }

  static void run_test_view_operator_c() {
#ifndef KOKKOS_ENABLE_OPENMPTARGET
    TestViewOperator_LeftAndRight<int[2][3][4][2][3][4][2][3], device> f8;
    f8.testit();
#endif
  }

  static void run_test_mirror() {
    using view_type   = Kokkos::View<int, host>;
    using mirror_type = typename view_type::HostMirror;

    static_assert(std::is_same<typename view_type::memory_space,
                               typename mirror_type::memory_space>::value,
                  "");

    view_type a("a");
    mirror_type am = Kokkos::create_mirror_view(a);
    mirror_type ax = Kokkos::create_mirror(a);
    ASSERT_EQ(&a(), &am());

    TestViewMirror<Kokkos::LayoutLeft, device>::testit();
    TestViewMirror<Kokkos::LayoutRight, device>::testit();
  }

  static void run_test_scalar() {
    using hView0 = typename dView0::HostMirror;

    dView0 dx, dy;
    hView0 hx, hy;

    dx = dView0("dx");
    dy = dView0("dy");

    hx = Kokkos::create_mirror(dx);
    hy = Kokkos::create_mirror(dy);

    hx() = 1;

    Kokkos::deep_copy(dx, hx);
    Kokkos::deep_copy(dy, dx);
    Kokkos::deep_copy(hy, dy);
#ifndef KOKKOS_ENABLE_OPENMPTARGET
    ASSERT_EQ(hx(), hy());
#endif
  }

  static void run_test_contruction_from_layout() {
    using hView0 = typename dView0::HostMirror;
    using hView1 = typename dView1::HostMirror;
    using hView2 = typename dView2::HostMirror;
    using hView3 = typename dView3::HostMirror;
    using hView4 = typename dView4::HostMirror;

    hView0 hv_0("dView0::HostMirror");
    hView1 hv_1("dView1::HostMirror", N0);
    hView2 hv_2("dView2::HostMirror", N0);
    hView3 hv_3("dView3::HostMirror", N0);
    hView4 hv_4("dView4::HostMirror", N0);

    dView0 dv_0_1(nullptr, 0);
    dView0 dv_0_2(hv_0.label(), hv_0.layout());

    dView1 dv_1_1(nullptr, 0);
    dView1 dv_1_2(hv_1.label(), hv_1.layout());

    dView2 dv_2_1(nullptr, 0);
    dView2 dv_2_2(hv_2.label(), hv_2.layout());

    dView3 dv_3_1(nullptr, 0);
    dView3 dv_3_2(hv_3.label(), hv_3.layout());

    dView4 dv_4_1(nullptr, 0);
    dView4 dv_4_2(hv_4.label(), hv_4.layout());
  }

  static void run_test_contruction_from_layout_2() {
    using dView3_0 = Kokkos::View<T ***, device>;
    using dView3_1 = Kokkos::View<T * * [N1], device>;
    using dView3_2 = Kokkos::View<T * [N1][N2], device>;
    using dView3_3 = Kokkos::View<T[N0][N1][N2], device>;

    dView3_0 v_0("v_0", N0, N1, N2);
    dView3_1 v_1("v_1", N0, N1);
    dView3_2 v_2("v_2", N0);
    dView3_3 v_3("v_2");

    dView3_1 v_1_a("v_1", N0, N1, N2);
    dView3_2 v_2_a("v_2", N0, N1, N2);
    dView3_3 v_3_a("v_2", N0, N1, N2);

    {
      dView3_0 dv_1(v_0.label(), v_0.layout());
      dView3_0 dv_2(v_1.label(), v_1.layout());
      dView3_0 dv_3(v_2.label(), v_2.layout());
      dView3_0 dv_4(v_3.label(), v_3.layout());
      dView3_0 dv_5(v_1_a.label(), v_1_a.layout());
      dView3_0 dv_6(v_2_a.label(), v_2_a.layout());
      dView3_0 dv_7(v_3_a.label(), v_3_a.layout());
    }

    {
      dView3_1 dv_1(v_0.label(), v_0.layout());
      dView3_1 dv_2(v_1.label(), v_1.layout());
      dView3_1 dv_3(v_2.label(), v_2.layout());
      dView3_1 dv_4(v_3.label(), v_3.layout());
      dView3_1 dv_5(v_1_a.label(), v_1_a.layout());
      dView3_1 dv_6(v_2_a.label(), v_2_a.layout());
      dView3_1 dv_7(v_3_a.label(), v_3_a.layout());
    }

    {
      dView3_2 dv_1(v_0.label(), v_0.layout());
      dView3_2 dv_2(v_1.label(), v_1.layout());
      dView3_2 dv_3(v_2.label(), v_2.layout());
      dView3_2 dv_4(v_3.label(), v_3.layout());
      dView3_2 dv_5(v_1_a.label(), v_1_a.layout());
      dView3_2 dv_6(v_2_a.label(), v_2_a.layout());
      dView3_2 dv_7(v_3_a.label(), v_3_a.layout());
    }

    {
      dView3_3 dv_1(v_0.label(), v_0.layout());
      dView3_3 dv_2(v_1.label(), v_1.layout());
      dView3_3 dv_3(v_2.label(), v_2.layout());
      dView3_3 dv_4(v_3.label(), v_3.layout());
      dView3_3 dv_5(v_1_a.label(), v_1_a.layout());
      dView3_3 dv_6(v_2_a.label(), v_2_a.layout());
      dView3_3 dv_7(v_3_a.label(), v_3_a.layout());
    }
  }

  static void run_test() {
    // mfh 14 Feb 2014: This test doesn't actually create instances of
    // these types.  In order to avoid "unused type alias"
    // warnings, we declare empty instances of these types, with the
    // usual "(void)" marker to avoid compiler warnings for unused
    // variables.

    using hView0 = typename dView0::HostMirror;
    using hView1 = typename dView1::HostMirror;
    using hView2 = typename dView2::HostMirror;
    using hView3 = typename dView3::HostMirror;
    using hView4 = typename dView4::HostMirror;

    {
      hView0 thing;
      (void)thing;
    }
    {
      hView1 thing;
      (void)thing;
    }
    {
      hView2 thing;
      (void)thing;
    }
    {
      hView3 thing;
      (void)thing;
    }
    {
      hView4 thing;
      (void)thing;
    }

    dView4 dx, dy, dz;
    hView4 hx, hy, hz;

    ASSERT_EQ(dx.data(), nullptr);
    ASSERT_EQ(dy.data(), nullptr);
    ASSERT_EQ(dz.data(), nullptr);
    ASSERT_EQ(hx.data(), nullptr);
    ASSERT_EQ(hy.data(), nullptr);
    ASSERT_EQ(hz.data(), nullptr);
    ASSERT_EQ(dx.extent(0), 0u);
    ASSERT_EQ(dy.extent(0), 0u);
    ASSERT_EQ(dz.extent(0), 0u);
    ASSERT_EQ(hx.extent(0), 0u);
    ASSERT_EQ(hy.extent(0), 0u);
    ASSERT_EQ(hz.extent(0), 0u);
    ASSERT_EQ(dx.extent(1), unsigned(N1));
    ASSERT_EQ(dy.extent(1), unsigned(N1));
    ASSERT_EQ(dz.extent(1), unsigned(N1));
    ASSERT_EQ(hx.extent(1), unsigned(N1));
    ASSERT_EQ(hy.extent(1), unsigned(N1));
    ASSERT_EQ(hz.extent(1), unsigned(N1));

    dx = dView4("dx", N0);
    dy = dView4("dy", N0);

    ASSERT_EQ(dx.use_count(), 1);

    dView4_unmanaged unmanaged_dx = dx;
    ASSERT_EQ(dx.use_count(), 1);

    // Test self assignment
#if defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wself-assign-overloaded"
#endif
    dx = dx;  // copy-assignment operator
#if defined(__clang__)
#pragma GCC diagnostic pop
#endif
    ASSERT_EQ(dx.use_count(), 1);
    dx = reinterpret_cast<typename dView4::uniform_type &>(
        dx);  // conversion assignment operator
    ASSERT_EQ(dx.use_count(), 1);

    dView4_unmanaged unmanaged_from_ptr_dx = dView4_unmanaged(
        dx.data(), dx.extent(0), dx.extent(1), dx.extent(2), dx.extent(3));

    {
      // Destruction of this view should be harmless.

      const_dView4 unmanaged_from_ptr_const_dx(dx.data(), dx.extent(0));
    }

    const_dView4 const_dx = dx;
    ASSERT_EQ(dx.use_count(), 2);

    {
      const_dView4 const_dx2;
      const_dx2 = const_dx;
      ASSERT_EQ(dx.use_count(), 3);

      const_dx2 = dy;
      ASSERT_EQ(dx.use_count(), 2);

      const_dView4 const_dx3(dx);
      ASSERT_EQ(dx.use_count(), 3);

      dView4_unmanaged dx4_unmanaged(dx);
      ASSERT_EQ(dx.use_count(), 3);
    }

    ASSERT_EQ(dx.use_count(), 2);

    ASSERT_NE(dx.data(), nullptr);
    ASSERT_NE(const_dx.data(), nullptr);
    ASSERT_NE(unmanaged_dx.data(), nullptr);
    ASSERT_NE(unmanaged_from_ptr_dx.data(), nullptr);
    ASSERT_NE(dy.data(), nullptr);
    ASSERT_NE(dx, dy);

    ASSERT_EQ(dx.extent(0), unsigned(N0));
    ASSERT_EQ(dx.extent(1), unsigned(N1));
    ASSERT_EQ(dx.extent(2), unsigned(N2));
    ASSERT_EQ(dx.extent(3), unsigned(N3));

    ASSERT_EQ(dy.extent(0), unsigned(N0));
    ASSERT_EQ(dy.extent(1), unsigned(N1));
    ASSERT_EQ(dy.extent(2), unsigned(N2));
    ASSERT_EQ(dy.extent(3), unsigned(N3));

    ASSERT_EQ(unmanaged_from_ptr_dx.span(),
              unsigned(N0) * unsigned(N1) * unsigned(N2) * unsigned(N3));
#ifdef KOKKOS_ENABLE_OPENMPTARGET
    return;
#endif
    hx = Kokkos::create_mirror(dx);
    hy = Kokkos::create_mirror(dy);

    // T v1 = hx();       // Generates compile error as intended.
    // T v2 = hx( 0, 0 ); // Generates compile error as intended.
    // hx( 0, 0 ) = v2;   // Generates compile error as intended.

    // Testing with asynchronous deep copy with respect to device
    {
      size_t count = 0;

      for (size_t ip = 0; ip < N0; ++ip)
        for (size_t i1 = 0; i1 < hx.extent(1); ++i1)
          for (size_t i2 = 0; i2 < hx.extent(2); ++i2)
            for (size_t i3 = 0; i3 < hx.extent(3); ++i3) {
              hx(ip, i1, i2, i3) = ++count;
            }

      Kokkos::deep_copy(typename hView4::execution_space(), dx, hx);
      Kokkos::deep_copy(typename hView4::execution_space(), dy, dx);
      Kokkos::deep_copy(typename hView4::execution_space(), hy, dy);
      typename hView4::execution_space().fence();

      for (size_t ip = 0; ip < N0; ++ip)
        for (size_t i1 = 0; i1 < N1; ++i1)
          for (size_t i2 = 0; i2 < N2; ++i2)
            for (size_t i3 = 0; i3 < N3; ++i3) {
              ASSERT_EQ(hx(ip, i1, i2, i3), hy(ip, i1, i2, i3));
            }

      Kokkos::deep_copy(typename hView4::execution_space(), dx, T(0));
      Kokkos::deep_copy(typename hView4::execution_space(), hx, dx);
      typename hView4::execution_space().fence();

      for (size_t ip = 0; ip < N0; ++ip)
        for (size_t i1 = 0; i1 < N1; ++i1)
          for (size_t i2 = 0; i2 < N2; ++i2)
            for (size_t i3 = 0; i3 < N3; ++i3) {
              ASSERT_EQ(hx(ip, i1, i2, i3), T(0));
            }
    }

    // Testing with asynchronous deep copy with respect to host.
    {
      size_t count = 0;

      for (size_t ip = 0; ip < N0; ++ip)
        for (size_t i1 = 0; i1 < hx.extent(1); ++i1)
          for (size_t i2 = 0; i2 < hx.extent(2); ++i2)
            for (size_t i3 = 0; i3 < hx.extent(3); ++i3) {
              hx(ip, i1, i2, i3) = ++count;
            }

      Kokkos::deep_copy(typename dView4::execution_space(), dx, hx);
      Kokkos::deep_copy(typename dView4::execution_space(), dy, dx);
      Kokkos::deep_copy(typename dView4::execution_space(), hy, dy);
      typename dView4::execution_space().fence();

      for (size_t ip = 0; ip < N0; ++ip)
        for (size_t i1 = 0; i1 < N1; ++i1)
          for (size_t i2 = 0; i2 < N2; ++i2)
            for (size_t i3 = 0; i3 < N3; ++i3) {
              ASSERT_EQ(hx(ip, i1, i2, i3), hy(ip, i1, i2, i3));
            }

      Kokkos::deep_copy(typename dView4::execution_space(), dx, T(0));
      Kokkos::deep_copy(typename dView4::execution_space(), hx, dx);
      typename dView4::execution_space().fence();

      for (size_t ip = 0; ip < N0; ++ip)
        for (size_t i1 = 0; i1 < N1; ++i1)
          for (size_t i2 = 0; i2 < N2; ++i2)
            for (size_t i3 = 0; i3 < N3; ++i3) {
              ASSERT_EQ(hx(ip, i1, i2, i3), T(0));
            }
    }

    // Testing with synchronous deep copy.
    {
      size_t count = 0;

      for (size_t ip = 0; ip < N0; ++ip)
        for (size_t i1 = 0; i1 < hx.extent(1); ++i1)
          for (size_t i2 = 0; i2 < hx.extent(2); ++i2)
            for (size_t i3 = 0; i3 < hx.extent(3); ++i3) {
              hx(ip, i1, i2, i3) = ++count;
            }

      Kokkos::deep_copy(dx, hx);
      Kokkos::deep_copy(dy, dx);
      Kokkos::deep_copy(hy, dy);

      for (size_t ip = 0; ip < N0; ++ip)
        for (size_t i1 = 0; i1 < N1; ++i1)
          for (size_t i2 = 0; i2 < N2; ++i2)
            for (size_t i3 = 0; i3 < N3; ++i3) {
              ASSERT_EQ(hx(ip, i1, i2, i3), hy(ip, i1, i2, i3));
            }

      Kokkos::deep_copy(dx, T(0));
      Kokkos::deep_copy(hx, dx);

      for (size_t ip = 0; ip < N0; ++ip)
        for (size_t i1 = 0; i1 < N1; ++i1)
          for (size_t i2 = 0; i2 < N2; ++i2)
            for (size_t i3 = 0; i3 < N3; ++i3) {
              ASSERT_EQ(hx(ip, i1, i2, i3), T(0));
            }
    }

    dz = dx;
    ASSERT_EQ(dx, dz);
    ASSERT_NE(dy, dz);

    dz = dy;
    ASSERT_EQ(dy, dz);
    ASSERT_NE(dx, dz);

    dx = dView4();
    ASSERT_EQ(dx.data(), nullptr);
    ASSERT_NE(dy.data(), nullptr);
    ASSERT_NE(dz.data(), nullptr);

    dy = dView4();
    ASSERT_EQ(dx.data(), nullptr);
    ASSERT_EQ(dy.data(), nullptr);
    ASSERT_NE(dz.data(), nullptr);

    dz = dView4();
    ASSERT_EQ(dx.data(), nullptr);
    ASSERT_EQ(dy.data(), nullptr);
    ASSERT_EQ(dz.data(), nullptr);
  }

  static void run_test_deep_copy_empty() {
    // Check Deep Copy of LayoutLeft to LayoutRight
    {
      Kokkos::View<double *, Kokkos::LayoutLeft> dll("dll", 10);
      Kokkos::View<double *, Kokkos::LayoutRight, Kokkos::HostSpace> hlr("hlr",
                                                                         10);
      Kokkos::deep_copy(dll, hlr);
      Kokkos::deep_copy(hlr, dll);
    }

    // Check Deep Copy of two empty 1D views
    {
      Kokkos::View<double *> d;
      Kokkos::View<double *, Kokkos::HostSpace> h;
      Kokkos::deep_copy(d, h);
      Kokkos::deep_copy(h, d);
    }

    // Check Deep Copy of two empty 2D views
    {
      Kokkos::View<double * [3], Kokkos::LayoutRight> d;
      Kokkos::View<double * [3], Kokkos::LayoutRight, Kokkos::HostSpace> h;
      Kokkos::deep_copy(d, h);
      Kokkos::deep_copy(h, d);
    }
  }

  using DataType = T[2];

  static void check_auto_conversion_to_const(
      const Kokkos::View<const DataType, device> &arg_const,
      const Kokkos::View<DataType, device> &arg) {
    ASSERT_EQ(arg_const, arg);
  }

  static void run_test_const() {
    using typeX       = Kokkos::View<DataType, device>;
    using const_typeX = Kokkos::View<const DataType, device>;
    using const_typeR =
        Kokkos::View<const DataType, device, Kokkos::MemoryRandomAccess>;

    typeX x("X");
    const_typeX xc = x;
    const_typeR xr = x;

    ASSERT_EQ(xc, x);
    ASSERT_EQ(x, xc);

    // For CUDA the constant random access View does not return
    // an lvalue reference due to retrieving through texture cache
    // therefore not allowed to query the underlying pointer.
#if defined(KOKKOS_ENABLE_CUDA)
    if (!std::is_same<typename device::execution_space, Kokkos::Cuda>::value)
#endif
    {
      ASSERT_EQ(x.data(), xr.data());
    }

    // typeX xf = xc; // Setting non-const from const must not compile.

    check_auto_conversion_to_const(x, x);
  }

  static void run_test_subview() {
    using sView = Kokkos::View<const T, device>;

    dView0 d0("d0");
    dView1 d1("d1", N0);
    dView2 d2("d2", N0);
    dView3 d3("d3", N0);
    dView4 d4("d4", N0);

    sView s0 = d0;
    sView s1 = Kokkos::subview(d1, 1);
    sView s2 = Kokkos::subview(d2, 1, 1);
    sView s3 = Kokkos::subview(d3, 1, 1, 1);
    sView s4 = Kokkos::subview(d4, 1, 1, 1, 1);
  }

  static void run_test_subview_strided() {
    using view_left_4  = Kokkos::View<int ****, Kokkos::LayoutLeft, host>;
    using view_right_4 = Kokkos::View<int ****, Kokkos::LayoutRight, host>;
    using view_left_2  = Kokkos::View<int **, Kokkos::LayoutLeft, host>;
    using view_right_2 = Kokkos::View<int **, Kokkos::LayoutRight, host>;

    using view_stride_1 = Kokkos::View<int *, Kokkos::LayoutStride, host>;
    using view_stride_2 = Kokkos::View<int **, Kokkos::LayoutStride, host>;

    view_left_2 xl2("xl2", 100, 200);
    view_right_2 xr2("xr2", 100, 200);
    view_stride_1 yl1 = Kokkos::subview(xl2, 0, Kokkos::ALL());
    view_stride_1 yl2 = Kokkos::subview(xl2, 1, Kokkos::ALL());
    view_stride_1 yr1 = Kokkos::subview(xr2, 0, Kokkos::ALL());
    view_stride_1 yr2 = Kokkos::subview(xr2, 1, Kokkos::ALL());

    ASSERT_EQ(yl1.extent(0), xl2.extent(1));
    ASSERT_EQ(yl2.extent(0), xl2.extent(1));
    ASSERT_EQ(yr1.extent(0), xr2.extent(1));
    ASSERT_EQ(yr2.extent(0), xr2.extent(1));

    ASSERT_EQ(&yl1(0) - &xl2(0, 0), 0);
    ASSERT_EQ(&yl2(0) - &xl2(1, 0), 0);
    ASSERT_EQ(&yr1(0) - &xr2(0, 0), 0);
    ASSERT_EQ(&yr2(0) - &xr2(1, 0), 0);

    view_left_4 xl4("xl4", 10, 20, 30, 40);
    view_right_4 xr4("xr4", 10, 20, 30, 40);

    view_stride_2 yl4 =
        Kokkos::subview(xl4, 1, Kokkos::ALL(), 2, Kokkos::ALL());
    view_stride_2 yr4 =
        Kokkos::subview(xr4, 1, Kokkos::ALL(), 2, Kokkos::ALL());

    ASSERT_EQ(yl4.extent(0), xl4.extent(1));
    ASSERT_EQ(yl4.extent(1), xl4.extent(3));
    ASSERT_EQ(yr4.extent(0), xr4.extent(1));
    ASSERT_EQ(yr4.extent(1), xr4.extent(3));

    ASSERT_EQ(&yl4(4, 4) - &xl4(1, 4, 2, 4), 0);
    ASSERT_EQ(&yr4(4, 4) - &xr4(1, 4, 2, 4), 0);
  }

  static void run_test_vector() {
    static const unsigned Length = 1000, Count = 8;

    using vector_type      = Kokkos::View<T *, Kokkos::LayoutLeft, host>;
    using multivector_type = Kokkos::View<T **, Kokkos::LayoutLeft, host>;

    using vector_right_type = Kokkos::View<T *, Kokkos::LayoutRight, host>;
    using multivector_right_type =
        Kokkos::View<T **, Kokkos::LayoutRight, host>;

    using const_vector_right_type =
        Kokkos::View<const T *, Kokkos::LayoutRight, host>;
    using const_vector_type = Kokkos::View<const T *, Kokkos::LayoutLeft, host>;
    using const_multivector_type =
        Kokkos::View<const T **, Kokkos::LayoutLeft, host>;

    multivector_type mv = multivector_type("mv", Length, Count);
    multivector_right_type mv_right =
        multivector_right_type("mv", Length, Count);

    vector_type v1 = Kokkos::subview(mv, Kokkos::ALL(), 0);
    vector_type v2 = Kokkos::subview(mv, Kokkos::ALL(), 1);
    vector_type v3 = Kokkos::subview(mv, Kokkos::ALL(), 2);

    vector_type rv1 = Kokkos::subview(mv_right, 0, Kokkos::ALL());
    vector_type rv2 = Kokkos::subview(mv_right, 1, Kokkos::ALL());
    vector_type rv3 = Kokkos::subview(mv_right, 2, Kokkos::ALL());

    multivector_type mv1 =
        Kokkos::subview(mv, std::make_pair(1, 998), std::make_pair(2, 5));

    multivector_right_type mvr1 =
        Kokkos::subview(mv_right, std::make_pair(1, 998), std::make_pair(2, 5));

    const_vector_type cv1 = Kokkos::subview(mv, Kokkos::ALL(), 0);
    const_vector_type cv2 = Kokkos::subview(mv, Kokkos::ALL(), 1);
    const_vector_type cv3 = Kokkos::subview(mv, Kokkos::ALL(), 2);

    vector_right_type vr1 = Kokkos::subview(mv, Kokkos::ALL(), 0);
    vector_right_type vr2 = Kokkos::subview(mv, Kokkos::ALL(), 1);
    vector_right_type vr3 = Kokkos::subview(mv, Kokkos::ALL(), 2);

    const_vector_right_type cvr1 = Kokkos::subview(mv, Kokkos::ALL(), 0);
    const_vector_right_type cvr2 = Kokkos::subview(mv, Kokkos::ALL(), 1);
    const_vector_right_type cvr3 = Kokkos::subview(mv, Kokkos::ALL(), 2);

    ASSERT_EQ(&v1[0], &v1(0));
    ASSERT_EQ(&v1[0], &mv(0, 0));
    ASSERT_EQ(&v2[0], &mv(0, 1));
    ASSERT_EQ(&v3[0], &mv(0, 2));

    ASSERT_EQ(&cv1[0], &mv(0, 0));
    ASSERT_EQ(&cv2[0], &mv(0, 1));
    ASSERT_EQ(&cv3[0], &mv(0, 2));

    ASSERT_EQ(&vr1[0], &mv(0, 0));
    ASSERT_EQ(&vr2[0], &mv(0, 1));
    ASSERT_EQ(&vr3[0], &mv(0, 2));

    ASSERT_EQ(&cvr1[0], &mv(0, 0));
    ASSERT_EQ(&cvr2[0], &mv(0, 1));
    ASSERT_EQ(&cvr3[0], &mv(0, 2));

    ASSERT_EQ(&mv1(0, 0), &mv(1, 2));
    ASSERT_EQ(&mv1(1, 1), &mv(2, 3));
    ASSERT_EQ(&mv1(3, 2), &mv(4, 4));
    ASSERT_EQ(&mvr1(0, 0), &mv_right(1, 2));
    ASSERT_EQ(&mvr1(1, 1), &mv_right(2, 3));
    ASSERT_EQ(&mvr1(3, 2), &mv_right(4, 4));

    const_vector_type c_cv1(v1);
    typename vector_type::const_type c_cv2(v2);
    typename const_vector_type::const_type c_ccv2(v2);

    const_multivector_type cmv(mv);
    typename multivector_type::const_type cmvX(cmv);
    typename const_multivector_type::const_type ccmvX(cmv);
  }

  static void run_test_error() {
#ifdef KOKKOS_ENABLE_OPENMPTARGET
    if (std::is_same<typename dView1::memory_space,
                     Kokkos::Experimental::OpenMPTargetSpace>::value)
      return;
#endif
// FIXME_MSVC_WITH_CUDA
// This test doesn't behave as expected on Windows with CUDA
#if defined(_WIN32) && defined(KOKKOS_ENABLE_CUDA)
    if (std::is_same<typename dView1::memory_space,
                     Kokkos::CudaUVMSpace>::value)
      return;
#endif
    auto alloc_size = std::numeric_limits<size_t>::max() - 42;
    try {
      auto should_always_fail = dView1("hello_world_failure", alloc_size);
    } catch (std::runtime_error const &error) {
      // TODO once we remove the conversion to std::runtime_error, catch the
      //      appropriate Kokkos error here
      std::string msg = error.what();
      ASSERT_PRED_FORMAT2(::testing::IsSubstring, "hello_world_failure", msg);
      ASSERT_PRED_FORMAT2(::testing::IsSubstring,
                          typename device::memory_space{}.name(), msg);
      // Can't figure out how to make assertions either/or, so we'll just use
      // an if statement here for now.  Test failure message will be a bit
      // misleading, but developers should figure out what's going on pretty
      // quickly.
      if (msg.find("is not a valid size") != std::string::npos) {
        ASSERT_PRED_FORMAT2(::testing::IsSubstring, "is not a valid size", msg);
      } else
#ifdef KOKKOS_ENABLE_SYCL
          if (msg.find("insufficient memory") != std::string::npos)
#endif
      {
        ASSERT_PRED_FORMAT2(::testing::IsSubstring, "insufficient memory", msg);
      }
      // SYCL cannot tell the reason why a memory allocation failed
#ifdef KOKKOS_ENABLE_SYCL
      else {
        // Otherwise, there has to be some sort of "unknown error" error
        ASSERT_PRED_FORMAT2(::testing::IsSubstring,
                            "because of an unknown error.", msg);
      }
#endif
    }
  }
};

}  // namespace Test
