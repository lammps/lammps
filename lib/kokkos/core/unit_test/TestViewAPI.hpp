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
#include <stdexcept>
#include <sstream>
#include <iostream>

namespace Test {

template <class T, class... P>
size_t allocation_count(const Kokkos::View<T, P...> &view) {
  const size_t card  = view.size();
  const size_t alloc = view.span();

  const int memory_span = Kokkos::View<int *>::required_allocation_size(100);

  return (card <= alloc && memory_span == 400) ? alloc : 0;
}

/*--------------------------------------------------------------------------*/

template <typename T, class DeviceType>
struct TestViewOperator {
  typedef typename DeviceType::execution_space execution_space;

  enum { N = 1000 };
  enum { D = 3 };

  typedef Kokkos::View<T * [D], execution_space> view_type;

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
  typedef typename DeviceType::execution_space execution_space;
  typedef typename DeviceType::memory_space memory_space;
  typedef typename execution_space::size_type size_type;

  typedef int value_type;

  KOKKOS_INLINE_FUNCTION
  static void join(volatile value_type &update,
                   const volatile value_type &input) {
    update |= input;
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type &update) { update = 0; }

  typedef Kokkos::View<DataType, Kokkos::LayoutLeft, execution_space> left_view;
  typedef Kokkos::View<DataType, Kokkos::LayoutRight, execution_space>
      right_view;
  typedef Kokkos::View<DataType, Kokkos::LayoutStride, execution_space>
      stride_view;

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
  typedef typename DeviceType::execution_space execution_space;
  typedef typename DeviceType::memory_space memory_space;
  typedef typename execution_space::size_type size_type;

  typedef int value_type;

  KOKKOS_INLINE_FUNCTION
  static void join(volatile value_type &update,
                   const volatile value_type &input) {
    update |= input;
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type &update) { update = 0; }

  typedef Kokkos::View<DataType, Kokkos::LayoutLeft, execution_space> left_view;
  typedef Kokkos::View<DataType, Kokkos::LayoutRight, execution_space>
      right_view;

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
  typedef typename DeviceType::execution_space execution_space;
  typedef typename DeviceType::memory_space memory_space;
  typedef typename execution_space::size_type size_type;

  typedef int value_type;

  KOKKOS_INLINE_FUNCTION
  static void join(volatile value_type &update,
                   const volatile value_type &input) {
    update |= input;
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type &update) { update = 0; }

  typedef Kokkos::View<DataType, Kokkos::LayoutLeft, execution_space> left_view;
  typedef Kokkos::View<DataType, Kokkos::LayoutRight, execution_space>
      right_view;

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
  typedef typename DeviceType::execution_space execution_space;
  typedef typename DeviceType::memory_space memory_space;
  typedef typename execution_space::size_type size_type;

  typedef int value_type;

  KOKKOS_INLINE_FUNCTION
  static void join(volatile value_type &update,
                   const volatile value_type &input) {
    update |= input;
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type &update) { update = 0; }

  typedef Kokkos::View<DataType, Kokkos::LayoutLeft, execution_space> left_view;
  typedef Kokkos::View<DataType, Kokkos::LayoutRight, execution_space>
      right_view;
  typedef Kokkos::View<DataType, Kokkos::LayoutStride, execution_space>
      stride_view;

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
  typedef typename DeviceType::execution_space execution_space;
  typedef typename DeviceType::memory_space memory_space;
  typedef typename execution_space::size_type size_type;

  typedef int value_type;

  KOKKOS_INLINE_FUNCTION
  static void join(volatile value_type &update,
                   const volatile value_type &input) {
    update |= input;
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type &update) { update = 0; }

  typedef Kokkos::View<DataType, Kokkos::LayoutLeft, execution_space> left_view;
  typedef Kokkos::View<DataType, Kokkos::LayoutRight, execution_space>
      right_view;

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
  typedef typename DeviceType::execution_space execution_space;
  typedef typename DeviceType::memory_space memory_space;
  typedef typename execution_space::size_type size_type;

  typedef int value_type;

  KOKKOS_INLINE_FUNCTION
  static void join(volatile value_type &update,
                   const volatile value_type &input) {
    update |= input;
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type &update) { update = 0; }

  typedef Kokkos::View<DataType, Kokkos::LayoutLeft, execution_space> left_view;
  typedef Kokkos::View<DataType, Kokkos::LayoutRight, execution_space>
      right_view;
  typedef Kokkos::View<DataType, Kokkos::LayoutStride, execution_space>
      stride_view;

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
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
          if (&left(i0, i1, i2) != &left(i0, i1, i2, 0, 0, 0, 0, 0)) {
            update |= 3;
          }
          if (&right(i0, i1, i2) != &right(i0, i1, i2, 0, 0, 0, 0, 0)) {
            update |= 3;
          }
#else
          if (&left(i0, i1, i2) != &left.access(i0, i1, i2, 0, 0, 0, 0, 0)) {
            update |= 3;
          }
          if (&right(i0, i1, i2) != &right.access(i0, i1, i2, 0, 0, 0, 0, 0)) {
            update |= 3;
          }
#endif
        }
  }
};

template <class DataType, class DeviceType>
struct TestViewOperator_LeftAndRight<DataType, DeviceType, 2> {
  typedef typename DeviceType::execution_space execution_space;
  typedef typename DeviceType::memory_space memory_space;
  typedef typename execution_space::size_type size_type;

  typedef int value_type;

  KOKKOS_INLINE_FUNCTION
  static void join(volatile value_type &update,
                   const volatile value_type &input) {
    update |= input;
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type &update) { update = 0; }

  typedef Kokkos::View<DataType, Kokkos::LayoutLeft, execution_space> left_view;
  typedef Kokkos::View<DataType, Kokkos::LayoutRight, execution_space>
      right_view;

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
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
        if (&left(i0, i1) != &left(i0, i1, 0, 0, 0, 0, 0, 0)) {
          update |= 3;
        }
        if (&right(i0, i1) != &right(i0, i1, 0, 0, 0, 0, 0, 0)) {
          update |= 3;
        }
#else
        if (&left(i0, i1) != &left.access(i0, i1, 0, 0, 0, 0, 0, 0)) {
          update |= 3;
        }
        if (&right(i0, i1) != &right.access(i0, i1, 0, 0, 0, 0, 0, 0)) {
          update |= 3;
        }
#endif
      }
  }
};

template <class DataType, class DeviceType>
struct TestViewOperator_LeftAndRight<DataType, DeviceType, 1> {
  typedef typename DeviceType::execution_space execution_space;
  typedef typename DeviceType::memory_space memory_space;
  typedef typename execution_space::size_type size_type;

  typedef int value_type;

  KOKKOS_INLINE_FUNCTION
  static void join(volatile value_type &update,
                   const volatile value_type &input) {
    update |= input;
  }

  KOKKOS_INLINE_FUNCTION
  static void init(value_type &update) { update = 0; }

  typedef Kokkos::View<DataType, Kokkos::LayoutLeft, execution_space> left_view;
  typedef Kokkos::View<DataType, Kokkos::LayoutRight, execution_space>
      right_view;
  typedef Kokkos::View<DataType, Kokkos::LayoutStride, execution_space>
      stride_view;

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
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
      if (&left(i0) != &left(i0, 0, 0, 0, 0, 0, 0, 0)) {
        update |= 3;
      }
      if (&right(i0) != &right(i0, 0, 0, 0, 0, 0, 0, 0)) {
        update |= 3;
      }
#else
      if (&left(i0) != &left.access(i0, 0, 0, 0, 0, 0, 0, 0)) {
        update |= 3;
      }
      if (&right(i0) != &right.access(i0, 0, 0, 0, 0, 0, 0, 0)) {
        update |= 3;
      }
#endif
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

  template <class MemoryTraits, class Space>
  struct CopyUnInit {
    typedef typename Kokkos::Impl::MirrorViewType<
        Space, double *, Layout, Kokkos::HostSpace, MemoryTraits>::view_type
        mirror_view_type;

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
    auto a_d = Kokkos::create_mirror_view(DeviceType(), a_h,
                                          Kokkos::WithoutInitializing);

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
    test_mirror_no_initialize<Kokkos::MemoryTraits<0> >();
    test_mirror_no_initialize<Kokkos::MemoryTraits<Kokkos::Unmanaged> >();
  }
};

/*--------------------------------------------------------------------------*/

template <typename T, class DeviceType>
class TestViewAPI {
 public:
  typedef DeviceType device;

  enum { N0 = 1000, N1 = 3, N2 = 5, N3 = 7 };

  typedef Kokkos::View<T, device> dView0;
  typedef Kokkos::View<T *, device> dView1;
  typedef Kokkos::View<T * [N1], device> dView2;
  typedef Kokkos::View<T * [N1][N2], device> dView3;
  typedef Kokkos::View<T * [N1][N2][N3], device> dView4;
  typedef Kokkos::View<const T * [N1][N2][N3], device> const_dView4;
  typedef Kokkos::View<T ****, device, Kokkos::MemoryUnmanaged>
      dView4_unmanaged;
  typedef typename dView0::host_mirror_space host;

  TestViewAPI() {}

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
    typedef Kokkos::View<int, host> view_type;
    typedef typename view_type::HostMirror mirror_type;

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
    typedef typename dView0::HostMirror hView0;

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

  static void run_test() {
    // mfh 14 Feb 2014: This test doesn't actually create instances of
    // these types.  In order to avoid "declared but unused typedef"
    // warnings, we declare empty instances of these types, with the
    // usual "(void)" marker to avoid compiler warnings for unused
    // variables.

    typedef typename dView0::HostMirror hView0;
    typedef typename dView1::HostMirror hView1;
    typedef typename dView2::HostMirror hView2;
    typedef typename dView3::HostMirror hView3;
    typedef typename dView4::HostMirror hView4;

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

    ASSERT_TRUE(dx.data() == 0);
    ASSERT_TRUE(dy.data() == 0);
    ASSERT_TRUE(dz.data() == 0);
    ASSERT_TRUE(hx.data() == 0);
    ASSERT_TRUE(hy.data() == 0);
    ASSERT_TRUE(hz.data() == 0);
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

    ASSERT_EQ(dx.use_count(), size_t(1));

    dView4_unmanaged unmanaged_dx = dx;
    ASSERT_EQ(dx.use_count(), size_t(1));

    dView4_unmanaged unmanaged_from_ptr_dx = dView4_unmanaged(
        dx.data(), dx.extent(0), dx.extent(1), dx.extent(2), dx.extent(3));

    {
      // Destruction of this view should be harmless.

      const_dView4 unmanaged_from_ptr_const_dx(dx.data(), dx.extent(0)
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
                                                              ,
                                               dx.extent(1), dx.extent(2),
                                               dx.extent(3)
#endif
      );
    }

    const_dView4 const_dx = dx;
    ASSERT_EQ(dx.use_count(), size_t(2));

    {
      const_dView4 const_dx2;
      const_dx2 = const_dx;
      ASSERT_EQ(dx.use_count(), size_t(3));

      const_dx2 = dy;
      ASSERT_EQ(dx.use_count(), size_t(2));

      const_dView4 const_dx3(dx);
      ASSERT_EQ(dx.use_count(), size_t(3));

      dView4_unmanaged dx4_unmanaged(dx);
      ASSERT_EQ(dx.use_count(), size_t(3));
    }

    ASSERT_EQ(dx.use_count(), size_t(2));

    ASSERT_FALSE(dx.data() == 0);
    ASSERT_FALSE(const_dx.data() == 0);
    ASSERT_FALSE(unmanaged_dx.data() == 0);
    ASSERT_FALSE(unmanaged_from_ptr_dx.data() == 0);
    ASSERT_FALSE(dy.data() == 0);
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

      for (size_t ip = 0; ip < N0; ++ip)
        for (size_t i1 = 0; i1 < N1; ++i1)
          for (size_t i2 = 0; i2 < N2; ++i2)
            for (size_t i3 = 0; i3 < N3; ++i3) {
              ASSERT_EQ(hx(ip, i1, i2, i3), hy(ip, i1, i2, i3));
            }

      Kokkos::deep_copy(typename hView4::execution_space(), dx, T(0));
      Kokkos::deep_copy(typename hView4::execution_space(), hx, dx);

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

      for (size_t ip = 0; ip < N0; ++ip)
        for (size_t i1 = 0; i1 < N1; ++i1)
          for (size_t i2 = 0; i2 < N2; ++i2)
            for (size_t i3 = 0; i3 < N3; ++i3) {
              ASSERT_EQ(hx(ip, i1, i2, i3), hy(ip, i1, i2, i3));
            }

      Kokkos::deep_copy(typename dView4::execution_space(), dx, T(0));
      Kokkos::deep_copy(typename dView4::execution_space(), hx, dx);

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
    ASSERT_TRUE(dx.data() == 0);
    ASSERT_FALSE(dy.data() == 0);
    ASSERT_FALSE(dz.data() == 0);

    dy = dView4();
    ASSERT_TRUE(dx.data() == 0);
    ASSERT_TRUE(dy.data() == 0);
    ASSERT_FALSE(dz.data() == 0);

    dz = dView4();
    ASSERT_TRUE(dx.data() == 0);
    ASSERT_TRUE(dy.data() == 0);
    ASSERT_TRUE(dz.data() == 0);
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

  typedef T DataType[2];

  static void check_auto_conversion_to_const(
      const Kokkos::View<const DataType, device> &arg_const,
      const Kokkos::View<DataType, device> &arg) {
    ASSERT_TRUE(arg_const == arg);
  }

  static void run_test_const() {
    typedef Kokkos::View<DataType, device> typeX;
    typedef Kokkos::View<const DataType, device> const_typeX;
    typedef Kokkos::View<const DataType, device, Kokkos::MemoryRandomAccess>
        const_typeR;

    typeX x("X");
    const_typeX xc = x;
    const_typeR xr = x;

    ASSERT_TRUE(xc == x);
    ASSERT_TRUE(x == xc);

    // For CUDA the constant random access View does not return
    // an lvalue reference due to retrieving through texture cache
    // therefore not allowed to query the underlying pointer.
#if defined(KOKKOS_ENABLE_CUDA)
    if (!std::is_same<typename device::execution_space, Kokkos::Cuda>::value)
#endif
    {
      ASSERT_TRUE(x.data() == xr.data());
    }

    // typeX xf = xc; // Setting non-const from const must not compile.

    check_auto_conversion_to_const(x, x);
  }

  static void run_test_subview() {
    typedef Kokkos::View<const T, device> sView;

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
    typedef Kokkos::View<int ****, Kokkos::LayoutLeft, host> view_left_4;
    typedef Kokkos::View<int ****, Kokkos::LayoutRight, host> view_right_4;
    typedef Kokkos::View<int **, Kokkos::LayoutLeft, host> view_left_2;
    typedef Kokkos::View<int **, Kokkos::LayoutRight, host> view_right_2;

    typedef Kokkos::View<int *, Kokkos::LayoutStride, host> view_stride_1;
    typedef Kokkos::View<int **, Kokkos::LayoutStride, host> view_stride_2;

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

    typedef Kokkos::View<T *, Kokkos::LayoutLeft, host> vector_type;
    typedef Kokkos::View<T **, Kokkos::LayoutLeft, host> multivector_type;

    typedef Kokkos::View<T *, Kokkos::LayoutRight, host> vector_right_type;
    typedef Kokkos::View<T **, Kokkos::LayoutRight, host>
        multivector_right_type;

    typedef Kokkos::View<const T *, Kokkos::LayoutRight, host>
        const_vector_right_type;
    typedef Kokkos::View<const T *, Kokkos::LayoutLeft, host> const_vector_type;
    typedef Kokkos::View<const T **, Kokkos::LayoutLeft, host>
        const_multivector_type;

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

    ASSERT_TRUE(&v1[0] == &v1(0));
    ASSERT_TRUE(&v1[0] == &mv(0, 0));
    ASSERT_TRUE(&v2[0] == &mv(0, 1));
    ASSERT_TRUE(&v3[0] == &mv(0, 2));

    ASSERT_TRUE(&cv1[0] == &mv(0, 0));
    ASSERT_TRUE(&cv2[0] == &mv(0, 1));
    ASSERT_TRUE(&cv3[0] == &mv(0, 2));

    ASSERT_TRUE(&vr1[0] == &mv(0, 0));
    ASSERT_TRUE(&vr2[0] == &mv(0, 1));
    ASSERT_TRUE(&vr3[0] == &mv(0, 2));

    ASSERT_TRUE(&cvr1[0] == &mv(0, 0));
    ASSERT_TRUE(&cvr2[0] == &mv(0, 1));
    ASSERT_TRUE(&cvr3[0] == &mv(0, 2));

    ASSERT_TRUE(&mv1(0, 0) == &mv(1, 2));
    ASSERT_TRUE(&mv1(1, 1) == &mv(2, 3));
    ASSERT_TRUE(&mv1(3, 2) == &mv(4, 4));
    ASSERT_TRUE(&mvr1(0, 0) == &mv_right(1, 2));
    ASSERT_TRUE(&mvr1(1, 1) == &mv_right(2, 3));
    ASSERT_TRUE(&mvr1(3, 2) == &mv_right(4, 4));

    const_vector_type c_cv1(v1);
    typename vector_type::const_type c_cv2(v2);
    typename const_vector_type::const_type c_ccv2(v2);

    const_multivector_type cmv(mv);
    typename multivector_type::const_type cmvX(cmv);
    typename const_multivector_type::const_type ccmvX(cmv);
  }

  static void run_test_error() {
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
      } else {
        // Otherwise, there has to be some sort of "insufficient memory" error
        ASSERT_PRED_FORMAT2(::testing::IsSubstring, "insufficient memory", msg);
      }
    }
  }
};

}  // namespace Test
