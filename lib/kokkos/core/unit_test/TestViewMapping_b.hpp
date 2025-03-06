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

#include <sstream>
#include <iostream>

#include <Kokkos_Core.hpp>

namespace Test {

/*--------------------------------------------------------------------------*/

template <class Space>
struct TestViewMappingAtomic {
  using ExecSpace = typename Space::execution_space;
  using MemSpace  = typename Space::memory_space;

  using mem_trait = Kokkos::MemoryTraits<Kokkos::Atomic>;

  using T      = Kokkos::View<int *, ExecSpace>;
  using T_atom = Kokkos::View<int *, ExecSpace, mem_trait>;

  T x;
  T_atom x_atom;

  enum { N = 100000 };

  struct TagInit {};
  struct TagUpdate {};
  struct TagVerify {};

  KOKKOS_INLINE_FUNCTION
  void operator()(const TagInit &, const int i) const { x(i) = i; }

  KOKKOS_INLINE_FUNCTION
  void operator()(const TagUpdate &, const int i) const { x_atom(i % 2) += 1; }

  KOKKOS_INLINE_FUNCTION
  void operator()(const TagVerify &, const int i, long &error_count) const {
    if (i < 2) {
      if (x(i) != int(i + N / 2)) ++error_count;
    } else {
      if (x(i) != int(i)) ++error_count;
    }
  }

  TestViewMappingAtomic() : x("x", N), x_atom(x) {}

  void run() {
    ASSERT_TRUE(T::reference_type_is_lvalue_reference);
    ASSERT_FALSE(T_atom::reference_type_is_lvalue_reference);

    Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace, TagInit>(0, N), *this);
    Kokkos::parallel_for(Kokkos::RangePolicy<ExecSpace, TagUpdate>(0, N),
                         *this);

    long error_count = -1;

    Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace, TagVerify>(0, N),
                            *this, error_count);

    ASSERT_EQ(0, error_count);

    typename T_atom::HostMirror x_host = Kokkos::create_mirror_view(x);
    Kokkos::deep_copy(x_host, x);

    error_count = -1;

    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace, TagVerify>(0, N),
        [=](const TagVerify &, const int i, long &tmp_error_count) {
          if (i < 2) {
            if (x_host(i) != int(i + N / 2)) ++tmp_error_count;
          } else {
            if (x_host(i) != int(i)) ++tmp_error_count;
          }
        },
        error_count);

    ASSERT_EQ(0, error_count);
    Kokkos::deep_copy(x, x_host);
  }
};

TEST(TEST_CATEGORY, view_mapping_atomic) {
  TestViewMappingAtomic<TEST_EXECSPACE> f;
  f.run();
}

}  // namespace Test

/*--------------------------------------------------------------------------*/

namespace Test {

struct MappingClassValueType {
  KOKKOS_INLINE_FUNCTION
  MappingClassValueType() {
#if 0
    KOKKOS_IF_ON_DEVICE(
        (printf("TestViewMappingClassValue construct on Device\n");))
    KOKKOS_IF_ON_HOST((printf("TestViewMappingClassValue construct on Host\n");))
#endif
  }
  KOKKOS_INLINE_FUNCTION
  ~MappingClassValueType() {
#if 0
    KOKKOS_IF_ON_DEVICE(
        (printf("TestViewMappingClassValue destruct on Device\n");))
    KOKKOS_IF_ON_HOST((printf("TestViewMappingClassValue destruct on Host\n");))
#endif
  }
};

template <class Space>
void test_view_mapping_class_value() {
  using ExecSpace = typename Space::execution_space;

  ExecSpace().fence();
  {
    Kokkos::View<MappingClassValueType, ExecSpace> a("a");
    ExecSpace().fence();
  }
  ExecSpace().fence();
}

TEST(TEST_CATEGORY, view_mapping_class_value) {
  test_view_mapping_class_value<TEST_EXECSPACE>();
}

}  // namespace Test

/*--------------------------------------------------------------------------*/

namespace Test {

TEST(TEST_CATEGORY, view_mapping_assignable) {
  using exec_space = TEST_EXECSPACE;

  {  // Assignment of rank-0 Left = Right
    using dst_traits = Kokkos::ViewTraits<int, Kokkos::LayoutLeft, exec_space>;
    using src_traits = Kokkos::ViewTraits<int, Kokkos::LayoutRight, exec_space>;
    using mapping    = Kokkos::Impl::ViewMapping<dst_traits, src_traits, void>;
    static_assert(mapping::is_assignable);

    Kokkos::View<int, Kokkos::LayoutRight, exec_space> src;
    Kokkos::View<int, Kokkos::LayoutLeft, exec_space> dst(src);
    dst = src;
  }

  {  // Assignment of rank-0 Right = Left
    using dst_traits = Kokkos::ViewTraits<int, Kokkos::LayoutRight, exec_space>;
    using src_traits = Kokkos::ViewTraits<int, Kokkos::LayoutLeft, exec_space>;
    using mapping    = Kokkos::Impl::ViewMapping<dst_traits, src_traits, void>;
    static_assert(mapping::is_assignable);

    Kokkos::View<int, Kokkos::LayoutLeft, exec_space> src;
    Kokkos::View<int, Kokkos::LayoutRight, exec_space> dst(src);
    dst = src;
  }

  {  // Assignment of rank-1 Left = Right
    using dst_traits =
        Kokkos::ViewTraits<int *, Kokkos::LayoutLeft, exec_space>;
    using src_traits =
        Kokkos::ViewTraits<int *, Kokkos::LayoutRight, exec_space>;
    using mapping = Kokkos::Impl::ViewMapping<dst_traits, src_traits, void>;
    static_assert(mapping::is_assignable);

    Kokkos::View<int *, Kokkos::LayoutRight, exec_space> src;
    Kokkos::View<int *, Kokkos::LayoutLeft, exec_space> dst(src);
    dst = src;
  }

  {  // Assignment of rank-1 Right = Left
    using dst_traits =
        Kokkos::ViewTraits<int *, Kokkos::LayoutRight, exec_space>;
    using src_traits =
        Kokkos::ViewTraits<int *, Kokkos::LayoutLeft, exec_space>;
    using mapping = Kokkos::Impl::ViewMapping<dst_traits, src_traits, void>;
    static_assert(mapping::is_assignable);

    Kokkos::View<int *, Kokkos::LayoutLeft, exec_space> src;
    Kokkos::View<int *, Kokkos::LayoutRight, exec_space> dst(src);
    dst = src;
  }

  {  // Assignment of rank-2 Left = Right
    using dst_traits =
        Kokkos::ViewTraits<int **, Kokkos::LayoutLeft, exec_space>;
    using src_traits =
        Kokkos::ViewTraits<int **, Kokkos::LayoutRight, exec_space>;
    using mapping = Kokkos::Impl::ViewMapping<dst_traits, src_traits, void>;
    static_assert(!mapping::is_assignable);
  }

  {  // Assignment of rank-2 Right = Left
    using dst_traits =
        Kokkos::ViewTraits<int **, Kokkos::LayoutRight, exec_space>;
    using src_traits =
        Kokkos::ViewTraits<int **, Kokkos::LayoutLeft, exec_space>;
    using mapping = Kokkos::Impl::ViewMapping<dst_traits, src_traits, void>;
    static_assert(!mapping::is_assignable);
  }
}

TEST(TEST_CATEGORY, view_mapping_trivially_copyable) {
  using exec_space = TEST_EXECSPACE;

  using dst_traits = Kokkos::ViewTraits<int *, exec_space>;
  using src_traits = dst_traits;
  using mapping    = Kokkos::Impl::ViewMapping<dst_traits, src_traits, void>;

  static_assert(std::is_trivially_copyable<mapping>{});
}

}  // namespace Test
