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

#ifndef KOKKOS_TEST_DUALVIEW_HPP
#define KOKKOS_TEST_DUALVIEW_HPP

#include <gtest/gtest.h>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <Kokkos_Timer.hpp>
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
    result =
        run_me<Kokkos::DualView<Scalar**, Kokkos::LayoutLeft, Device>>(size, 3);
  }
};

template <typename Scalar, class Device>
struct test_dualview_copy_construction_and_assignment {
  using scalar_type     = Scalar;
  using execution_space = Device;

  void operator()() {
    constexpr unsigned int n = 10;
    constexpr unsigned int m = 5;

    using SrcViewType = Kokkos::DualView<Scalar**, Kokkos::LayoutLeft, Device>;
    using DstViewType =
        Kokkos::DualView<const Scalar* [m], Kokkos::LayoutLeft, Device>;

    SrcViewType a("A", n, m);

    // Copy construction
    DstViewType b(a);

    // Copy assignment
    DstViewType c = a;

    // Check equality (shallow) of the host and device views
    ASSERT_EQ(a.view_host(), b.view_host());
    ASSERT_EQ(a.view_device(), b.view_device());

    ASSERT_EQ(a.view_host(), c.view_host());
    ASSERT_EQ(a.view_device(), c.view_device());

    // We can't test shallow equality of modified_flags because it's protected.
    // So we test it indirectly through sync state behavior.
    if (!SrcViewType::impl_dualview_is_single_device) {
      a.clear_sync_state();
      a.modify_host();
      ASSERT_TRUE(a.need_sync_device());
      ASSERT_TRUE(b.need_sync_device());
      ASSERT_TRUE(c.need_sync_device());
      a.clear_sync_state();
    }
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
    Kokkos::deep_copy(a.view_device(), 1);

    a.template modify<typename ViewType::execution_space>();
    a.template sync<typename ViewType::host_mirror_space>();
    a.template sync<typename ViewType::host_mirror_space>(
        Kokkos::DefaultExecutionSpace{});

    a.view_host()(5, 1) = 3;
    a.view_host()(6, 1) = 4;
    a.view_host()(7, 2) = 5;
    a.template modify<typename ViewType::host_mirror_space>();
    ViewType b = Kokkos::subview(a, std::pair<unsigned int, unsigned int>(6, 9),
                                 std::pair<unsigned int, unsigned int>(0, 1));
    a.template sync<typename ViewType::execution_space>();
    a.template sync<typename ViewType::execution_space>(
        Kokkos::DefaultExecutionSpace{});
    b.template modify<typename ViewType::execution_space>();

    Kokkos::deep_copy(b.view_device(), 2);

    a.template sync<typename ViewType::host_mirror_space>();
    a.template sync<typename ViewType::host_mirror_space>(
        Kokkos::DefaultExecutionSpace{});
    Scalar count = 0;
    for (unsigned int i = 0; i < a.view_device().extent(0); i++)
      for (unsigned int j = 0; j < a.view_device().extent(1); j++)
        count += a.view_host()(i, j);
    return count - a.view_device().extent(0) * a.view_device().extent(1) - 2 -
           4 - 3 * 2;
  }

  test_dualview_combinations(unsigned int size, bool with_init) {
    result = run_me<Kokkos::DualView<Scalar**, Kokkos::LayoutLeft, Device>>(
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

    Kokkos::deep_copy(a.view_device(), 1);

    if (use_templ_sync) {
      a.template modify<typename ViewType::execution_space>();
      a.template sync<typename ViewType::host_mirror_space>();
    } else {
      a.modify_device();
      a.sync_host();
      a.sync_host(Kokkos::DefaultExecutionSpace{});
    }

    // Check device view is initialized as expected
    scalar_type a_d_sum = 0;
    // Execute on the execution_space associated with t_dev's memory space
    using t_dev_exec_space =
        typename ViewType::t_dev::memory_space::execution_space;
    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<t_dev_exec_space>(0, n),
        SumViewEntriesFunctor<scalar_type, typename ViewType::t_dev>(
            a.view_device()),
        a_d_sum);
    ASSERT_EQ(a_d_sum, sum_total);

    // Check host view is synced as expected
    scalar_type a_h_sum = 0;
    for (size_t i = 0; i < a.view_host().extent(0); ++i)
      for (size_t j = 0; j < a.view_host().extent(1); ++j) {
        a_h_sum += a.view_host()(i, j);
      }

    ASSERT_EQ(a_h_sum, sum_total);

    // Test deep_copy
    Kokkos::deep_copy(b, a);
    if (use_templ_sync) {
      b.template sync<typename ViewType::host_mirror_space>();
    } else {
      b.sync_host();
      b.sync_host(Kokkos::DefaultExecutionSpace{});
    }

    // Perform same checks on b as done on a
    // Check device view is initialized as expected
    scalar_type b_d_sum = 0;
    // Execute on the execution_space associated with t_dev's memory space
    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<t_dev_exec_space>(0, n),
        SumViewEntriesFunctor<scalar_type, typename ViewType::t_dev>(
            b.view_device()),
        b_d_sum);
    ASSERT_EQ(b_d_sum, sum_total);

    // Check host view is synced as expected
    scalar_type b_h_sum = 0;
    for (size_t i = 0; i < b.view_host().extent(0); ++i)
      for (size_t j = 0; j < b.view_host().extent(1); ++j) {
        b_h_sum += b.view_host()(i, j);
      }

    ASSERT_EQ(b_h_sum, sum_total);

  }  // end run_me

  test_dual_view_deep_copy() {
    run_me<Kokkos::DualView<Scalar**, Kokkos::LayoutLeft, Device>>(10, 5, true);
    run_me<Kokkos::DualView<Scalar**, Kokkos::LayoutLeft, Device>>(10, 5,
                                                                   false);
    // Test zero length but allocated (a.view_device().data() != nullptr but
    // a.view_device().span() == 0)
    run_me<Kokkos::DualView<Scalar**, Kokkos::LayoutLeft, Device>>(0, 5, true);
    run_me<Kokkos::DualView<Scalar**, Kokkos::LayoutLeft, Device>>(0, 5, false);

    // Test default constructed view
    run_me<Kokkos::DualView<Scalar**, Kokkos::LayoutLeft, Device>>(-1, 5, true);
    run_me<Kokkos::DualView<Scalar**, Kokkos::LayoutLeft, Device>>(-1, 5,
                                                                   false);
  }
};

template <typename Scalar, class Device, bool Initialize>
struct test_dualview_resize {
  using scalar_type     = Scalar;
  using execution_space = Device;

  template <typename ViewType>
  void run_me() {
    const unsigned int n      = 10;
    const unsigned int m      = 5;
    const unsigned int factor = 2;

    ViewType a;
    if constexpr (Initialize)
      a = ViewType("A", n, m);
    else
      a = ViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing, "A"), n, m);

    Kokkos::deep_copy(a.view_device(), 1);

    /* Covers case "Resize on Device" */
    a.modify_device();
    if constexpr (Initialize)
      Kokkos::resize(a, factor * n, factor * m);
    else
      Kokkos::resize(Kokkos::WithoutInitializing, a, factor * n, factor * m);
    ASSERT_EQ(a.extent(0), n * factor);
    ASSERT_EQ(a.extent(1), m * factor);

    Kokkos::deep_copy(a.view_device(), 1);
    a.sync_host();

    // Check device view is initialized as expected
    // Execute on the execution_space associated with t_dev's memory space
    using t_dev_exec_space =
        typename ViewType::t_dev::memory_space::execution_space;
    Kokkos::View<int, typename ViewType::t_dev::memory_space> errors_d(
        "errors");
    Kokkos::parallel_for(
        Kokkos::MDRangePolicy<t_dev_exec_space, Kokkos::Rank<2>>(
            {0, 0}, {a.view_device().extent(0), a.view_device().extent(1)}),
        KOKKOS_LAMBDA(int i, int j) {
          if (a.view_device()(i, j) != 1) Kokkos::atomic_inc(errors_d.data());
        });
    int errors_d_scalar;
    Kokkos::deep_copy(errors_d_scalar, errors_d);

    // Check host view is synced as expected
    int errors_h_scalar = 0;
    for (size_t i = 0; i < a.view_host().extent(0); ++i)
      for (size_t j = 0; j < a.view_host().extent(1); ++j) {
        if (a.view_host()(i, j) != 1) ++errors_h_scalar;
      }

    // Check
    ASSERT_EQ(errors_d_scalar, 0);
    ASSERT_EQ(errors_h_scalar, 0);

    /* Covers case "Resize on Host" */
    a.modify_host();

    if constexpr (Initialize)
      Kokkos::resize(a, n / factor, m / factor);
    else
      Kokkos::resize(Kokkos::WithoutInitializing, a, n / factor, m / factor);
    ASSERT_EQ(a.extent(0), n / factor);
    ASSERT_EQ(a.extent(1), m / factor);

    a.sync_device();
    a.sync_device(Kokkos::DefaultExecutionSpace{});

    // Check device view is initialized as expected
    Kokkos::deep_copy(errors_d, 0);
    // Execute on the execution_space associated with t_dev's memory space
    using t_dev_exec_space =
        typename ViewType::t_dev::memory_space::execution_space;
    Kokkos::parallel_for(
        Kokkos::MDRangePolicy<t_dev_exec_space, Kokkos::Rank<2>>(
            {0, 0}, {a.view_device().extent(0), a.view_device().extent(1)}),
        KOKKOS_LAMBDA(int i, int j) {
          if (a.view_device()(i, j) != 1) Kokkos::atomic_inc(errors_d.data());
        });
    Kokkos::deep_copy(errors_d_scalar, errors_d);

    // Check host view is synced as expected
    errors_h_scalar = 0;
    for (size_t i = 0; i < a.view_host().extent(0); ++i)
      for (size_t j = 0; j < a.view_host().extent(1); ++j) {
        if (a.view_host()(i, j) != 1) ++errors_h_scalar;
      }

    // Check
    ASSERT_EQ(errors_d_scalar, 0);
    ASSERT_EQ(errors_h_scalar, 0);

  }  // end run_me

  test_dualview_resize() {
    run_me<Kokkos::DualView<Scalar**, Kokkos::LayoutLeft, Device>>();
  }
};

template <typename Scalar, class Device, bool Initialize>
struct test_dualview_realloc {
  using scalar_type     = Scalar;
  using execution_space = Device;

  template <typename ViewType>
  void run_me() {
    const unsigned int n = 10;
    const unsigned int m = 5;

    ViewType a;
    if constexpr (Initialize) {
      a = ViewType("A", n, m);
      Kokkos::realloc(a, n, m);
    } else {
      a = ViewType(Kokkos::view_alloc(Kokkos::WithoutInitializing, "A"), n, m);
      Kokkos::realloc(Kokkos::WithoutInitializing, a, n, m);
    }
    ASSERT_EQ(a.extent(0), n);
    ASSERT_EQ(a.extent(1), m);

    Kokkos::deep_copy(a.view_device(), 1);

    a.modify_device();
    a.sync_host();

    // Check device view is initialized as expected
    // Execute on the execution_space associated with t_dev's memory space
    using t_dev_exec_space =
        typename ViewType::t_dev::memory_space::execution_space;
    Kokkos::View<int, typename ViewType::t_dev::memory_space> errors_d(
        "errors");
    Kokkos::parallel_for(
        Kokkos::MDRangePolicy<t_dev_exec_space, Kokkos::Rank<2>>(
            {0, 0}, {a.view_device().extent(0), a.view_device().extent(1)}),
        KOKKOS_LAMBDA(int i, int j) {
          if (a.view_device()(i, j) != 1) Kokkos::atomic_inc(errors_d.data());
        });
    int errors_d_scalar;
    Kokkos::deep_copy(errors_d_scalar, errors_d);

    // Check host view is synced as expected
    int errors_h_scalar = 0;
    for (size_t i = 0; i < a.view_host().extent(0); ++i)
      for (size_t j = 0; j < a.view_host().extent(1); ++j) {
        if (a.view_host()(i, j) != 1) ++errors_h_scalar;
      }

    // Check
    ASSERT_EQ(errors_d_scalar, 0);
    ASSERT_EQ(errors_h_scalar, 0);
  }  // end run_me

  test_dualview_realloc() {
    run_me<Kokkos::DualView<Scalar**, Kokkos::LayoutLeft, Device>>();
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
void test_dualview_copy_construction_and_assignment() {
  Impl::test_dualview_copy_construction_and_assignment<Scalar, Device>()();
}

template <typename Scalar, typename Device>
void test_dualview_deep_copy() {
  Impl::test_dual_view_deep_copy<Scalar, Device>();
}

template <typename Scalar, typename Device>
void test_dualview_realloc() {
  Impl::test_dualview_realloc<Scalar, Device, false>();
  Impl::test_dualview_realloc<Scalar, Device, true>();
}

template <typename Scalar, typename Device>
void test_dualview_resize() {
  Impl::test_dualview_resize<Scalar, Device, false>();
  Impl::test_dualview_resize<Scalar, Device, true>();
}

TEST(TEST_CATEGORY, dualview_combination) {
  test_dualview_combinations<int, TEST_EXECSPACE>(10, true);
}

TEST(TEST_CATEGORY, dualview_alloc) {
  test_dualview_alloc<int, TEST_EXECSPACE>(10);
}

TEST(TEST_CATEGORY, test_dualview_copy_construction_and_assignment) {
  test_dualview_copy_construction_and_assignment<int, TEST_EXECSPACE>();
}

TEST(TEST_CATEGORY, dualview_combinations_without_init) {
  test_dualview_combinations<int, TEST_EXECSPACE>(10, false);
}

TEST(TEST_CATEGORY, dualview_deep_copy) {
  test_dualview_deep_copy<int, TEST_EXECSPACE>();
  test_dualview_deep_copy<double, TEST_EXECSPACE>();
}

template <typename ExecutionSpace>
void test_dualview_sync_should_fence() {
  using DualViewType = Kokkos::DualView<int, ExecutionSpace>;
  {
    DualViewType dv("test_dual_view");
    dv.modify_device();
    Kokkos::parallel_for(
        Kokkos::RangePolicy<ExecutionSpace>(0, 10000),
        KOKKOS_LAMBDA(int) { Kokkos::atomic_add(dv.view_device().data(), 1); });
    dv.sync_host();
    ASSERT_EQ(dv.view_host()(), 10000);
  }
  {
    DualViewType dv("test_dual_view");
    dv.modify_device();
    Kokkos::parallel_for(
        Kokkos::RangePolicy<ExecutionSpace>(0, 10000),
        KOKKOS_LAMBDA(int) { Kokkos::atomic_add(dv.view_device().data(), 1); });
    dv.template sync<typename DualViewType::t_host::device_type>();
    ASSERT_EQ(dv.view_host()(), 10000);
  }
  {
    DualViewType dv("test_dual_view");
    dv.modify_host();
    Kokkos::parallel_for(
        Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, 10000),
        KOKKOS_LAMBDA(int) { Kokkos::atomic_add(dv.view_host().data(), 1); });
    dv.sync_device();
    int result;
    auto device_exec =
        Kokkos::Experimental::partition_space(ExecutionSpace{}, 1);
    Kokkos::deep_copy(device_exec[0], result, dv.view_device());
    device_exec[0].fence();
    ASSERT_EQ(result, 10000);
  }
  {
    DualViewType dv("test_dual_view");
    dv.modify_host();
    Kokkos::parallel_for(
        Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace>(0, 10000),
        KOKKOS_LAMBDA(int) { Kokkos::atomic_add(dv.view_host().data(), 1); });
    dv.template sync<typename DualViewType::t_dev::device_type>();
    int result;
    auto device_exec =
        Kokkos::Experimental::partition_space(ExecutionSpace{}, 1);
    Kokkos::deep_copy(device_exec[0], result, dv.view_device());
    device_exec[0].fence();
    ASSERT_EQ(result, 10000);
  }
}

TEST(TEST_CATEGORY, dualview_sync_should_fence) {
#ifdef KOKKOS_ENABLE_HPX  // FIXME_DUALVIEW_ASYNCHRONOUS_BACKENDS
  GTEST_SKIP() << "Known to fail with HPX";
#endif
  test_dualview_sync_should_fence<TEST_EXECSPACE>();
}

struct NoDefaultConstructor {
  NoDefaultConstructor(int i_) : i(i_) {}
  KOKKOS_FUNCTION operator int() const { return i; }

  int i;
};

TEST(TEST_CATEGORY, dualview_realloc) {
  test_dualview_realloc<int, TEST_EXECSPACE>();
  Impl::test_dualview_realloc<NoDefaultConstructor, TEST_EXECSPACE,
                              /* Initialize */ false>();
}

TEST(TEST_CATEGORY, dualview_resize) {
  test_dualview_resize<int, TEST_EXECSPACE>();
  Impl::test_dualview_resize<NoDefaultConstructor, TEST_EXECSPACE,
                             /* Initialize */ false>();
}

template <typename ExecutionSpace>
void check_dualview_external_view_construction() {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  Kokkos::View<int*, ExecutionSpace> view1("view1", 10);
  Kokkos::View<int*, ExecutionSpace> view2("view2", 10);

  Kokkos::DualView<int*, ExecutionSpace> v_dual(view1, view1);
  ASSERT_DEATH(
      (Kokkos::DualView<int*, ExecutionSpace>(view1, view2)),
      "DualView storing one View constructed from two different Views");
}

// FIXME_MSVC+CUDA error C2094: label 'gtest_label_520' was undefined
#if !(defined(KOKKOS_COMPILER_MSVC) && defined(KOKKOS_ENABLE_CUDA))
TEST(TEST_CATEGORY_DEATH, dualview_external_view_construction) {
  if constexpr (!Kokkos::SpaceAccessibility<
                    Kokkos::HostSpace,
                    TEST_EXECSPACE::memory_space>::accessible) {
    GTEST_SKIP() << "test only relevant if DualView uses one allocation";
  } else {
    // FIXME_CLANG We can't inline the function because recent clang versions
    // would deduce that a static_assert isn't satisfied for TEST_EXECSPACE.
    // Thus, we need to template the function on the execution space.
    check_dualview_external_view_construction<TEST_EXECSPACE>();
  }
}
#endif

namespace {
/**
 *
 * The following tests are a response to
 * https://github.com/kokkos/kokkos/issues/3850
 * and
 * https://github.com/kokkos/kokkos/pull/3857
 *
 * DualViews were returning incorrect view types and taking
 * inappropriate actions based on the templated view methods.
 *
 * Specifically, template view methods were always returning
 * a device view if the memory space was UVM and a Kokkos::Device was passed.
 * Sync/modify methods completely broke down So these tests exist to make sure
 * that we keep the semantics of UVM DualViews intact.
 */
// modify if we have other UVM enabled backends

#ifdef KOKKOS_HAS_SHARED_SPACE
template <typename ExecutionSpace>
using TestSharedSpace = Kokkos::SharedSpace;
#else
template <typename ExecutionSpace>
using TestSharedSpace = typename ExecutionSpace::memory_space;
#endif

using ExecSpace  = Kokkos::DefaultExecutionSpace;
using MemSpace   = TestSharedSpace<Kokkos::DefaultExecutionSpace>;
using DeviceType = Kokkos::Device<ExecSpace, MemSpace>;

using DualViewType = Kokkos::DualView<double*, Kokkos::LayoutLeft, DeviceType>;
using ConstDualViewType =
    Kokkos::DualView<const double*, Kokkos::LayoutLeft, DeviceType>;
using d_device = DeviceType;
using h_device =
    Kokkos::Device<Kokkos::DefaultHostExecutionSpace,
                   TestSharedSpace<Kokkos::DefaultHostExecutionSpace>>;

TEST(TEST_CATEGORY, dualview_device_correct_kokkos_device) {
  DualViewType dv("myView", 100);
  dv.clear_sync_state();
  auto v_d      = dv.template view<d_device>();
  using vdt     = decltype(v_d);
  using vdt_d   = vdt::device_type;
  using vdt_d_e = vdt_d::execution_space;
  ASSERT_STREQ(vdt_d_e::name(), Kokkos::DefaultExecutionSpace::name());
}
TEST(TEST_CATEGORY, dualview_host_correct_kokkos_device) {
  DualViewType dv("myView", 100);
  dv.clear_sync_state();
  auto v_h      = dv.template view<h_device>();
  using vht     = decltype(v_h);
  using vht_d   = vht::device_type;
  using vht_d_e = vht_d::execution_space;
  ASSERT_STREQ(vht_d_e::name(), Kokkos::DefaultHostExecutionSpace::name());
}

TEST(TEST_CATEGORY, dualview_host_modify_template_device_sync) {
  DualViewType dv("myView", 100);
  dv.clear_sync_state();
  dv.modify_host();
  dv.template sync<d_device>();
  EXPECT_TRUE(!dv.need_sync_device());
  EXPECT_TRUE(!dv.need_sync_host());
  dv.clear_sync_state();
}

TEST(TEST_CATEGORY, dualview_host_modify_template_device_execspace_sync) {
  DualViewType dv("myView", 100);
  dv.clear_sync_state();
  dv.modify_host();
  dv.template sync<d_device::execution_space>();
  EXPECT_TRUE(!dv.need_sync_device());
  EXPECT_TRUE(!dv.need_sync_host());
  dv.clear_sync_state();
}

TEST(TEST_CATEGORY, dualview_device_modify_template_host_sync) {
  DualViewType dv("myView", 100);
  dv.clear_sync_state();
  dv.modify_device();
  dv.template sync<h_device>();
  EXPECT_TRUE(!dv.need_sync_device());
  EXPECT_TRUE(!dv.need_sync_host());
  dv.clear_sync_state();
}
TEST(TEST_CATEGORY, dualview_device_modify_template_host_execspace_sync) {
  DualViewType dv("myView", 100);
  dv.clear_sync_state();
  dv.modify_device();
  dv.template sync<h_device::execution_space>();
  EXPECT_TRUE(!dv.need_sync_device());
  EXPECT_TRUE(!dv.need_sync_host());
  dv.clear_sync_state();
}

TEST(TEST_CATEGORY,
     dualview_template_views_return_correct_executionspace_views) {
  DualViewType dv("myView", 100);
  dv.clear_sync_state();
  using hvt = decltype(dv.view<Kokkos::DefaultHostExecutionSpace>());
  using dvt = decltype(dv.view<Kokkos::DefaultExecutionSpace>());
  ASSERT_STREQ(Kokkos::DefaultExecutionSpace::name(),
               dvt::device_type::execution_space::name());
  ASSERT_STREQ(Kokkos::DefaultHostExecutionSpace::name(),
               hvt::device_type::execution_space::name());
}

TEST(TEST_CATEGORY,
     dualview_template_views_return_correct_views_from_const_dual_view) {
  DualViewType dv("myView", 100);
  ConstDualViewType const_dv = dv;
  dv.clear_sync_state();
  ASSERT_EQ(dv.view<Kokkos::DefaultHostExecutionSpace>(),
            const_dv.view<Kokkos::DefaultHostExecutionSpace>());
  ASSERT_EQ(dv.view<Kokkos::DefaultExecutionSpace>(),
            const_dv.view<Kokkos::DefaultExecutionSpace>());
}

// User-defined types with a View data member, only host-constructible
template <class V>
class S {
  V v_;

 public:
  template <class... Extents>
  S(std::string label, Extents... extents) : v_(std::move(label), extents...) {}
  S() : v_("v", 10) {}
};

template <typename V>
auto initialize_view_of_views() {
  Kokkos::DualView<V*, TEST_EXECSPACE> dv_v(
      Kokkos::view_alloc("myView", Kokkos::SequentialHostInit), 3u);

  V v("v", 2);
  V w("w", 2);
  dv_v.view_host()(0) = v;
  dv_v.view_host()(1) = w;

  dv_v.modify_host();
  dv_v.sync_device();

  return dv_v;
}

TEST(TEST_CATEGORY, dualview_sequential_host_init) {
  auto dv_v = initialize_view_of_views<Kokkos::View<double*, TEST_EXECSPACE>>();
  dv_v.resize(Kokkos::view_alloc(Kokkos::SequentialHostInit), 2u);
  ASSERT_EQ(dv_v.view_device().size(), 2u);
  ASSERT_EQ(dv_v.view_host().size(), 2u);

  initialize_view_of_views<S<Kokkos::View<double*, TEST_EXECSPACE>>>();

  Kokkos::DualView<double*> dv(
      Kokkos::view_alloc("myView", Kokkos::SequentialHostInit), 1u);
  dv.resize(Kokkos::view_alloc(Kokkos::SequentialHostInit), 2u);
  ASSERT_EQ(dv.view_device().size(), 2u);
  ASSERT_EQ(dv.view_host().size(), 2u);
  dv.realloc(Kokkos::view_alloc(Kokkos::SequentialHostInit), 3u);
  ASSERT_EQ(dv.view_device().size(), 3u);
  ASSERT_EQ(dv.view_host().size(), 3u);
}

TEST(TEST_CATEGORY, dualview_default_constructed) {
  DualViewType dv;

  dv.modify<DualViewType::t_dev>();
  ASSERT_FALSE(dv.need_sync_host());
  ASSERT_FALSE(dv.need_sync_device());
  dv.sync<DualViewType::t_dev>();

  dv.modify_host();
  ASSERT_FALSE(dv.need_sync_host());
  ASSERT_FALSE(dv.need_sync_device());
  dv.sync_host();

  dv.modify_device();
  ASSERT_FALSE(dv.need_sync_host());
  ASSERT_FALSE(dv.need_sync_device());
  dv.sync_device();
}
}  // anonymous namespace
}  // namespace Test

#endif  // KOKKOS_TEST_DUALVIEW_HPP
