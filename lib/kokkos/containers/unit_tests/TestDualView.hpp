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
    result = run_me<Kokkos::DualView<Scalar**, Kokkos::LayoutLeft, Device> >(
        size, 3);
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
        Kokkos::DualView<const Scalar * [m], Kokkos::LayoutLeft, Device>;

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
    if (!std::decay_t<SrcViewType>::impl_dualview_is_single_device::value) {
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
    Kokkos::deep_copy(a.d_view, 1);

    a.template modify<typename ViewType::execution_space>();
    a.template sync<typename ViewType::host_mirror_space>();
    a.template sync<typename ViewType::host_mirror_space>(
        Kokkos::DefaultExecutionSpace{});

    a.h_view(5, 1) = 3;
    a.h_view(6, 1) = 4;
    a.h_view(7, 2) = 5;
    a.template modify<typename ViewType::host_mirror_space>();
    ViewType b = Kokkos::subview(a, std::pair<unsigned int, unsigned int>(6, 9),
                                 std::pair<unsigned int, unsigned int>(0, 1));
    a.template sync<typename ViewType::execution_space>();
    a.template sync<typename ViewType::execution_space>(
        Kokkos::DefaultExecutionSpace{});
    b.template modify<typename ViewType::execution_space>();

    Kokkos::deep_copy(b.d_view, 2);

    a.template sync<typename ViewType::host_mirror_space>();
    a.template sync<typename ViewType::host_mirror_space>(
        Kokkos::DefaultExecutionSpace{});
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
      a.sync_host(Kokkos::DefaultExecutionSpace{});
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
      b.sync_host(Kokkos::DefaultExecutionSpace{});
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

template <typename Scalar, class Device, bool Initialize>
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
    if (Initialize)
      Kokkos::resize(Kokkos::WithoutInitializing, a, factor * n, factor * m);
    else
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
    ASSERT_EQ(a_h_sum, scalar_type(a.extent(0) * a.extent(1)));

    /* Covers case "Resize on Host" */
    a.modify_host();

    if (Initialize)
      Kokkos::resize(Kokkos::WithoutInitializing, a, n / factor, m / factor);
    else
      Kokkos::resize(a, n / factor, m / factor);
    ASSERT_EQ(a.extent(0), n / factor);
    ASSERT_EQ(a.extent(1), m / factor);

    a.sync_device();
    a.sync_device(Kokkos::DefaultExecutionSpace{});

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
    ASSERT_EQ(a_h_sum, scalar_type(a.extent(0) * a.extent(1)));
    ASSERT_EQ(a_h_sum, a_d_sum);

  }  // end run_me

  test_dualview_resize() {
    run_me<Kokkos::DualView<Scalar**, Kokkos::LayoutLeft, Device> >();
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

    ViewType a("A", n, m);
    if (Initialize)
      Kokkos::realloc(Kokkos::WithoutInitializing, a, n, m);
    else
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
    ASSERT_EQ(a_h_sum, scalar_type(a.extent(0) * a.extent(1)));
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

TEST(TEST_CATEGORY, dualview_realloc) {
  test_dualview_realloc<int, TEST_EXECSPACE>();
}

TEST(TEST_CATEGORY, dualview_resize) {
  test_dualview_resize<int, TEST_EXECSPACE>();
}

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
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_SYCL) || \
    defined(KOKKOS_ENABLE_HIP)  // OR other UVM builds
#define UVM_ENABLED_BUILD
#endif

#ifdef UVM_ENABLED_BUILD
template <typename ExecSpace>
struct UVMSpaceFor;
#endif

#ifdef KOKKOS_ENABLE_CUDA  // specific to CUDA
template <>
struct UVMSpaceFor<Kokkos::Cuda> {
  using type = Kokkos::CudaUVMSpace;
};
#endif

#ifdef KOKKOS_ENABLE_SYCL  // specific to SYCL
template <>
struct UVMSpaceFor<Kokkos::Experimental::SYCL> {
  using type = Kokkos::Experimental::SYCLSharedUSMSpace;
};
#endif

#ifdef KOKKOS_ENABLE_HIP  // specific to HIP
template <>
struct UVMSpaceFor<Kokkos::HIP> {
  using type = Kokkos::HIPManagedSpace;
};
#endif

#ifdef UVM_ENABLED_BUILD
template <>
struct UVMSpaceFor<Kokkos::DefaultHostExecutionSpace> {
  using type = typename UVMSpaceFor<Kokkos::DefaultExecutionSpace>::type;
};
#else
template <typename ExecSpace>
struct UVMSpaceFor {
  using type = typename ExecSpace::memory_space;
};
#endif

using ExecSpace  = Kokkos::DefaultExecutionSpace;
using MemSpace   = typename UVMSpaceFor<Kokkos::DefaultExecutionSpace>::type;
using DeviceType = Kokkos::Device<ExecSpace, MemSpace>;

using DualViewType = Kokkos::DualView<double*, Kokkos::LayoutLeft, DeviceType>;
using d_device     = DeviceType;
using h_device     = Kokkos::Device<
    Kokkos::DefaultHostExecutionSpace,
    typename UVMSpaceFor<Kokkos::DefaultHostExecutionSpace>::type>;

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
  using hvt = decltype(dv.view<typename Kokkos::DefaultHostExecutionSpace>());
  using dvt = decltype(dv.view<typename Kokkos::DefaultExecutionSpace>());
  ASSERT_STREQ(Kokkos::DefaultExecutionSpace::name(),
               dvt::device_type::execution_space::name());
  ASSERT_STREQ(Kokkos::DefaultHostExecutionSpace::name(),
               hvt::device_type::execution_space::name());
}

}  // anonymous namespace
}  // namespace Test

#endif  // KOKKOS_TEST_DUALVIEW_HPP
