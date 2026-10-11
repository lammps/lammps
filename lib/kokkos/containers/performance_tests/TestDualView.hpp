// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_TEST_DUALVIEW_HPP
#define KOKKOS_TEST_DUALVIEW_HPP

#include <iostream>
#include <cstdlib>
#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
import kokkos.dual_view;
#else
#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#include <Kokkos_Abort.hpp>
#endif
#include <Kokkos_Timer.hpp>

namespace Performance {

namespace Impl {

struct times_t {
  double t_sync_to_host;
  double t_compute;
  double t_sync_to_device;
};

template <typename Scalar, class ViewType>
struct SumViewEntriesFunctor {
  using value_type = Scalar;
  ViewType view;
  SumViewEntriesFunctor(const ViewType& view_) : view(view_) {}

  template <typename I, typename J>
  KOKKOS_INLINE_FUNCTION void operator()(const I i, const J j,
                                         value_type& sum) const {
    sum += view(i, j);
  }
};

template <typename Scalar, class ViewType>
struct IncrViewEntriesFunctor {
  using value_type = Scalar;
  ViewType view;
  IncrViewEntriesFunctor(const ViewType& view_) : view(view_) {}

  template <typename I, typename J>
  KOKKOS_INLINE_FUNCTION void operator()(const I i, const J j) const {
    view(i, j)++;
  }
};

template <typename Scalar, class Device>
struct test_dualview_with_datacheck {
  using scalar_type     = Scalar;
  using execution_space = Device;

  template <typename ViewType>
  void run(const int n, const int m, times_t& times) {
    ViewType a, b;
    a = ViewType("A", n, m);
    b = ViewType("B", n, m);

    const scalar_type sum_total =
        static_cast<scalar_type>(n) * static_cast<scalar_type>(m);

    Kokkos::deep_copy(a.view_device(), 1);

    Kokkos::Timer timer;

    using device_space = typename ViewType::t_dev::execution_space;
    using host_space   = typename ViewType::t_host::execution_space;

    timer.reset();
    a.template modify<device_space>();
    a.template sync<host_space>();
    times.t_sync_to_host = timer.seconds();

    // Check device view is initialized as expected
    scalar_type a_d_sum = 0;

    timer.reset();
    Kokkos::parallel_reduce(
        Kokkos::MDRangePolicy<device_space, Kokkos::Rank<2>>({{0, 0}},
                                                             {{n, m}}),
        SumViewEntriesFunctor<scalar_type, typename ViewType::t_dev>(
            a.view_device()),
        a_d_sum);
    times.t_compute = timer.seconds();
    if (a_d_sum != sum_total)
      Kokkos::abort("test_dualview_with_datacheck: a_d_sum != sum_total");

    // Use deep_copy
    Kokkos::deep_copy(b, a);
    timer.reset();
    b.template sync<host_space>();
    times.t_sync_to_host += timer.seconds();

    // Perform same checks on b as done on a
    // Check device view is initialized as expected
    scalar_type b_d_sum = 0;
    timer.reset();
    Kokkos::parallel_reduce(
        Kokkos::MDRangePolicy<device_space, Kokkos::Rank<2>>({{0, 0}},
                                                             {{n, m}}),
        SumViewEntriesFunctor<scalar_type, typename ViewType::t_dev>(
            b.view_device()),
        b_d_sum);
    times.t_compute += timer.seconds();
    if (b_d_sum != sum_total)
      Kokkos::abort("test_dualview_with_datacheck: b_d_sum != sum_total");
  }
  test_dualview_with_datacheck(const int n, const int m) {
    times_t elapsed_time = {0.0, 0.0, 0.0};
    run<Kokkos::DualView<Scalar**, Kokkos::LayoutLeft, Device>>(n, m,
                                                                elapsed_time);
    std::cout << " DualView (test_dualview_with_datacheck) timing (sec): "
              << "t_sync_to_host: " << elapsed_time.t_sync_to_host
              << ", t_compute: " << elapsed_time.t_compute
              << ", t_sync_to_device: " << elapsed_time.t_sync_to_device
              << ", dim0: " << n << ", dim1: " << m << std::endl;
  }
};

template <typename Scalar, class Device>
struct test_dualview_sync {
  using scalar_type     = Scalar;
  using execution_space = Device;

  template <typename ViewType>
  void run(const int n, const int m, const int iters, times_t& elapsed_time) {
    ViewType a = ViewType("A", n, m);

    using device_space = typename ViewType::t_dev::execution_space;
    using host_space   = typename ViewType::t_host::execution_space;

    Kokkos::deep_copy(a.view_device(), 1);

    int i = iters;

    Kokkos::Timer timer;
    while (i--) {
      // Sync to host
      timer.reset();
      a.template modify<device_space>();
      a.template sync<host_space>();
      elapsed_time.t_sync_to_host += timer.seconds();

      // Update on host
      timer.reset();
      Kokkos::parallel_for(
          Kokkos::MDRangePolicy<
              host_space,
              Kokkos::Rank<2, Kokkos::Iterate::Left, Kokkos::Iterate::Left>>(
              {{0, 0}}, {{n, m}}),
          IncrViewEntriesFunctor<scalar_type, typename ViewType::t_host>(
              a.view_host()));
      Kokkos::fence();
      elapsed_time.t_compute += timer.seconds();

      // Sync to device
      timer.reset();
      a.template modify<host_space>();
      a.template sync<device_space>();
      elapsed_time.t_sync_to_device += timer.seconds();

      // Update on device
      timer.reset();
      Kokkos::parallel_for(
          Kokkos::MDRangePolicy<
              device_space,
              Kokkos::Rank<2, Kokkos::Iterate::Left, Kokkos::Iterate::Left>>(
              {{0, 0}}, {{n, m}}),
          IncrViewEntriesFunctor<scalar_type, typename ViewType::t_dev>(
              a.view_device()));
      Kokkos::fence();
      elapsed_time.t_compute += timer.seconds();
    }
    Kokkos::fence();

    const scalar_type sum_total =
        static_cast<scalar_type>(n) * static_cast<scalar_type>(m);
    scalar_type a_d_sum = 0;
    Kokkos::parallel_reduce(
        Kokkos::MDRangePolicy<
            device_space,
            Kokkos::Rank<2, Kokkos::Iterate::Left, Kokkos::Iterate::Left>>(
            {{0, 0}}, {{n, m}}),
        SumViewEntriesFunctor<scalar_type, typename ViewType::t_dev>(
            a.view_device()),
        a_d_sum);
    if (a_d_sum != sum_total + iters * 2 * sum_total)
      Kokkos::abort(
          "test_dualview_sync: a_d_sum != sum_total + iters * 2 * sum_total");
  }
  test_dualview_sync(const int n, const int m) {
    times_t elapsed_time = {0.0, 0.0, 0.0};
    const int iters      = 10;
    run<Kokkos::DualView<Scalar**, Kokkos::LayoutLeft, Device>>(n, m, iters,
                                                                elapsed_time);
    std::cout << " DualView (test_dualview_sync) timing (sec, avg per iter): "
              << "t_sync_to_host: "
              << elapsed_time.t_sync_to_host / static_cast<double>(iters)
              << ", t_compute: "
              << elapsed_time.t_compute / static_cast<double>(iters)
              << ", t_sync_to_device: "
              << elapsed_time.t_sync_to_device / static_cast<double>(iters)
              << ", dim0: " << n << ", dim1: " << m << std::endl;
  }
};
}  // namespace Impl

template <typename Scalar, typename Device>
void test_dualview() {
  Impl::test_dualview_with_datacheck<Scalar, Device>(128, 128);
  Impl::test_dualview_with_datacheck<Scalar, Device>(512, 512);
  Impl::test_dualview_with_datacheck<Scalar, Device>(1024, 1024);
  Impl::test_dualview_sync<Scalar, Device>(128, 128);
  Impl::test_dualview_sync<Scalar, Device>(512, 512);
  Impl::test_dualview_sync<Scalar, Device>(1024, 1024);
}
}  // namespace Performance

#endif  // KOKKOS_TEST_DUALVIEW_HPP
