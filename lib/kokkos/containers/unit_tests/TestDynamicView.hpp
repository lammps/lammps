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

#ifndef KOKKOS_TEST_DYNAMICVIEW_HPP
#define KOKKOS_TEST_DYNAMICVIEW_HPP

#include <gtest/gtest.h>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <Kokkos_Core.hpp>

#include <Kokkos_DynamicView.hpp>
#include <Kokkos_Timer.hpp>

namespace Test {

template <typename Scalar, class Space>
struct TestDynamicView {
  using execution_space = typename Space::execution_space;
  using memory_space    = typename Space::memory_space;

  using view_type = Kokkos::Experimental::DynamicView<Scalar*, Space>;

  using value_type = double;

  static void run(unsigned arg_total_size) {
    // Test: Create DynamicView, initialize size (via resize), run through
    // parallel_for to set values, check values (via parallel_reduce); resize
    // values and repeat
    //   Case 1: min_chunk_size is a power of 2
    {
      {
        view_type d1;
        ASSERT_FALSE(d1.is_allocated());

        d1 = view_type("d1", 1024, arg_total_size);
        view_type d2(d1);
        view_type d3("d3", 1024, arg_total_size);

        ASSERT_FALSE(d1.is_allocated());
        ASSERT_FALSE(d2.is_allocated());
        ASSERT_FALSE(d3.is_allocated());

        unsigned d_size = arg_total_size / 8;
        d1.resize_serial(d_size);
        d2.resize_serial(d_size);
        d3.resize_serial(d_size);

        ASSERT_TRUE(d1.is_allocated());
        ASSERT_TRUE(d2.is_allocated());
        ASSERT_TRUE(d3.is_allocated());
      }
      view_type da("da", 1024, arg_total_size);
      ASSERT_EQ(da.size(), 0u);
      // Init
      unsigned da_size = arg_total_size / 8;
      da.resize_serial(da_size);
      ASSERT_EQ(da.size(), da_size);

#if defined(KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA)
      Kokkos::parallel_for(
          Kokkos::RangePolicy<execution_space>(0, da_size),
          KOKKOS_LAMBDA(const int i) { da(i) = Scalar(i); });

      value_type result_sum = 0.0;
      Kokkos::parallel_reduce(
          Kokkos::RangePolicy<execution_space>(0, da_size),
          KOKKOS_LAMBDA(const int i, value_type& partial_sum) {
            partial_sum += (value_type)da(i);
          },
          result_sum);

      ASSERT_EQ(result_sum, (value_type)(da_size * (da_size - 1) / 2));
#endif

      // add 3x more entries i.e. 4x larger than previous size
      // the first 1/4 should remain the same
      unsigned da_resize = arg_total_size / 2;
      da.resize_serial(da_resize);
      ASSERT_EQ(da.size(), da_resize);

#if defined(KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA)
      Kokkos::parallel_for(
          Kokkos::RangePolicy<execution_space>(da_size, da_resize),
          KOKKOS_LAMBDA(const int i) { da(i) = Scalar(i); });

      value_type new_result_sum = 0.0;
      Kokkos::parallel_reduce(
          Kokkos::RangePolicy<execution_space>(da_size, da_resize),
          KOKKOS_LAMBDA(const int i, value_type& partial_sum) {
            partial_sum += (value_type)da(i);
          },
          new_result_sum);

      ASSERT_EQ(new_result_sum + result_sum,
                (value_type)(da_resize * (da_resize - 1) / 2));
#endif
    }  // end scope

    // Test: Create DynamicView, initialize size (via resize), run through
    // parallel_for to set values, check values (via parallel_reduce); resize
    // values and repeat
    //   Case 2: min_chunk_size is NOT a power of 2
    {
      view_type da("da", 1023, arg_total_size);
      ASSERT_EQ(da.size(), 0u);
      // Init
      unsigned da_size = arg_total_size / 8;
      da.resize_serial(da_size);
      ASSERT_EQ(da.size(), da_size);

#if defined(KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA)
      Kokkos::parallel_for(
          Kokkos::RangePolicy<execution_space>(0, da_size),
          KOKKOS_LAMBDA(const int i) { da(i) = Scalar(i); });

      value_type result_sum = 0.0;
      Kokkos::parallel_reduce(
          Kokkos::RangePolicy<execution_space>(0, da_size),
          KOKKOS_LAMBDA(const int i, value_type& partial_sum) {
            partial_sum += (value_type)da(i);
          },
          result_sum);

      ASSERT_EQ(result_sum, (value_type)(da_size * (da_size - 1) / 2));
#endif

      // add 3x more entries i.e. 4x larger than previous size
      // the first 1/4 should remain the same
      unsigned da_resize = arg_total_size / 2;
      da.resize_serial(da_resize);
      ASSERT_EQ(da.size(), da_resize);

#if defined(KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA)
      Kokkos::parallel_for(
          Kokkos::RangePolicy<execution_space>(da_size, da_resize),
          KOKKOS_LAMBDA(const int i) { da(i) = Scalar(i); });

      value_type new_result_sum = 0.0;
      Kokkos::parallel_reduce(
          Kokkos::RangePolicy<execution_space>(da_size, da_resize),
          KOKKOS_LAMBDA(const int i, value_type& partial_sum) {
            partial_sum += (value_type)da(i);
          },
          new_result_sum);

      ASSERT_EQ(new_result_sum + result_sum,
                (value_type)(da_resize * (da_resize - 1) / 2));
#endif
    }  // end scope

    // Test: Create DynamicView, initialize size (via resize), run through
    // parallel_for to set values, check values (via parallel_reduce); resize
    // values and repeat
    //   Case 3: resize reduces the size
    {
      view_type da("da", 1023, arg_total_size);
      ASSERT_EQ(da.size(), 0u);
      // Init
      unsigned da_size = arg_total_size / 2;
      da.resize_serial(da_size);
      ASSERT_EQ(da.size(), da_size);

#if defined(KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA)
      Kokkos::parallel_for(
          Kokkos::RangePolicy<execution_space>(0, da_size),
          KOKKOS_LAMBDA(const int i) { da(i) = Scalar(i); });

      value_type result_sum = 0.0;
      Kokkos::parallel_reduce(
          Kokkos::RangePolicy<execution_space>(0, da_size),
          KOKKOS_LAMBDA(const int i, value_type& partial_sum) {
            partial_sum += (value_type)da(i);
          },
          result_sum);

      ASSERT_EQ(result_sum, (value_type)(da_size * (da_size - 1) / 2));
#endif

      // remove the final 3/4 entries i.e. first 1/4 remain
      unsigned da_resize = arg_total_size / 8;
      da.resize_serial(da_resize);
      ASSERT_EQ(da.size(), da_resize);

#if defined(KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA)
      Kokkos::parallel_for(
          Kokkos::RangePolicy<execution_space>(0, da_resize),
          KOKKOS_LAMBDA(const int i) { da(i) = Scalar(i); });

      value_type new_result_sum = 0.0;
      Kokkos::parallel_reduce(
          Kokkos::RangePolicy<execution_space>(0, da_resize),
          KOKKOS_LAMBDA(const int i, value_type& partial_sum) {
            partial_sum += (value_type)da(i);
          },
          new_result_sum);

      ASSERT_EQ(new_result_sum, (value_type)(da_resize * (da_resize - 1) / 2));
#endif
    }  // end scope

    // Test: Reproducer to demonstrate compile-time error of deep_copy
    // of DynamicView to/from on-host View.
    //   Case 4:
    {
      using device_view_type = Kokkos::View<Scalar*, Space>;
      using host_view_type = typename Kokkos::View<Scalar*, Space>::HostMirror;

      view_type device_dynamic_view("on-device DynamicView", 1024,
                                    arg_total_size);
      device_view_type device_view("on-device View", arg_total_size);
      host_view_type host_view("on-host View", arg_total_size);

      unsigned da_size = arg_total_size / 8;
      device_dynamic_view.resize_serial(da_size);

      // Use parallel_for to populate device_dynamic_view and verify values
#if defined(KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA)
      Kokkos::parallel_for(
          Kokkos::RangePolicy<execution_space>(0, da_size),
          KOKKOS_LAMBDA(const int i) { device_dynamic_view(i) = Scalar(i); });

      value_type result_sum = 0.0;
      Kokkos::parallel_reduce(
          Kokkos::RangePolicy<execution_space>(0, da_size),
          KOKKOS_LAMBDA(const int i, value_type& partial_sum) {
            partial_sum += (value_type)device_dynamic_view(i);
          },
          result_sum);

      ASSERT_EQ(result_sum, (value_type)(da_size * (da_size - 1) / 2));
#endif

      // Use an on-device View as intermediate to deep_copy the
      // device_dynamic_view to host, zero out the device_dynamic_view,
      // deep_copy from host back to the device_dynamic_view and verify
      Kokkos::deep_copy(device_view, device_dynamic_view);
      Kokkos::deep_copy(host_view, device_view);
      Kokkos::deep_copy(device_view, host_view);
#if defined(KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA)
      Kokkos::parallel_for(
          Kokkos::RangePolicy<execution_space>(0, da_size),
          KOKKOS_LAMBDA(const int i) { device_dynamic_view(i) = Scalar(0); });
#endif
      Kokkos::deep_copy(device_dynamic_view, device_view);
#if defined(KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA)
      value_type new_result_sum = 0.0;
      Kokkos::parallel_reduce(
          Kokkos::RangePolicy<execution_space>(0, da_size),
          KOKKOS_LAMBDA(const int i, value_type& partial_sum) {
            partial_sum += (value_type)device_dynamic_view(i);
          },
          new_result_sum);

      ASSERT_EQ(new_result_sum, (value_type)(da_size * (da_size - 1) / 2));
#endif

      // Try to deep_copy device_dynamic_view directly to/from host.
      // host-to-device currently fails to compile because DP and SP are
      // swapped in the deep_copy implementation.
      // Once that's fixed, both deep_copy's will fail at runtime because the
      // destination execution space cannot access the source memory space.
      // Check if the memory spaces are different before testing the deep_copy.
      if (!Kokkos::SpaceAccessibility<Kokkos::HostSpace,
                                      memory_space>::accessible) {
        ASSERT_THROW(Kokkos::deep_copy(host_view, device_dynamic_view),
                     std::runtime_error);
        ASSERT_THROW(Kokkos::deep_copy(device_dynamic_view, host_view),
                     std::runtime_error);
      }
    }
  }
};

TEST(TEST_CATEGORY, dynamic_view) {
  using TestDynView = TestDynamicView<double, TEST_EXECSPACE>;

  for (int i = 0; i < 10; ++i) {
    TestDynView::run(100000 + 100 * i);
  }
}

}  // namespace Test

#endif /* #ifndef KOKKOS_TEST_DYNAMICVIEW_HPP */
