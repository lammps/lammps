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

#ifndef KOKKOS_TEST_DYNAMICVIEW_HPP
#define KOKKOS_TEST_DYNAMICVIEW_HPP

#include <gtest/gtest.h>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <Kokkos_Core.hpp>

#include <Kokkos_DynamicView.hpp>
#include <impl/Kokkos_Timer.hpp>

namespace Test {

template <typename Scalar, class Space>
struct TestDynamicView {
  typedef typename Space::execution_space execution_space;
  typedef typename Space::memory_space memory_space;

  typedef Kokkos::Experimental::DynamicView<Scalar*, Space> view_type;

  typedef double value_type;

  static void run(unsigned arg_total_size) {
    // Test: Create DynamicView, initialize size (via resize), run through
    // parallel_for to set values, check values (via parallel_reduce); resize
    // values and repeat
    //   Case 1: min_chunk_size is a power of 2
    {
      view_type da("da", 1024, arg_total_size);
      ASSERT_EQ(da.size(), 0);
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
      ASSERT_EQ(da.size(), 0);
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
      ASSERT_EQ(da.size(), 0);
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
  }
};

TEST(TEST_CATEGORY, dynamic_view) {
  typedef TestDynamicView<double, TEST_EXECSPACE> TestDynView;

  for (int i = 0; i < 10; ++i) {
    TestDynView::run(100000 + 100 * i);
  }
}

}  // namespace Test

#endif /* #ifndef KOKKOS_TEST_DYNAMICVIEW_HPP */
