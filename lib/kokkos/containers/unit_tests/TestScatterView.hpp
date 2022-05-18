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

#ifndef KOKKOS_TEST_SCATTER_VIEW_HPP
#define KOKKOS_TEST_SCATTER_VIEW_HPP

#include <Kokkos_ScatterView.hpp>
#include <gtest/gtest.h>

namespace Test {

template <typename DeviceType, typename Layout, typename Duplication,
          typename Contribution, typename Op, typename NumberType>
struct test_scatter_view_impl_cls;

template <typename DeviceType, typename Layout, typename Duplication,
          typename Contribution, typename NumberType>
struct test_scatter_view_impl_cls<DeviceType, Layout, Duplication, Contribution,
                                  Kokkos::Experimental::ScatterSum,
                                  NumberType> {
 public:
  using scatter_view_type =
      Kokkos::Experimental::ScatterView<NumberType * [12], Layout, DeviceType,
                                        Kokkos::Experimental::ScatterSum,
                                        Duplication, Contribution>;

  using orig_view_type = Kokkos::View<NumberType * [12], Layout, DeviceType>;

  using size_type = typename Kokkos::HostSpace::size_type;

  scatter_view_type scatter_view;
  int scatterSize;

  test_scatter_view_impl_cls(const scatter_view_type& view) {
    scatter_view = view;
    scatterSize  = 0;
  }

  void initialize(orig_view_type orig) {
    auto host_view =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), orig);
    Kokkos::fence();
    Kokkos::deep_copy(host_view, 0);
    Kokkos::fence();
    Kokkos::deep_copy(orig, host_view);
  }

  void run_parallel(int n) {
    scatterSize = n;
    auto policy =
        Kokkos::RangePolicy<typename DeviceType::execution_space, int>(0, n);
    Kokkos::parallel_for(policy, *this, "scatter_view_test: Sum");
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(int i) const {
    auto scatter_access = scatter_view.access();
    auto scatter_access_atomic =
        scatter_view.template access<Kokkos::Experimental::ScatterAtomic>();
    for (int j = 0; j < 10; ++j) {
      auto k = (i + j) % scatterSize;
      scatter_access(k, 0) += 4;
      ++scatter_access(k, 1);
      --scatter_access(k, 2);
      scatter_access(k, 3)++;
      scatter_access(k, 4)--;
      scatter_access(k, 5) -= 5;
// Workaround Intel 17 compiler bug which sometimes add random
// instruction alignment which makes the lock instruction
// illegal. Seems to be mostly just for unsigned int atomics.
// Looking at the assembly the compiler
// appears to insert cache line alignment for the instruction.
// Isn't restricted to specific archs. Seen it on SNB and SKX, but for
// different code. Another occurrence was with Desul atomics in
// a different unit test. This one here happens without desul atomics.
// Inserting an assembly nop instruction changes the alignment and
// works round this.
#ifdef KOKKOS_COMPILER_INTEL
#if (KOKKOS_COMPILER_INTEL < 1800)
      asm volatile("nop\n");
#endif
#endif
      scatter_access_atomic(k, 6) += 2;
#ifdef KOKKOS_COMPILER_INTEL
#if (KOKKOS_COMPILER_INTEL < 1800)
      asm volatile("nop\n");
#endif
#endif
      scatter_access_atomic(k, 7)++;
#ifdef KOKKOS_COMPILER_INTEL
#if (KOKKOS_COMPILER_INTEL < 1800)
      asm volatile("nop\n");
#endif
#endif
      scatter_access_atomic(k, 8)--;
#ifdef KOKKOS_COMPILER_INTEL
#if (KOKKOS_COMPILER_INTEL < 1800)
      asm volatile("nop\n");
#endif
#endif
      --scatter_access_atomic(k, 9);
#ifdef KOKKOS_COMPILER_INTEL
#if (KOKKOS_COMPILER_INTEL < 1800)
      asm volatile("nop\n");
#endif
#endif
      ++scatter_access_atomic(k, 10);
#ifdef KOKKOS_COMPILER_INTEL
#if (KOKKOS_COMPILER_INTEL < 1800)
      asm volatile("nop\n");
#endif
#endif
      scatter_access(k, 11) -= 3;
    }
  }

  void validateResults(orig_view_type orig) {
    auto host_view =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), orig);
    Kokkos::fence();
    for (size_type i = 0; i < host_view.extent(0); ++i) {
      for (size_type j = 0; j < host_view.extent(1); ++j) {
        EXPECT_NEAR(host_view(i, j), NumberType(ref[j]), 1e-14)
            << "Data differs at indices " << i << ", " << j;
      }
    }
  }

  // check for correct padding
  void validateResultsForSubview(
      orig_view_type orig, std::pair<size_type, size_type>& subRangeDim0,
      std::pair<size_type, size_type>& subRangeDim1) {
    auto host_view =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), orig);
    Kokkos::fence();
    for (size_type i = 0; i < host_view.extent(0); ++i) {
      for (size_type j = 0; j < host_view.extent(1); ++j) {
        auto val = host_view(i, j);
        if ((i >= std::get<0>(subRangeDim0) && i < std::get<1>(subRangeDim0)) &&
            (j >= std::get<0>(subRangeDim1) && j < std::get<1>(subRangeDim1))) {
          // is in subview
          EXPECT_NEAR(val, NumberType(ref[j]), 1e-14)
              << "Data differs at indices " << i << ", " << j;
        } else {
          // is outside of subview
          EXPECT_NEAR(val, NumberType(0), 1e-14)
              << "Data differs at indices " << i << ", " << j;
        }
      }
    }
  }

 private:
  NumberType ref[12] = {80, 20, -20, 20, -20, -100, 40, 20, -20, -20, 20, -60};
};

template <typename DeviceType, typename Layout, typename Duplication,
          typename Contribution, typename NumberType>
struct test_scatter_view_impl_cls<DeviceType, Layout, Duplication, Contribution,
                                  Kokkos::Experimental::ScatterProd,
                                  NumberType> {
 public:
  using scatter_view_type =
      Kokkos::Experimental::ScatterView<NumberType * [3], Layout, DeviceType,
                                        Kokkos::Experimental::ScatterProd,
                                        Duplication, Contribution>;

  using orig_view_type = Kokkos::View<NumberType * [3], Layout, DeviceType>;

  using size_type = typename Kokkos::HostSpace::size_type;

  scatter_view_type scatter_view;
  int scatterSize;

  test_scatter_view_impl_cls(const scatter_view_type& view) {
    scatter_view = view;
    scatterSize  = 0;
  }

  void initialize(orig_view_type orig) {
    auto host_view =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), orig);
    Kokkos::fence();
    for (size_type i = 0; i < host_view.extent(0); ++i) {
      host_view(i, 0) = 1.0;
      host_view(i, 1) = 1.0;
      host_view(i, 2) = 1.0;
    }
    Kokkos::fence();
    Kokkos::deep_copy(orig, host_view);
  }

  void run_parallel(int n) {
    scatterSize = n;
    auto policy =
        Kokkos::RangePolicy<typename DeviceType::execution_space, int>(0, n);
    Kokkos::parallel_for(policy, *this, "scatter_view_test: Prod");
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(int i) const {
    auto scatter_access = scatter_view.access();
    auto scatter_access_atomic =
        scatter_view.template access<Kokkos::Experimental::ScatterAtomic>();
    for (int j = 0; j < 4; ++j) {
      auto k = (i + j) % scatterSize;
      scatter_access(k, 0) *= 4.0;
      scatter_access_atomic(k, 1) *= 2.0;
      scatter_access(k, 2) *= 1.0;
    }
  }

  void validateResults(orig_view_type orig) {
    auto host_view =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), orig);
    Kokkos::fence();
    for (size_type i = 0; i < host_view.extent(0); ++i) {
      auto val0 = host_view(i, 0);
      auto val1 = host_view(i, 1);
      auto val2 = host_view(i, 2);
      EXPECT_TRUE(std::fabs((val0 - 65536.0) / 65536.0) < 1e-14)
          << "Data differs at index " << i;
      EXPECT_TRUE(std::fabs((val1 - 256.0) / 256.0) < 1e-14)
          << "Data differs at index " << i;
      EXPECT_TRUE(std::fabs((val2 - 1.0) / 1.0) < 1e-14)
          << "Data differs at index " << i;
    }
  }

  // check for correct padding
  void validateResultsForSubview(
      orig_view_type orig, std::pair<size_type, size_type>& subRangeDim0,
      std::pair<size_type, size_type>& subRangeDim1) {
    (void)subRangeDim1;
    auto host_view =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), orig);
    Kokkos::fence();
    for (size_type i = 0; i < host_view.extent(0); ++i) {
      auto val0 = host_view(i, 0);
      auto val1 = host_view(i, 1);
      auto val2 = host_view(i, 2);
      if (i >= std::get<0>(subRangeDim0) && i < std::get<1>(subRangeDim0)) {
        // is in subview
        EXPECT_TRUE(std::fabs((val0 - 65536.0) / 65536.0) < 1e-14);
        EXPECT_TRUE(std::fabs((val1 - 256.0) / 256.0) < 1e-14);
        EXPECT_TRUE(std::fabs((val2 - 1.0) / 1.0) < 1e-14);
      } else {
        // is outside of subview
        EXPECT_NEAR(val0, NumberType(1), 1e-14)
            << "Data differs at index " << i;
        EXPECT_NEAR(val1, NumberType(1), 1e-14)
            << "Data differs at index " << i;
        EXPECT_NEAR(val2, NumberType(1), 1e-14)
            << "Data differs at index " << i;
      }
    }
  }
};

template <typename DeviceType, typename Layout, typename Duplication,
          typename Contribution, typename NumberType>
struct test_scatter_view_impl_cls<DeviceType, Layout, Duplication, Contribution,
                                  Kokkos::Experimental::ScatterMin,
                                  NumberType> {
 public:
  using scatter_view_type =
      Kokkos::Experimental::ScatterView<NumberType * [3], Layout, DeviceType,
                                        Kokkos::Experimental::ScatterMin,
                                        Duplication, Contribution>;

  using orig_view_type = Kokkos::View<NumberType * [3], Layout, DeviceType>;

  using size_type = typename Kokkos::HostSpace::size_type;

  scatter_view_type scatter_view;
  int scatterSize;

  test_scatter_view_impl_cls(const scatter_view_type& view) {
    scatter_view = view;
    scatterSize  = 0;
  }

  void initialize(orig_view_type orig) {
    auto host_view =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), orig);
    Kokkos::fence();
    for (size_type i = 0; i < host_view.extent(0); ++i) {
      host_view(i, 0) = 999999.0;
      host_view(i, 1) = 999999.0;
      host_view(i, 2) = 999999.0;
    }
    Kokkos::fence();
    Kokkos::deep_copy(orig, host_view);
  }

  void run_parallel(int n) {
    scatterSize = n;
    auto policy =
        Kokkos::RangePolicy<typename DeviceType::execution_space, int>(0, n);
    Kokkos::parallel_for(policy, *this, "scatter_view_test: Prod");
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(int i) const {
    auto scatter_access = scatter_view.access();
    auto scatter_access_atomic =
        scatter_view.template access<Kokkos::Experimental::ScatterAtomic>();
    for (int j = 0; j < 4; ++j) {
      auto k = (i + j) % scatterSize;
      scatter_access(k, 0).update((NumberType)(j + 1) * 4);
      scatter_access_atomic(k, 1).update((NumberType)(j + 1) * 2.0);
      scatter_access(k, 2).update((NumberType)(j + 1) * 1.0);
    }
  }

  void validateResults(orig_view_type orig) {
    auto host_view =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), orig);
    Kokkos::fence();
    for (size_type i = 0; i < host_view.extent(0); ++i) {
      auto val0 = host_view(i, 0);
      auto val1 = host_view(i, 1);
      auto val2 = host_view(i, 2);
      EXPECT_TRUE(std::fabs((val0 - 4.0) / 4.0) < 1e-14)
          << "Data differs at index " << i;
      EXPECT_TRUE(std::fabs((val1 - 2.0) / 2.0) < 1e-14)
          << "Data differs at index " << i;
      EXPECT_TRUE(std::fabs((val2 - 1.0) / 1.0) < 1e-14)
          << "Data differs at index " << i;
    }
  }

  // check for correct padding
  void validateResultsForSubview(
      orig_view_type orig, std::pair<size_type, size_type>& subRangeDim0,
      std::pair<size_type, size_type>& subRangeDim1) {
    (void)subRangeDim1;
    auto host_view =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), orig);
    Kokkos::fence();
    for (size_type i = 0; i < host_view.extent(0); ++i) {
      auto val0 = host_view(i, 0);
      auto val1 = host_view(i, 1);
      auto val2 = host_view(i, 2);
      if (i >= std::get<0>(subRangeDim0) && i < std::get<1>(subRangeDim0)) {
        // is in subview
        EXPECT_TRUE(std::fabs((val0 - 4.0) / 4.0) < 1e-14)
            << "Data differs at index " << i;
        EXPECT_TRUE(std::fabs((val1 - 2.0) / 2.0) < 1e-14)
            << "Data differs at index " << i;
        EXPECT_TRUE(std::fabs((val2 - 1.0) / 1.0) < 1e-14)
            << "Data differs at index " << i;
      } else {
        // is outside of subview
        EXPECT_NEAR(val0, NumberType(999999), 1e-14)
            << "Data differs at index " << i;
        EXPECT_NEAR(val1, NumberType(999999), 1e-14)
            << "Data differs at index " << i;
        EXPECT_NEAR(val2, NumberType(999999), 1e-14)
            << "Data differs at index " << i;
      }
    }
  }
};

template <typename DeviceType, typename Layout, typename Duplication,
          typename Contribution, typename NumberType>
struct test_scatter_view_impl_cls<DeviceType, Layout, Duplication, Contribution,
                                  Kokkos::Experimental::ScatterMax,
                                  NumberType> {
 public:
  using scatter_view_type =
      Kokkos::Experimental::ScatterView<NumberType * [3], Layout, DeviceType,
                                        Kokkos::Experimental::ScatterMax,
                                        Duplication, Contribution>;

  using orig_view_type = Kokkos::View<NumberType * [3], Layout, DeviceType>;

  using size_type = typename Kokkos::HostSpace::size_type;

  scatter_view_type scatter_view;
  int scatterSize;

  test_scatter_view_impl_cls(const scatter_view_type& view) {
    scatter_view = view;
    scatterSize  = 0;
  }

  void initialize(orig_view_type orig) {
    auto host_view =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), orig);
    Kokkos::fence();
    for (size_type i = 0; i < host_view.extent(0); ++i) {
      host_view(i, 0) = 0.0;
      host_view(i, 1) = 0.0;
      host_view(i, 2) = 0.0;
    }
    Kokkos::fence();
    Kokkos::deep_copy(orig, host_view);
  }

  void run_parallel(int n) {
    scatterSize = n;
    Kokkos::RangePolicy<typename DeviceType::execution_space, int> policy(0, n);
    Kokkos::parallel_for(policy, *this, "scatter_view_test: Prod");
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(int i) const {
    auto scatter_access = scatter_view.access();
    auto scatter_access_atomic =
        scatter_view.template access<Kokkos::Experimental::ScatterAtomic>();
    for (int j = 0; j < 4; ++j) {
      auto k = (i + j) % scatterSize;
      scatter_access(k, 0).update((NumberType)(j + 1) * 4);
      scatter_access_atomic(k, 1).update((NumberType)(j + 1) * 2.0);
      scatter_access(k, 2).update((NumberType)(j + 1) * 1.0);
    }
  }

  void validateResults(orig_view_type orig) {
    auto host_view =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), orig);
    Kokkos::fence();
    for (size_type i = 0; i < host_view.extent(0); ++i) {
      auto val0 = host_view(i, 0);
      auto val1 = host_view(i, 1);
      auto val2 = host_view(i, 2);
      EXPECT_TRUE(std::fabs((val0 - 16.0) / 16.0) < 1e-14)
          << "Data differs at index " << i;
      EXPECT_TRUE(std::fabs((val1 - 8.0) / 8.0) < 1e-14)
          << "Data differs at index " << i;
      EXPECT_TRUE(std::fabs((val2 - 4.0) / 4.0) < 1e-14)
          << "Data differs at index " << i;
    }
  }

  // check for correct padding
  void validateResultsForSubview(
      orig_view_type orig, std::pair<size_type, size_type>& subRangeDim0,
      std::pair<size_type, size_type>& subRangeDim1) {
    (void)subRangeDim1;
    auto host_view =
        Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), orig);
    Kokkos::fence();
    for (size_type i = 0; i < host_view.extent(0); ++i) {
      auto val0 = host_view(i, 0);
      auto val1 = host_view(i, 1);
      auto val2 = host_view(i, 2);
      if (i >= std::get<0>(subRangeDim0) && i < std::get<1>(subRangeDim0)) {
        // is in subview
        EXPECT_TRUE(std::fabs((val0 - 16.0) / 16.0) < 1e-14)
            << "Data differs at index " << i;
        EXPECT_TRUE(std::fabs((val1 - 8.0) / 8.0) < 1e-14)
            << "Data differs at index " << i;
        EXPECT_TRUE(std::fabs((val2 - 4.0) / 4.0) < 1e-14)
            << "Data differs at index " << i;
      } else {
        // is outside of subview
        EXPECT_NEAR(val0, NumberType(0), 1e-14)
            << "Data differs at index " << i;
        EXPECT_NEAR(val1, NumberType(0), 1e-14)
            << "Data differs at index " << i;
        EXPECT_NEAR(val2, NumberType(0), 1e-14)
            << "Data differs at index " << i;
      }
    }
  }
};

// Test ScatterView on subview
template <typename DeviceType, typename Layout, typename Op,
          typename NumberType>
struct test_default_scatter_sub_view {
 public:
  using default_duplication = Kokkos::Impl::Experimental::DefaultDuplication<
      typename DeviceType::execution_space>;
  using Duplication  = typename default_duplication::type;
  using Contribution = typename Kokkos::Impl::Experimental::DefaultContribution<
      typename DeviceType::execution_space, Duplication>::type;
  using scatter_view_def =
      typename test_scatter_view_impl_cls<DeviceType, Layout, Duplication,
                                          Contribution, Op,
                                          NumberType>::scatter_view_type;
  using orig_view_def =
      typename test_scatter_view_impl_cls<DeviceType, Layout, Duplication,
                                          Contribution, Op,
                                          NumberType>::orig_view_type;

  using size_type = typename Kokkos::HostSpace::size_type;

  void run_test(int n) {
    // Test creation via create_scatter_view overload 1
    {
      orig_view_def original_view("original_view", n);

      auto rangeDim0 = std::pair<size_type, size_type>(0 + 1, n - 1);
      auto rangeDim1 =
          std::pair<size_type, size_type>(0, original_view.extent(1));

      auto original_sub_view =
          Kokkos::subview(original_view, rangeDim0, rangeDim1);

      scatter_view_def scatter_view =
          Kokkos::Experimental::create_scatter_view(Op{}, original_sub_view);

      test_scatter_view_impl_cls<DeviceType, Layout, Duplication, Contribution,
                                 Op, NumberType>
          scatter_view_test_impl(scatter_view);
      scatter_view_test_impl.initialize(original_view);
      scatter_view_test_impl.run_parallel(original_sub_view.extent(0));

      Kokkos::Experimental::contribute(original_sub_view, scatter_view);
      scatter_view.reset_except(original_sub_view);

      scatter_view_test_impl.run_parallel(original_sub_view.extent(0));

      Kokkos::Experimental::contribute(original_sub_view, scatter_view);
      Kokkos::fence();

      scatter_view_test_impl.validateResultsForSubview(original_view, rangeDim0,
                                                       rangeDim1);
    }
  }
};

template <typename DeviceType, typename Layout, typename Op,
          typename NumberType>
struct test_default_scatter_view {
 public:
  using default_duplication = Kokkos::Impl::Experimental::DefaultDuplication<
      typename DeviceType::execution_space>;
  using Duplication  = typename default_duplication::type;
  using Contribution = typename Kokkos::Impl::Experimental::DefaultContribution<
      typename DeviceType::execution_space, Duplication>::type;
  using scatter_view_def =
      typename test_scatter_view_impl_cls<DeviceType, Layout, Duplication,
                                          Contribution, Op,
                                          NumberType>::scatter_view_type;
  using orig_view_def =
      typename test_scatter_view_impl_cls<DeviceType, Layout, Duplication,
                                          Contribution, Op,
                                          NumberType>::orig_view_type;

  void run_test(int n) {
    // Test creation via create_scatter_view overload 1
    {
      orig_view_def original_view("original_view", n);
      scatter_view_def scatter_view =
          Kokkos::Experimental::create_scatter_view(Op{}, original_view);

      test_scatter_view_impl_cls<DeviceType, Layout, Duplication, Contribution,
                                 Op, NumberType>
          scatter_view_test_impl(scatter_view);
      scatter_view_test_impl.initialize(original_view);
      scatter_view_test_impl.run_parallel(n);

      Kokkos::Experimental::contribute(original_view, scatter_view);
      scatter_view.reset_except(original_view);

      scatter_view_test_impl.run_parallel(n);

      Kokkos::Experimental::contribute(original_view, scatter_view);
      Kokkos::fence();

      scatter_view_test_impl.validateResults(original_view);

      {
        scatter_view_def persistent_view("persistent", n);
        auto result_view = persistent_view.subview();
        contribute(result_view, persistent_view);
        Kokkos::fence();
      }
    }
  }
};

template <typename DeviceType, typename Layout, typename Duplication,
          typename Contribution, typename Op, typename NumberType>
struct test_scatter_view_config {
 public:
  using scatter_view_def =
      typename test_scatter_view_impl_cls<DeviceType, Layout, Duplication,
                                          Contribution, Op,
                                          NumberType>::scatter_view_type;
  using orig_view_def =
      typename test_scatter_view_impl_cls<DeviceType, Layout, Duplication,
                                          Contribution, Op,
                                          NumberType>::orig_view_type;

  void compile_constructor() {
    auto sv = scatter_view_def(Kokkos::view_alloc(DeviceType{}, "label"), 10);
  }

  void run_test(int n) {
    // test allocation
    {
      orig_view_def ov1("ov1", n);
      scatter_view_def sv1;

      ASSERT_FALSE(sv1.is_allocated());

      sv1 = Kokkos::Experimental::create_scatter_view<Op, Duplication,
                                                      Contribution>(ov1);

      scatter_view_def sv2(sv1);
      scatter_view_def sv3("sv3", n);

      ASSERT_TRUE(sv1.is_allocated());
      ASSERT_TRUE(sv2.is_allocated());
      ASSERT_TRUE(sv3.is_allocated());
    }

    // Test creation via create_scatter_view
    {
      orig_view_def original_view("original_view", n);
      scatter_view_def scatter_view = Kokkos::Experimental::create_scatter_view<
          Op, Duplication, Contribution>(original_view);

      test_scatter_view_impl_cls<DeviceType, Layout, Duplication, Contribution,
                                 Op, NumberType>
          scatter_view_test_impl(scatter_view);
      scatter_view_test_impl.initialize(original_view);
      scatter_view_test_impl.run_parallel(n);

      Kokkos::Experimental::contribute(original_view, scatter_view);
      scatter_view.reset_except(original_view);

      scatter_view_test_impl.run_parallel(n);

      Kokkos::Experimental::contribute(original_view, scatter_view);
      Kokkos::fence();

      scatter_view_test_impl.validateResults(original_view);

      {
        scatter_view_def persistent_view("persistent", n);
        auto result_view = persistent_view.subview();
        contribute(result_view, persistent_view);
        Kokkos::fence();
      }
    }
    // Test creation via create_scatter_view overload 2
    {
      orig_view_def original_view("original_view", n);
      scatter_view_def scatter_view = Kokkos::Experimental::create_scatter_view(
          Op{}, Duplication{}, Contribution{}, original_view);

      test_scatter_view_impl_cls<DeviceType, Layout, Duplication, Contribution,
                                 Op, NumberType>
          scatter_view_test_impl(scatter_view);
      scatter_view_test_impl.initialize(original_view);
      scatter_view_test_impl.run_parallel(n);

      Kokkos::Experimental::contribute(original_view, scatter_view);
      scatter_view.reset_except(original_view);

      scatter_view_test_impl.run_parallel(n);

      Kokkos::Experimental::contribute(original_view, scatter_view);
      Kokkos::fence();

      scatter_view_test_impl.validateResults(original_view);

      {
        scatter_view_def persistent_view("persistent", n);
        auto result_view = persistent_view.subview();
        contribute(result_view, persistent_view);
        Kokkos::fence();
      }
    }
    // Test creation via constructor
    {
      orig_view_def original_view("original_view", n);
      scatter_view_def scatter_view(original_view);

      test_scatter_view_impl_cls<DeviceType, Layout, Duplication, Contribution,
                                 Op, NumberType>
          scatter_view_test_impl(scatter_view);
      scatter_view_test_impl.initialize(original_view);
      scatter_view_test_impl.run_parallel(n);

      Kokkos::Experimental::contribute(original_view, scatter_view);
      scatter_view.reset_except(original_view);

      scatter_view_test_impl.run_parallel(n);

      Kokkos::Experimental::contribute(original_view, scatter_view);
      Kokkos::fence();

      scatter_view_test_impl.validateResults(original_view);

      {
        scatter_view_def persistent_view("persistent", n);
        auto result_view = persistent_view.subview();
        contribute(result_view, persistent_view);
        Kokkos::fence();
      }
    }
  }
};

template <typename DeviceType, typename ScatterType, typename NumberType>
struct TestDuplicatedScatterView {
  TestDuplicatedScatterView(int n) {
    // ScatterSum test
    test_scatter_view_config<DeviceType, Kokkos::LayoutRight,
                             Kokkos::Experimental::ScatterDuplicated,
                             Kokkos::Experimental::ScatterNonAtomic,
                             ScatterType, NumberType>
        test_sv_right_config;
    test_sv_right_config.run_test(n);
    test_scatter_view_config<
        DeviceType, Kokkos::LayoutLeft, Kokkos::Experimental::ScatterDuplicated,
        Kokkos::Experimental::ScatterNonAtomic, ScatterType, NumberType>
        test_sv_left_config;
    test_sv_left_config.run_test(n);
  }
};

#ifdef KOKKOS_ENABLE_CUDA
// disable duplicated instantiation with CUDA until
// UniqueToken can support it
template <typename ScatterType, typename NumberType>
struct TestDuplicatedScatterView<Kokkos::Cuda, ScatterType, NumberType> {
  TestDuplicatedScatterView(int) {}
};
template <typename ScatterType, typename NumberType>
struct TestDuplicatedScatterView<
    Kokkos::Device<Kokkos::Cuda, Kokkos::CudaSpace>, ScatterType, NumberType> {
  TestDuplicatedScatterView(int) {}
};
template <typename ScatterType, typename NumberType>
struct TestDuplicatedScatterView<
    Kokkos::Device<Kokkos::Cuda, Kokkos::CudaUVMSpace>, ScatterType,
    NumberType> {
  TestDuplicatedScatterView(int) {}
};
#endif

template <typename DeviceType, typename ScatterType,
          typename NumberType = double>
void test_scatter_view(int64_t n) {
  using execution_space = typename DeviceType::execution_space;

  // no atomics or duplication is only sensible if the execution space
  // is running essentially in serial (doesn't have to be Serial though,
  // we also test OpenMP with one thread: LAMMPS cares about that)
  if (execution_space().concurrency() == 1) {
    test_scatter_view_config<DeviceType, Kokkos::LayoutRight,
                             Kokkos::Experimental::ScatterNonDuplicated,
                             Kokkos::Experimental::ScatterNonAtomic,
                             ScatterType, NumberType>
        test_sv_config;
    test_sv_config.run_test(n);
  }
#ifdef KOKKOS_ENABLE_SERIAL
  if (!std::is_same<DeviceType, Kokkos::Serial>::value) {
#endif
    test_scatter_view_config<DeviceType, Kokkos::LayoutRight,
                             Kokkos::Experimental::ScatterNonDuplicated,
                             Kokkos::Experimental::ScatterAtomic, ScatterType,
                             NumberType>
        test_sv_config;
    test_sv_config.run_test(n);
#ifdef KOKKOS_ENABLE_SERIAL
  }
#endif
  // with hundreds of threads we were running out of memory.
  // limit (n) so that duplication doesn't exceed 4GB
  constexpr std::size_t maximum_allowed_total_bytes =
      4ull * 1024ull * 1024ull * 1024ull;
  std::size_t const maximum_allowed_copy_bytes =
      maximum_allowed_total_bytes /
      std::size_t(execution_space().concurrency());
  constexpr std::size_t bytes_per_value = sizeof(NumberType) * 12;
  std::size_t const maximum_allowed_copy_values =
      maximum_allowed_copy_bytes / bytes_per_value;
  n = std::min(n, int64_t(maximum_allowed_copy_values));

  // if the default is duplicated, this needs to follow the limit
  {
    test_default_scatter_view<DeviceType, Kokkos::LayoutRight, ScatterType,
                              NumberType>
        test_default_sv;
    test_default_sv.run_test(n);
  }

  // run same test but on a subview (this covers support for padded
  // ScatterViews)
  {
    test_default_scatter_sub_view<DeviceType, Kokkos::LayoutRight, ScatterType,
                                  NumberType>
        test_default_scatter_view_subview;
    test_default_scatter_view_subview.run_test(n);
  }

  TestDuplicatedScatterView<DeviceType, ScatterType, NumberType> duptest(n);
}

TEST(TEST_CATEGORY, scatterview) {
  test_scatter_view<TEST_EXECSPACE, Kokkos::Experimental::ScatterSum, double>(
      10);

  test_scatter_view<TEST_EXECSPACE, Kokkos::Experimental::ScatterSum, int>(10);
  test_scatter_view<TEST_EXECSPACE, Kokkos::Experimental::ScatterProd>(10);
  test_scatter_view<TEST_EXECSPACE, Kokkos::Experimental::ScatterMin>(10);
  test_scatter_view<TEST_EXECSPACE, Kokkos::Experimental::ScatterMax>(10);
  // tests were timing out in DEBUG mode, reduce the amount of work
#ifdef KOKKOS_ENABLE_DEBUG
  int big_n = 100 * 1000;
#else

#ifdef KOKKOS_ENABLE_SERIAL
  bool is_serial = std::is_same<TEST_EXECSPACE, Kokkos::Serial>::value;
  int big_n      = is_serial ? 100 * 1000 : 10000 * 1000;
#else
  int big_n = 10000 * 1000;
#endif

#endif

  test_scatter_view<TEST_EXECSPACE, Kokkos::Experimental::ScatterSum, double>(
      big_n);
  test_scatter_view<TEST_EXECSPACE, Kokkos::Experimental::ScatterSum, int>(
      big_n);
  test_scatter_view<TEST_EXECSPACE, Kokkos::Experimental::ScatterProd>(big_n);
  test_scatter_view<TEST_EXECSPACE, Kokkos::Experimental::ScatterMin>(big_n);
  test_scatter_view<TEST_EXECSPACE, Kokkos::Experimental::ScatterMax>(big_n);
}

TEST(TEST_CATEGORY, scatterview_devicetype) {
  using device_type =
      Kokkos::Device<TEST_EXECSPACE, typename TEST_EXECSPACE::memory_space>;

  test_scatter_view<device_type, Kokkos::Experimental::ScatterSum, double>(10);
  test_scatter_view<device_type, Kokkos::Experimental::ScatterSum, int>(10);
  test_scatter_view<device_type, Kokkos::Experimental::ScatterProd>(10);
  test_scatter_view<device_type, Kokkos::Experimental::ScatterMin>(10);
  test_scatter_view<device_type, Kokkos::Experimental::ScatterMax>(10);

#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
#ifdef KOKKOS_ENABLE_CUDA
  using device_execution_space = Kokkos::Cuda;
  using device_memory_space    = Kokkos::CudaSpace;
  using host_accessible_space  = Kokkos::CudaUVMSpace;
#else
  using device_execution_space = Kokkos::Experimental::HIP;
  using device_memory_space    = Kokkos::Experimental::HIPSpace;
  using host_accessible_space  = Kokkos::Experimental::HIPHostPinnedSpace;
#endif
  if (std::is_same<TEST_EXECSPACE, device_execution_space>::value) {
    using device_device_type =
        Kokkos::Device<device_execution_space, device_memory_space>;
    test_scatter_view<device_device_type, Kokkos::Experimental::ScatterSum,
                      double>(10);
    test_scatter_view<device_device_type, Kokkos::Experimental::ScatterSum,
                      int>(10);
    test_scatter_view<device_device_type, Kokkos::Experimental::ScatterProd>(
        10);
    test_scatter_view<device_device_type, Kokkos::Experimental::ScatterMin>(10);
    test_scatter_view<device_device_type, Kokkos::Experimental::ScatterMax>(10);
    using host_device_type =
        Kokkos::Device<device_execution_space, host_accessible_space>;
    test_scatter_view<host_device_type, Kokkos::Experimental::ScatterSum,
                      double>(10);
    test_scatter_view<host_device_type, Kokkos::Experimental::ScatterSum, int>(
        10);
    test_scatter_view<host_device_type, Kokkos::Experimental::ScatterProd>(10);
    test_scatter_view<host_device_type, Kokkos::Experimental::ScatterMin>(10);
    test_scatter_view<host_device_type, Kokkos::Experimental::ScatterMax>(10);
  }
#endif
}

}  // namespace Test

#endif  // KOKKOS_TEST_SCATTER_VIEW_HPP
