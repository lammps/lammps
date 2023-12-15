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

#include <Kokkos_Core.hpp>
#include <TestHPX_Category.hpp>

#include <hpx/config.hpp>
#include <hpx/future.hpp>

#ifndef HPX_COMPUTE_DEVICE_CODE

namespace {

struct FunctorInitConstant {
  Kokkos::View<int *, Kokkos::Experimental::HPX> a;
  int c;
  FunctorInitConstant(Kokkos::View<int *, Kokkos::Experimental::HPX> a_, int c_)
      : a(a_), c(c_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const { a(i) = c; }
};

struct FunctorAdd {
  Kokkos::View<int *, Kokkos::Experimental::HPX> a;
  Kokkos::View<int *, Kokkos::Experimental::HPX> b;
  int c;
  FunctorAdd(Kokkos::View<int *, Kokkos::Experimental::HPX> a_,
             Kokkos::View<int *, Kokkos::Experimental::HPX> b_, int c_)
      : a(a_), b(b_), c(c_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const { b(i) += a(i) + c; }
};

struct FunctorAddIndex {
  Kokkos::View<int *, Kokkos::Experimental::HPX> a;
  Kokkos::View<int *, Kokkos::Experimental::HPX> b;
  FunctorAddIndex(Kokkos::View<int *, Kokkos::Experimental::HPX> a_,
                  Kokkos::View<int *, Kokkos::Experimental::HPX> b_)
      : a(a_), b(b_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const { b(i) += a(i) + i; }
};

struct FunctorPointwiseSum {
  Kokkos::View<int *, Kokkos::Experimental::HPX> a;
  Kokkos::View<int *, Kokkos::Experimental::HPX> b;
  Kokkos::View<int *, Kokkos::Experimental::HPX> c;
  FunctorPointwiseSum(Kokkos::View<int *, Kokkos::Experimental::HPX> a_,
                      Kokkos::View<int *, Kokkos::Experimental::HPX> b_,
                      Kokkos::View<int *, Kokkos::Experimental::HPX> c_)
      : a(a_), b(b_), c(c_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i) const { c(i) = a(i) + b(i); }
};

struct FunctorReduce {
  Kokkos::View<int *, Kokkos::Experimental::HPX> a;
  FunctorReduce(Kokkos::View<int *, Kokkos::Experimental::HPX> a_) : a(a_) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, int &lsum) const { lsum += a(i); }
};

TEST(hpx, independent_instances) {
  const int n = 100;
  const int c = 1;
  const int d = 3;

  Kokkos::View<int *, Kokkos::Experimental::HPX> v1("v1", n);
  Kokkos::View<int *, Kokkos::Experimental::HPX> v2("v2", n);
  Kokkos::View<int *, Kokkos::Experimental::HPX> v3("v3", n);
  Kokkos::View<int *, Kokkos::Experimental::HPX> v4("v4", n);
  Kokkos::View<int, Kokkos::Experimental::HPX> sum_v("sum_v");

  Kokkos::Experimental::HPX hpx1(
      Kokkos::Experimental::HPX::instance_mode::independent);
  Kokkos::parallel_for(
      "Test::hpx::independent_instances::init",
      Kokkos::Experimental::require(
          Kokkos::RangePolicy<Kokkos::Experimental::HPX>(hpx1, 0, n),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      FunctorInitConstant(v1, c));

  Kokkos::Experimental::HPX hpx2(hpx1.get_sender());
  Kokkos::parallel_for(
      "Test::hpx::independent_instances::add",
      Kokkos::Experimental::require(
          Kokkos::RangePolicy<Kokkos::Experimental::HPX>(hpx2, 0, n),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      FunctorAdd(v1, v2, d));

  Kokkos::Experimental::HPX hpx3(hpx1.get_sender());
  Kokkos::parallel_for(
      "Test::hpx::independent_instances::add_index",
      Kokkos::Experimental::require(
          Kokkos::RangePolicy<Kokkos::Experimental::HPX>(hpx3, 0, n),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      FunctorAddIndex(v1, v3));

  Kokkos::Experimental::HPX hpx4(hpx::execution::experimental::when_all(
      hpx2.get_sender(), hpx3.get_sender()));
  Kokkos::parallel_for(
      "Test::hpx::independent_instances::pointwise_sum",
      Kokkos::Experimental::require(
          Kokkos::RangePolicy<Kokkos::Experimental::HPX>(hpx4, 0, n),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      FunctorPointwiseSum(v2, v3, v4));

  Kokkos::parallel_reduce(
      "Test::hpx::independent_instances::reduce",
      Kokkos::Experimental::require(
          Kokkos::RangePolicy<Kokkos::Experimental::HPX>(hpx4, 0, n),
          Kokkos::Experimental::WorkItemProperty::HintLightWeight),
      FunctorReduce(v4), Kokkos::Sum<int>(sum_v));

  hpx4.fence();

  const int expected_sum = n * (2 * c + d) + (n * (n - 1) / 2);
  ASSERT_EQ(expected_sum, sum_v());
}

}  // namespace

#endif
