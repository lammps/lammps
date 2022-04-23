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

#include <Kokkos_Core.hpp>
#include <TestHPX_Category.hpp>

#include <hpx/config.hpp>
#include <hpx/include/lcos.hpp>

#ifdef KOKKOS_ENABLE_HPX_ASYNC_DISPATCH
#ifndef HPX_COMPUTE_DEVICE_CODE

namespace Test {

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
}  // namespace

TEST(hpx, independent_instances) {
  Kokkos::InitArguments arguments{-1, -1, -1, false};
  Kokkos::initialize(arguments);

  const int n = 100;
  const int c = 1;
  const int d = 3;

  {
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

    Kokkos::Experimental::HPX hpx2(hpx1.impl_get_future());
    Kokkos::parallel_for(
        "Test::hpx::independent_instances::add",
        Kokkos::Experimental::require(
            Kokkos::RangePolicy<Kokkos::Experimental::HPX>(hpx2, 0, n),
            Kokkos::Experimental::WorkItemProperty::HintLightWeight),
        FunctorAdd(v1, v2, d));

    Kokkos::Experimental::HPX hpx3(hpx1.impl_get_future());
    Kokkos::parallel_for(
        "Test::hpx::independent_instances::add_index",
        Kokkos::Experimental::require(
            Kokkos::RangePolicy<Kokkos::Experimental::HPX>(hpx3, 0, n),
            Kokkos::Experimental::WorkItemProperty::HintLightWeight),
        FunctorAddIndex(v1, v3));

    // NOTE: This monstrosity is used to collapse a future<tuple<future<void>,
    // future<void>>> (return type of when_all) into a future<void> which is
    // ready whenever the un-collapsed future would've been ready. HPX does not
    // currently have the functionality to collapse this automatically.
    Kokkos::Experimental::HPX hpx4(hpx::util::get<0>(hpx::split_future(
        hpx::when_all(hpx2.impl_get_future(), hpx3.impl_get_future()))));
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

    ASSERT_EQ(true, hpx1.impl_get_future().is_ready());
    ASSERT_EQ(true, hpx2.impl_get_future().is_ready());
    ASSERT_EQ(true, hpx3.impl_get_future().is_ready());
    ASSERT_EQ(true, hpx4.impl_get_future().is_ready());

    const int expected_sum = n * (2 * c + d) + (n * (n - 1) / 2);
    ASSERT_EQ(expected_sum, sum_v());
  }

  Kokkos::finalize();
}
}  // namespace Test

#endif
#endif
