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

template <class Space>
struct NestedView {
  Kokkos::View<int *, Space> member;

 public:
  KOKKOS_INLINE_FUNCTION
  NestedView() : member() {}

  KOKKOS_INLINE_FUNCTION
  NestedView &operator=(const Kokkos::View<int *, Space> &lhs) {
    member = lhs;
    if (member.extent(0)) Kokkos::atomic_add(&member(0), 1);
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  ~NestedView() {
    if (member.extent(0)) {
      Kokkos::atomic_add(&member(0), -1);
    }
  }
};

template <class Space>
struct NestedViewFunctor {
  Kokkos::View<NestedView<Space> *, Space> nested;
  Kokkos::View<int *, Space> array;

  NestedViewFunctor(const Kokkos::View<NestedView<Space> *, Space> &arg_nested,
                    const Kokkos::View<int *, Space> &arg_array)
      : nested(arg_nested), array(arg_array) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(int i) const { nested[i] = array; }
};

template <class Space>
void view_nested_view() {
  Kokkos::View<int *, Space> tracking("tracking", 1);

  typename Kokkos::View<int *, Space>::HostMirror host_tracking =
      Kokkos::create_mirror(tracking);

  {
    Kokkos::View<NestedView<Space> *, Space> a("a_nested_view", 2);

    Kokkos::parallel_for(Kokkos::RangePolicy<Space>(0, 2),
                         NestedViewFunctor<Space>(a, tracking));
    Kokkos::deep_copy(host_tracking, tracking);
    ASSERT_EQ(2, host_tracking(0));

    Kokkos::View<NestedView<Space> *, Space> b("b_nested_view", 2);
    Kokkos::parallel_for(Kokkos::RangePolicy<Space>(0, 2),
                         NestedViewFunctor<Space>(b, tracking));
    Kokkos::deep_copy(host_tracking, tracking);
    ASSERT_EQ(4, host_tracking(0));
  }

  Kokkos::deep_copy(host_tracking, tracking);

  ASSERT_EQ(0, host_tracking(0));
}

TEST(TEST_CATEGORY, view_nested_view) { view_nested_view<TEST_EXECSPACE>(); }

}  // namespace Test
