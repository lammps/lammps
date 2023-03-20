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

#include <TestStdAlgorithmsCommon.hpp>

namespace Test {
namespace stdalgos {
namespace ModOps {

namespace KE = Kokkos::Experimental;

// ------------
// move
// ------------
struct MyMovableType {
  int m_value = 11;

  MyMovableType() = default;
  MyMovableType(MyMovableType&& other) {
    if (this != &other) {
      m_value       = other.m_value;
      other.m_value = -2;
    }
  }

  MyMovableType& operator=(MyMovableType&& other) {
    if (this != &other) {
      m_value       = other.m_value;
      other.m_value = -4;
    }
    return *this;
  }
};

TEST(std_algorithms_mod_ops_test, move) {
  MyMovableType a;
  using move_t = decltype(std::move(a));
  static_assert(std::is_rvalue_reference<move_t>::value, "");

  // move constr
  MyMovableType b(std::move(a));
  EXPECT_EQ(b.m_value, 11);
  EXPECT_EQ(a.m_value, -2);

  // move assign
  MyMovableType c;
  c = std::move(b);
  EXPECT_EQ(c.m_value, 11);
  EXPECT_EQ(b.m_value, -4);
}

template <class ViewType>
struct StdAlgoModSeqOpsTestMove {
  ViewType m_view;

  KOKKOS_INLINE_FUNCTION
  void operator()(const int index) const {
    typename ViewType::value_type a{11};
    using move_t = decltype(std::move(a));
    static_assert(std::is_rvalue_reference<move_t>::value, "");
    m_view(index) = std::move(a);
  }

  StdAlgoModSeqOpsTestMove(ViewType view) : m_view(view) {}
};

TEST(std_algorithms_mod_ops_test, move_within_parfor) {
  using view_t = Kokkos::View<double*>;
  view_t a("a", 10);

  StdAlgoModSeqOpsTestMove<view_t> fnc(a);
  Kokkos::parallel_for(a.extent(0), fnc);
  auto a_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), a);
  for (std::size_t i = 0; i < a.extent(0); ++i) {
    EXPECT_DOUBLE_EQ(a_h(0), 11.);
  }
}

// ------------
// swap
// ------------
TEST(std_algorithms_mod_ops_test, swap) {
  {
    int a = 1;
    int b = 2;
    KE::swap(a, b);
    EXPECT_EQ(a, 2);
    EXPECT_EQ(b, 1);
  }

  {
    double a = 3.;
    double b = 1.;
    KE::swap(a, b);
    EXPECT_DOUBLE_EQ(a, 1.);
    EXPECT_DOUBLE_EQ(b, 3.);
  }
}

template <class ViewType>
struct StdAlgoModSeqOpsTestSwap {
  ViewType m_view;

  KOKKOS_INLINE_FUNCTION
  void operator()(const int index) const {
    typename ViewType::value_type newval{11};
    KE::swap(m_view(index), newval);
  }

  StdAlgoModSeqOpsTestSwap(ViewType aIn) : m_view(aIn) {}
};

TEST(std_algorithms_mod_ops_test, swap_within_parfor) {
  auto a = create_view<double>(stdalgos::DynamicTag{}, 10, "a");
  StdAlgoModSeqOpsTestSwap<decltype(a)> fnc(a);
  Kokkos::parallel_for(a.extent(0), fnc);
  auto a_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), a);
  for (std::size_t i = 0; i < a.extent(0); ++i) {
    EXPECT_DOUBLE_EQ(a_h(0), 11.);
  }
}

// ------------
// iter_swap
// ------------
template <class ViewType>
void test_iter_swap(ViewType view) {
  /* fill view */
  auto F = AssignIndexFunctor<ViewType>(view);
  Kokkos::parallel_for(view.extent(0), std::move(F));

  /* call iter_swap */
  auto it1 = KE::begin(view);
  KE::iter_swap(it1, it1 + 3);
  KE::iter_swap(it1 + 4, it1 + 6);

  /* check result */
  using value_type = typename ViewType::value_type;
  auto a_dc        = create_deep_copyable_compatible_clone(view);
  auto a_h         = create_mirror_view_and_copy(Kokkos::HostSpace(), a_dc);
  EXPECT_EQ(view.extent_int(0), 10);
  EXPECT_EQ(a_h(0), value_type(3));
  EXPECT_EQ(a_h(1), value_type(1));
  EXPECT_EQ(a_h(2), value_type(2));
  EXPECT_EQ(a_h(3), value_type(0));
  EXPECT_EQ(a_h(4), value_type(6));
  EXPECT_EQ(a_h(5), value_type(5));
  EXPECT_EQ(a_h(6), value_type(4));
  EXPECT_EQ(a_h(7), value_type(7));
  EXPECT_EQ(a_h(8), value_type(8));
  EXPECT_EQ(a_h(9), value_type(9));
}

TEST(std_algorithms_mod_ops_test, iter_swap_static_view) {
  auto a = create_view<double>(stdalgos::DynamicTag{}, 10, "a");
  test_iter_swap(a);

  auto a1 = create_view<double>(stdalgos::StridedTwoTag{}, 10, "a1");
  test_iter_swap(a1);

  auto a2 = create_view<double>(stdalgos::StridedThreeTag{}, 10, "a2");
  test_iter_swap(a2);
}

}  // namespace ModOps
}  // namespace stdalgos
}  // namespace Test
