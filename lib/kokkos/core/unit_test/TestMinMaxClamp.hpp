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

// FIXME C++17
#define STATIC_ASSERT(cond) static_assert(cond, "");

namespace Test {
template <class T>
struct Greater {
  KOKKOS_FUNCTION constexpr bool operator()(T const& lhs, T const& rhs) {
    return lhs > rhs;
  }
};

struct PairIntCompareFirst {
  int first;
  int second;

 private:
  friend KOKKOS_FUNCTION constexpr bool operator<(
      PairIntCompareFirst const& lhs, PairIntCompareFirst const& rhs) {
    return lhs.first < rhs.first;
  }
};
}  // namespace Test

// ----------------------------------------------------------
// test max()
// ----------------------------------------------------------
TEST(TEST_CATEGORY, max) {
  int a = 1;
  int b = 2;
  EXPECT_EQ(Kokkos::max(a, b), 2);

  a = 3;
  b = 1;
  EXPECT_EQ(Kokkos::max(a, b), 3);

  STATIC_ASSERT(Kokkos::max(1, 2) == 2);
  STATIC_ASSERT(Kokkos::max(1, 2, ::Test::Greater<int>{}) == 1);

  EXPECT_EQ(Kokkos::max({3.f, -1.f, 0.f}), 3.f);

  STATIC_ASSERT(Kokkos::max({3, -1, 0}) == 3);
  STATIC_ASSERT(Kokkos::max({3, -1, 0}, ::Test::Greater<int>{}) == -1);

  STATIC_ASSERT(Kokkos::max({
                                ::Test::PairIntCompareFirst{255, 0},
                                ::Test::PairIntCompareFirst{255, 1},
                                ::Test::PairIntCompareFirst{0, 2},
                                ::Test::PairIntCompareFirst{0, 3},
                                ::Test::PairIntCompareFirst{255, 4},
                                ::Test::PairIntCompareFirst{0, 5},
                            })
                    .second == 0);  // leftmost element
}

template <class ViewType>
struct StdAlgoMinMaxOpsTestMax {
  ViewType m_view;

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& ind) const {
    auto v1 = 10.;
    if (Kokkos::max(v1, m_view(ind)) == 10.) {
      m_view(ind) = 6.;
    }
  }

  KOKKOS_INLINE_FUNCTION
  StdAlgoMinMaxOpsTestMax(ViewType aIn) : m_view(aIn) {}
};

TEST(TEST_CATEGORY, max_within_parfor) {
  namespace KE = Kokkos::Experimental;

  using view_t = Kokkos::View<double*>;
  view_t a("a", 10);

  StdAlgoMinMaxOpsTestMax<view_t> fnc(a);
  Kokkos::parallel_for(a.extent(0), fnc);
  auto a_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), a);
  for (int i = 0; i < 10; ++i) {
    EXPECT_DOUBLE_EQ(a_h(0), 6.);
  }
}

// ----------------------------------------------------------
// test min()
// ----------------------------------------------------------
TEST(TEST_CATEGORY, min) {
  int a = 1;
  int b = 2;
  EXPECT_EQ(Kokkos::min(a, b), 1);

  a = 3;
  b = 2;
  EXPECT_EQ(Kokkos::min(a, b), 2);

  STATIC_ASSERT(Kokkos::min(3.f, 2.f) == 2.f);
  STATIC_ASSERT(Kokkos::min(3.f, 2.f, ::Test::Greater<int>{}) == 3.f);

  EXPECT_EQ(Kokkos::min({3.f, -1.f, 0.f}), -1.f);

  STATIC_ASSERT(Kokkos::min({3, -1, 0}) == -1);
  STATIC_ASSERT(Kokkos::min({3, -1, 0}, ::Test::Greater<int>{}) == 3);

  STATIC_ASSERT(Kokkos::min({
                                ::Test::PairIntCompareFirst{255, 0},
                                ::Test::PairIntCompareFirst{255, 1},
                                ::Test::PairIntCompareFirst{0, 2},
                                ::Test::PairIntCompareFirst{0, 3},
                                ::Test::PairIntCompareFirst{255, 4},
                                ::Test::PairIntCompareFirst{0, 5},
                            })
                    .second == 2);  // leftmost element
}

template <class ViewType>
struct StdAlgoMinMaxOpsTestMin {
  ViewType m_view;

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& ind) const {
    auto v1 = 10.;
    if (Kokkos::min(v1, m_view(ind)) == 0.) {
      m_view(ind) = 8.;
    }
  }

  KOKKOS_INLINE_FUNCTION
  StdAlgoMinMaxOpsTestMin(ViewType aIn) : m_view(aIn) {}
};

TEST(TEST_CATEGORY, min_within_parfor) {
  namespace KE = Kokkos::Experimental;
  using view_t = Kokkos::View<double*>;
  view_t a("a", 10);

  StdAlgoMinMaxOpsTestMin<view_t> fnc(a);
  Kokkos::parallel_for(a.extent(0), fnc);
  auto a_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), a);
  for (int i = 0; i < 10; ++i) {
    EXPECT_DOUBLE_EQ(a_h(0), 8.);
  }
}

// ----------------------------------------------------------
// test minmax()
// ----------------------------------------------------------
TEST(TEST_CATEGORY, minmax) {
  int a         = 1;
  int b         = 2;
  const auto& r = Kokkos::minmax(a, b);
  EXPECT_EQ(r.first, 1);
  EXPECT_EQ(r.second, 2);

  a              = 3;
  b              = 2;
  const auto& r2 = Kokkos::minmax(a, b);
  EXPECT_EQ(r2.first, 2);
  EXPECT_EQ(r2.second, 3);

#ifndef KOKKOS_COMPILER_NVHPC  // FIXME_NVHPC nvhpc can't deal with device side
                               // constexpr constructors so I removed the
                               // constexpr in pair, which makes STATIC_ASSERT
                               // here fail
  STATIC_ASSERT((Kokkos::pair<float, float>(Kokkos::minmax(3.f, 2.f)) ==
                 Kokkos::make_pair(2.f, 3.f)));
  STATIC_ASSERT(
      (Kokkos::pair<float, float>(Kokkos::minmax(
           3.f, 2.f, ::Test::Greater<int>{})) == Kokkos::make_pair(3.f, 2.f)));

  EXPECT_EQ(Kokkos::minmax({3.f, -1.f, 0.f}), Kokkos::make_pair(-1.f, 3.f));

  STATIC_ASSERT(Kokkos::minmax({3, -1, 0}) == Kokkos::make_pair(-1, 3));
  STATIC_ASSERT(Kokkos::minmax({3, -1, 0}, ::Test::Greater<int>{}) ==
                Kokkos::make_pair(3, -1));

  STATIC_ASSERT(Kokkos::minmax({
                                   ::Test::PairIntCompareFirst{255, 0},
                                   ::Test::PairIntCompareFirst{255, 1},
                                   ::Test::PairIntCompareFirst{0, 2},
                                   ::Test::PairIntCompareFirst{0, 3},
                                   ::Test::PairIntCompareFirst{255, 4},
                                   ::Test::PairIntCompareFirst{0, 5},
                               })
                    .first.second == 2);  // leftmost
  STATIC_ASSERT(Kokkos::minmax({
                                   ::Test::PairIntCompareFirst{255, 0},
                                   ::Test::PairIntCompareFirst{255, 1},
                                   ::Test::PairIntCompareFirst{0, 2},
                                   ::Test::PairIntCompareFirst{0, 3},
                                   ::Test::PairIntCompareFirst{255, 4},
                                   ::Test::PairIntCompareFirst{0, 5},
                               })
                    .second.second == 4);  // rightmost
#endif
}

template <class ViewType>
struct StdAlgoMinMaxOpsTestMinMax {
  ViewType m_view;

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& ind) const {
    auto v1       = 7.;
    const auto& r = Kokkos::minmax(v1, m_view(ind));
    m_view(ind)   = (double)(r.first - r.second);
  }

  KOKKOS_INLINE_FUNCTION
  StdAlgoMinMaxOpsTestMinMax(ViewType aIn) : m_view(aIn) {}
};

TEST(TEST_CATEGORY, minmax_within_parfor) {
  using view_t = Kokkos::View<double*>;
  view_t a("a", 10);

  StdAlgoMinMaxOpsTestMinMax<view_t> fnc(a);
  Kokkos::parallel_for(a.extent(0), fnc);
  auto a_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), a);
  for (int i = 0; i < 10; ++i) {
    EXPECT_DOUBLE_EQ(a_h(0), -7.);
  }
}

// ----------------------------------------------------------
// test clamp()
// ----------------------------------------------------------
TEST(TEST_CATEGORY, clamp) {
  int a         = 1;
  int b         = 2;
  int c         = 19;
  const auto& r = Kokkos::clamp(a, b, c);
  EXPECT_EQ(&r, &b);
  EXPECT_EQ(r, b);

  a              = 5;
  b              = -2;
  c              = 3;
  const auto& r2 = Kokkos::clamp(a, b, c);
  EXPECT_EQ(&r2, &c);
  EXPECT_EQ(r2, c);

  a              = 5;
  b              = -2;
  c              = 7;
  const auto& r3 = Kokkos::clamp(a, b, c);
  EXPECT_EQ(&r3, &a);
  EXPECT_EQ(r3, a);
}

template <class ViewType>
struct StdAlgoMinMaxOpsTestClamp {
  ViewType m_view;

  KOKKOS_INLINE_FUNCTION
  void operator()(const int& ind) const {
    m_view(ind)   = 10.;
    const auto b  = -2.;
    const auto c  = 3.;
    const auto& r = Kokkos::clamp(m_view(ind), b, c);
    m_view(ind)   = (double)(r);
  }

  KOKKOS_INLINE_FUNCTION
  StdAlgoMinMaxOpsTestClamp(ViewType aIn) : m_view(aIn) {}
};

TEST(TEST_CATEGORY, clamp_within_parfor) {
  using view_t = Kokkos::View<double*>;
  view_t a("a", 10);

  StdAlgoMinMaxOpsTestClamp<view_t> fnc(a);
  Kokkos::parallel_for(a.extent(0), fnc);
  auto a_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), a);
  for (std::size_t i = 0; i < a.extent(0); ++i) {
    EXPECT_DOUBLE_EQ(a_h(0), 3.);
  }
}
