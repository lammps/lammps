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
#include <algorithm>

namespace Test {
namespace stdalgos {
namespace ForEach {

namespace KE = Kokkos::Experimental;

template <class ViewType>
void test_for_each(const ViewType view) {
  using value_t           = typename ViewType::value_type;
  using view_host_space_t = Kokkos::View<value_t*, Kokkos::HostSpace>;

  view_host_space_t expected("for_each_expected", view.extent(0));
  compare_views(expected, view);

  const auto mod_functor = IncrementElementWiseFunctor<value_t>();

  // pass view, functor takes non-const ref
  KE::for_each("label", exespace(), view, mod_functor);
  std::for_each(KE::begin(expected), KE::end(expected), mod_functor);
  compare_views(expected, view);

  // pass iterators, functor takes non-const ref
  KE::for_each(exespace(), KE::begin(view), KE::end(view), mod_functor);
  std::for_each(KE::begin(expected), KE::end(expected), mod_functor);
  compare_views(expected, view);

  const auto non_mod_functor = NoOpNonMutableFunctor<value_t>();

  // pass view, functor takes const ref
  KE::for_each(exespace(), view, non_mod_functor);
  std::for_each(KE::begin(expected), KE::end(expected), non_mod_functor);
  compare_views(expected, view);

  // pass const iterators, functor takes const ref
  KE::for_each(exespace(), KE::cbegin(view), KE::cend(view), non_mod_functor);
  std::for_each(KE::begin(expected), KE::end(expected), non_mod_functor);
  compare_views(expected, view);

#if defined(KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA)
  const auto mod_lambda = KOKKOS_LAMBDA(value_t & i) { ++i; };

  // pass view, lambda takes non-const ref
  KE::for_each(exespace(), view, mod_lambda);
  std::for_each(KE::begin(expected), KE::end(expected), mod_lambda);
  compare_views(expected, view);

  // pass iterators, lambda takes non-const ref
  KE::for_each(exespace(), KE::begin(view), KE::end(view), mod_lambda);
  std::for_each(KE::begin(expected), KE::end(expected), mod_lambda);
  compare_views(expected, view);

  const auto non_mod_lambda = KOKKOS_LAMBDA(const value_t& i) { (void)i; };

  // pass view, lambda takes const ref
  KE::for_each(exespace(), view, non_mod_lambda);
  std::for_each(KE::cbegin(expected), KE::cend(expected), non_mod_lambda);
  compare_views(expected, view);

  // pass const iterators, lambda takes const ref
  KE::for_each(exespace(), KE::cbegin(view), KE::cend(view), non_mod_lambda);
  std::for_each(KE::cbegin(expected), KE::cend(expected), non_mod_lambda);
  compare_views(expected, view);
#endif
}

// std::for_each_n is C++17, so we cannot compare results directly
template <class ViewType>
void test_for_each_n(const ViewType view) {
  using value_t       = typename ViewType::value_type;
  const std::size_t n = view.extent(0);

  const auto non_mod_functor = NoOpNonMutableFunctor<value_t>();

  // pass const iterators, functor takes const ref
  EXPECT_EQ(KE::cbegin(view) + n,
            KE::for_each_n(exespace(), KE::cbegin(view), n, non_mod_functor));
  verify_values(value_t{0}, view);

  // pass view, functor takes const ref
  EXPECT_EQ(KE::begin(view) + n,
            KE::for_each_n(exespace(), view, n, non_mod_functor));
  verify_values(value_t{0}, view);

  // pass iterators, functor takes non-const ref
  const auto mod_functor = IncrementElementWiseFunctor<value_t>();
  EXPECT_EQ(KE::begin(view) + n,
            KE::for_each_n(exespace(), KE::begin(view), n, mod_functor));
  verify_values(value_t{1}, view);

  // pass view, functor takes non-const ref
  EXPECT_EQ(KE::begin(view) + n,
            KE::for_each_n("label", exespace(), view, n, mod_functor));
  verify_values(value_t{2}, view);
}

template <class Tag, class ValueType>
void run_all_scenarios() {
  for (const auto& scenario : default_scenarios) {
    {
      auto view = create_view<ValueType>(Tag{}, scenario.second, "for_each");
      test_for_each(view);
    }
    {
      auto view = create_view<ValueType>(Tag{}, scenario.second, "for_each_n");
      test_for_each_n(view);
    }
  }
}

TEST(std_algorithms_for_each_test, test) {
  run_all_scenarios<DynamicTag, double>();
  run_all_scenarios<StridedTwoTag, int>();
  run_all_scenarios<StridedThreeTag, unsigned>();
}

}  // namespace ForEach
}  // namespace stdalgos
}  // namespace Test
