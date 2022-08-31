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
#include <std_algorithms/Kokkos_BeginEnd.hpp>
#include <std_algorithms/Kokkos_ModifyingSequenceOperations.hpp>
#include <utility>

namespace Test {
namespace stdalgos {
namespace Search {

namespace KE = Kokkos::Experimental;

template <class ValueType>
struct UnifDist;

template <>
struct UnifDist<int> {
  using dist_type = std::uniform_int_distribution<int>;
  std::mt19937 m_gen;
  dist_type m_dist;

  UnifDist() : m_dist(0, 20) { m_gen.seed(1034343); }
  UnifDist(int a, int b) : m_dist(a, b) { m_gen.seed(234343); }

  int operator()() { return m_dist(m_gen); }
};

template <class ViewType>
void fill_view(ViewType dest_view, const std::string& name) {
  using value_type      = typename ViewType::value_type;
  using exe_space       = typename ViewType::execution_space;
  const std::size_t ext = dest_view.extent(0);
  using aux_view_t      = Kokkos::View<value_type*, exe_space>;
  aux_view_t aux_view("aux_view", ext);
  auto v_h = create_mirror_view(Kokkos::HostSpace(), aux_view);

  if (name == "empty") {
    // no op
  }

  else if (name == "one-element-a") {
    v_h(0) = static_cast<value_type>(1);
  }

  else if (name == "one-element-b") {
    v_h(0) = static_cast<value_type>(2);
  }

  else if (name == "two-elements-a") {
    v_h(0) = static_cast<value_type>(1);
    v_h(1) = static_cast<value_type>(2);
  }

  else if (name == "two-elements-b") {
    v_h(0) = static_cast<value_type>(2);
    v_h(1) = static_cast<value_type>(-1);
  }

  else if (name == "three-elements-a") {
    v_h(0) = static_cast<value_type>(-1);
    v_h(1) = static_cast<value_type>(2);
    v_h(2) = static_cast<value_type>(2);
  }

  else if (name == "three-elements-b") {
    v_h(0) = static_cast<value_type>(3);
    v_h(1) = static_cast<value_type>(1);
    v_h(2) = static_cast<value_type>(3);
  }

  else if (name == "four-elements-a") {
    v_h(0) = static_cast<value_type>(-1);
    v_h(1) = static_cast<value_type>(2);
    v_h(2) = static_cast<value_type>(2);
    v_h(3) = static_cast<value_type>(4);
  }

  else if (name == "four-elements-b") {
    v_h(0) = static_cast<value_type>(1);
    v_h(1) = static_cast<value_type>(1);
    v_h(2) = static_cast<value_type>(1);
    v_h(3) = static_cast<value_type>(1);
  }

  else if (name == "small-a") {
    v_h(0)  = static_cast<value_type>(0);
    v_h(1)  = static_cast<value_type>(4);
    v_h(2)  = static_cast<value_type>(1);
    v_h(3)  = static_cast<value_type>(2);
    v_h(4)  = static_cast<value_type>(-1);
    v_h(5)  = static_cast<value_type>(4);
    v_h(6)  = static_cast<value_type>(1);
    v_h(7)  = static_cast<value_type>(2);
    v_h(8)  = static_cast<value_type>(2);
    v_h(9)  = static_cast<value_type>(4);
    v_h(10) = static_cast<value_type>(1);
  }

  else if (name == "small-b") {
    v_h(0)  = static_cast<value_type>(1);
    v_h(1)  = static_cast<value_type>(2);
    v_h(2)  = static_cast<value_type>(3);
    v_h(3)  = static_cast<value_type>(1);
    v_h(4)  = static_cast<value_type>(-1);
    v_h(5)  = static_cast<value_type>(-2);
    v_h(6)  = static_cast<value_type>(0);
    v_h(7)  = static_cast<value_type>(1);
    v_h(8)  = static_cast<value_type>(2);
    v_h(9)  = static_cast<value_type>(2);
    v_h(10) = static_cast<value_type>(5);
    v_h(11) = static_cast<value_type>(9);
    v_h(12) = static_cast<value_type>(8);
  }

  else {
    UnifDist<value_type> randObj;
    for (std::size_t i = 0; i < ext; ++i) {
      v_h(i) = randObj();
    }
  }

  Kokkos::deep_copy(aux_view, v_h);
  CopyFunctor<aux_view_t, ViewType> F1(aux_view, dest_view);
  Kokkos::parallel_for("copy", dest_view.extent(0), F1);
}

template <class ViewType>
auto create_seq_to_search(ViewType data_view, std::size_t seq_extent) {
  // for the search, we need to specify a sequence that we search for
  // within the original view/range.
  // to do this, rather than doing something purely random,
  // we use the view with the data, and select a subsequence.

  auto data_view_h            = create_host_space_copy(data_view);
  const auto data_view_extent = data_view.extent(0);

  using value_type = typename ViewType::value_type;
  using exe_space  = typename ViewType::execution_space;
  using seq_view_t = Kokkos::View<value_type*, exe_space>;
  seq_view_t seq_view("seq_view", seq_extent);
  auto seq_view_h = create_mirror_view(Kokkos::HostSpace(), seq_view);

  // when the target sequence is of same size as view, just fill
  // sequeunce with all values of the view
  if (seq_extent == data_view_extent) {
    for (std::size_t i = 0; i < seq_extent; ++i) {
      seq_view_h(i) = data_view_h(i);
    }
  } else {
    // if target sequence to fill is smaller, then we need to pick
    // a starting point to copy data from to make the the sequence.
    // we pick randomly between 0 and data_view_extent - seq_extent.
    // and fill the sequeunce data with the values copied from data view.

    using dist_type = std::uniform_int_distribution<int>;
    std::random_device r;
    // from this:
    // https://stackoverflow.com/questions/34490599/c11-how-to-set-seed-using-random
    std::seed_seq seed{r(), r(), r(), r(), r(), r()};
    std::mt19937 gen(seed);
    dist_type dist(0, data_view_extent - seq_extent);
    const auto start = dist(gen);
    // std::cout << "start= " << start << "\n";
    for (std::size_t i = 0; i < seq_extent; ++i) {
      seq_view_h(i) = data_view_h(start + i);
      // std::cout << "i= " << i << " " << seq_view_h(i) << "\n";
    }
  }

  Kokkos::deep_copy(seq_view, seq_view_h);
  return seq_view;
}

// search is only avai from c++17, so I have to put it here
template <class ForwardIt1, class ForwardIt2, class BinaryPredicate>
ForwardIt1 my_std_search(ForwardIt1 first, ForwardIt1 last, ForwardIt2 s_first,
                         ForwardIt2 s_last, BinaryPredicate p) {
  for (;; ++first) {
    ForwardIt1 it = first;
    for (ForwardIt2 s_it = s_first;; ++it, ++s_it) {
      if (s_it == s_last) {
        return first;
      }
      if (it == last) {
        return last;
      }
      if (!p(*it, *s_it)) {
        break;
      }
    }
  }
}

// search is only avai from c++17, so I have to put it here
template <class ForwardIt1, class ForwardIt2>
ForwardIt1 my_std_search(ForwardIt1 first, ForwardIt1 last, ForwardIt2 s_first,
                         ForwardIt2 s_last) {
  using value_type1 = typename ForwardIt1::value_type;
  using value_type2 = typename ForwardIt2::value_type;

  using pred_t = IsEqualFunctor<value_type1, value_type2>;
  return my_std_search(first, last, s_first, s_last, pred_t());
}

std::string value_type_to_string(int) { return "int"; }
std::string value_type_to_string(double) { return "double"; }

template <class Tag, class ValueType>
void print_scenario_details(const std::string& name, std::size_t seq_ext) {
  std::cout << "search: default predicate: " << name << ", "
            << "search_seq_ext = " << seq_ext << ", "
            << view_tag_to_string(Tag{}) << " "
            << value_type_to_string(ValueType()) << std::endl;
}

template <class Tag, class ValueType, class Predicate>
void print_scenario_details(const std::string& name, std::size_t seq_ext,
                            Predicate pred) {
  (void)pred;
  std::cout << "search: custom  predicate: " << name << ", "
            << "search_seq_ext = " << seq_ext << ", "
            << view_tag_to_string(Tag{}) << " "
            << value_type_to_string(ValueType()) << std::endl;
}

template <class Tag, class ValueType, class InfoType, class... Args>
void run_single_scenario(const InfoType& scenario_info, std::size_t seq_ext,
                         Args... args) {
  const auto name            = std::get<0>(scenario_info);
  const std::size_t view_ext = std::get<1>(scenario_info);
  // print_scenario_details<Tag, ValueType>(name, seq_ext, args...);

  auto view = create_view<ValueType>(Tag{}, view_ext, "search_test_view");
  fill_view(view, name);
  auto s_view = create_seq_to_search(view, seq_ext);

  // run std
  auto view_h   = create_host_space_copy(view);
  auto s_view_h = create_host_space_copy(s_view);
  auto stdrit =
      my_std_search(KE::cbegin(view_h), KE::cend(view_h), KE::cbegin(s_view_h),
                    KE::cend(s_view_h), args...);

  {
    auto myrit        = KE::search(exespace(), KE::cbegin(view), KE::cend(view),
                            KE::cbegin(s_view), KE::cend(s_view), args...);
    const auto mydiff = myrit - KE::cbegin(view);
    const auto stddiff = stdrit - KE::cbegin(view_h);
    EXPECT_TRUE(mydiff == stddiff);
  }

  {
    auto myrit =
        KE::search("label", exespace(), KE::cbegin(view), KE::cend(view),
                   KE::cbegin(s_view), KE::cend(s_view), args...);
    const auto mydiff  = myrit - KE::cbegin(view);
    const auto stddiff = stdrit - KE::cbegin(view_h);
    EXPECT_TRUE(mydiff == stddiff);
  }

  {
    auto myrit         = KE::search(exespace(), view, s_view, args...);
    const auto mydiff  = myrit - KE::begin(view);
    const auto stddiff = stdrit - KE::cbegin(view_h);
    EXPECT_TRUE(mydiff == stddiff);
  }

  {
    auto myrit         = KE::search("label", exespace(), view, s_view, args...);
    const auto mydiff  = myrit - KE::begin(view);
    const auto stddiff = stdrit - KE::cbegin(view_h);
    EXPECT_TRUE(mydiff == stddiff);
  }

  Kokkos::fence();
}

template <class Tag, class ValueType>
void run_all_scenarios() {
  const std::map<std::string, std::size_t> scenarios = {{"empty", 0},
                                                        {"one-element-a", 1},
                                                        {"one-element-b", 1},
                                                        {"two-elements-a", 2},
                                                        {"two-elements-b", 2},
                                                        {"three-elements-a", 3},
                                                        {"three-elements-b", 3},
                                                        {"four-elements-a", 4},
                                                        {"four-elements-b", 4},
                                                        {"small-a", 11},
                                                        {"small-b", 13},
                                                        {"medium-a", 11103},
                                                        {"medium-b", 21103},
                                                        {"large-a", 101513},
                                                        {"large-b", 100111}};

  const std::vector<std::size_t> seq_extents = {
      0, 1, 2, 3, 4, 5, 8, 11, 15, 31, 113, 523, 1035, 11103};

  // for each scenario we want to run "search"
  // for a set of sequences of various extents
  for (const auto& it : scenarios) {
    for (const auto& it2 : seq_extents) {
      // only run if view is larger or equal than sequence to search for
      if (it.second >= it2) {
        run_single_scenario<Tag, ValueType>(it, it2);

        using func_t = IsEqualFunctor<ValueType>;
        run_single_scenario<Tag, ValueType>(it, it2, func_t());
      }
    }
  }
}

TEST(std_algorithms_non_mod_seq_ops, search) {
  run_all_scenarios<DynamicTag, int>();
  run_all_scenarios<StridedThreeTag, int>();
}

}  // namespace Search
}  // namespace stdalgos
}  // namespace Test
