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
#include <utility>

namespace Test {
namespace stdalgos {
namespace UniqueCopy {

namespace KE = Kokkos::Experimental;

// impl is here for std because it is only avail from c++>=17
template <class InputIterator, class OutputIterator, class BinaryPredicate>
auto my_unique_copy(InputIterator first, InputIterator last,
                    OutputIterator result, BinaryPredicate pred) {
  if (first != last) {
    typename OutputIterator::value_type t(*first);
    *result = t;
    ++result;
    while (++first != last) {
      if (!pred(t, *first)) {
        t       = *first;
        *result = t;
        ++result;
      }
    }
  }
  return result;
}

template <class InputIterator, class OutputIterator>
auto my_unique_copy(InputIterator first, InputIterator last,
                    OutputIterator result) {
  using value_type = typename OutputIterator::value_type;
  using func_t     = IsEqualFunctor<value_type>;
  return my_unique_copy(first, last, result, func_t());
}

template <class ValueType>
struct UnifDist;

template <>
struct UnifDist<int> {
  using dist_type = std::uniform_int_distribution<int>;
  std::mt19937 m_gen;
  dist_type m_dist;

  UnifDist() : m_dist(2, 7) { m_gen.seed(1034343); }

  int operator()() { return m_dist(m_gen); }
};

template <class ViewType>
std::size_t fill_view(ViewType dest_view, const std::string& name) {
  using value_type      = typename ViewType::value_type;
  using exe_space       = typename ViewType::execution_space;
  const std::size_t ext = dest_view.extent(0);
  using aux_view_t      = Kokkos::View<value_type*, exe_space>;
  aux_view_t aux_view("aux_view", ext);
  auto v_h = create_mirror_view(Kokkos::HostSpace(), aux_view);

  std::size_t count = 0;
  if (name == "empty") {
    // no op
  }

  else if (name == "one-element-a") {
    v_h(0) = static_cast<value_type>(1);
    count  = 1;
  }

  else if (name == "one-element-b") {
    v_h(0) = static_cast<value_type>(2);
    count  = 1;
  }

  else if (name == "two-elements-a") {
    v_h(0) = static_cast<value_type>(1);
    v_h(1) = static_cast<value_type>(2);
    count  = 2;
  }

  else if (name == "two-elements-b") {
    v_h(0) = static_cast<value_type>(2);
    v_h(1) = static_cast<value_type>(-1);
    count  = 2;
  }

  else if (name == "small-a") {
    v_h(0)  = static_cast<value_type>(0);
    v_h(1)  = static_cast<value_type>(1);
    v_h(2)  = static_cast<value_type>(1);
    v_h(3)  = static_cast<value_type>(2);
    v_h(4)  = static_cast<value_type>(3);
    v_h(5)  = static_cast<value_type>(4);
    v_h(6)  = static_cast<value_type>(4);
    v_h(7)  = static_cast<value_type>(4);
    v_h(8)  = static_cast<value_type>(5);
    v_h(9)  = static_cast<value_type>(6);
    v_h(10) = static_cast<value_type>(6);
    count   = 7;
  }

  else if (name == "small-b") {
    v_h(0)  = static_cast<value_type>(1);
    v_h(1)  = static_cast<value_type>(1);
    v_h(2)  = static_cast<value_type>(1);
    v_h(3)  = static_cast<value_type>(2);
    v_h(4)  = static_cast<value_type>(3);
    v_h(5)  = static_cast<value_type>(4);
    v_h(6)  = static_cast<value_type>(4);
    v_h(7)  = static_cast<value_type>(4);
    v_h(8)  = static_cast<value_type>(5);
    v_h(9)  = static_cast<value_type>(6);
    v_h(10) = static_cast<value_type>(8);
    v_h(11) = static_cast<value_type>(9);
    v_h(12) = static_cast<value_type>(8);
    count   = 9;
  }

  else if (name == "medium" || name == "large") {
    UnifDist<value_type> randObj;
    for (std::size_t i = 0; i < ext; ++i) {
      v_h(i) = randObj();
    }
    std::vector<value_type> tmp(v_h.extent(0));
    std::fill(tmp.begin(), tmp.end(), static_cast<value_type>(0));
    using func_t = IsEqualFunctor<value_type>;
    auto std_r =
        my_unique_copy(KE::cbegin(v_h), KE::cend(v_h), tmp.begin(), func_t());
    count = (std::size_t)(std_r - tmp.begin());
  }

  else {
    throw std::runtime_error("invalid choice");
  }

  Kokkos::deep_copy(aux_view, v_h);
  CopyFunctor<aux_view_t, ViewType> F1(aux_view, dest_view);
  Kokkos::parallel_for("copy", dest_view.extent(0), F1);
  return count;
}

template <class ViewTypeFrom, class ViewTypeTest, class... Args>
void verify_data(const std::string& name, ViewTypeFrom view_from,
                 ViewTypeTest view_test, Args... args) {
  using value_type = typename ViewTypeTest::value_type;

  //! always careful because views might not be deep copyable
  auto view_test_dc = create_deep_copyable_compatible_clone(view_test);
  auto view_test_h =
      create_mirror_view_and_copy(Kokkos::HostSpace(), view_test_dc);

  auto view_from_dc = create_deep_copyable_compatible_clone(view_from);
  auto view_from_h =
      create_mirror_view_and_copy(Kokkos::HostSpace(), view_from_dc);

  if (name == "empty") {
    // no op
  }

  else if (name == "one-element-a") {
    EXPECT_EQ(view_test_h(0), static_cast<value_type>(1));
  }

  else if (name == "one-element-b") {
    EXPECT_EQ(view_test_h(0), static_cast<value_type>(2));
  }

  else if (name == "two-elements-a") {
    EXPECT_EQ(view_test_h(0), static_cast<value_type>(1));
    EXPECT_EQ(view_test_h(1), static_cast<value_type>(2));
  }

  else if (name == "two-elements-b") {
    EXPECT_EQ(view_test_h(0), static_cast<value_type>(2));
    EXPECT_EQ(view_test_h(1), static_cast<value_type>(-1));
  }

  else if (name == "small-a") {
    EXPECT_EQ(view_test_h(0), static_cast<value_type>(0));
    EXPECT_EQ(view_test_h(1), static_cast<value_type>(1));
    EXPECT_EQ(view_test_h(2), static_cast<value_type>(2));
    EXPECT_EQ(view_test_h(3), static_cast<value_type>(3));
    EXPECT_EQ(view_test_h(4), static_cast<value_type>(4));
    EXPECT_EQ(view_test_h(5), static_cast<value_type>(5));
    EXPECT_EQ(view_test_h(6), static_cast<value_type>(6));
    EXPECT_EQ(view_test_h(7), static_cast<value_type>(0));
    EXPECT_EQ(view_test_h(8), static_cast<value_type>(0));
    EXPECT_EQ(view_test_h(9), static_cast<value_type>(0));
    EXPECT_EQ(view_test_h(10), static_cast<value_type>(0));
  }

  else if (name == "small-b") {
    EXPECT_EQ(view_test_h(0), static_cast<value_type>(1));
    EXPECT_EQ(view_test_h(1), static_cast<value_type>(2));
    EXPECT_EQ(view_test_h(2), static_cast<value_type>(3));
    EXPECT_EQ(view_test_h(3), static_cast<value_type>(4));
    EXPECT_EQ(view_test_h(4), static_cast<value_type>(5));
    EXPECT_EQ(view_test_h(5), static_cast<value_type>(6));
    EXPECT_EQ(view_test_h(6), static_cast<value_type>(8));
    EXPECT_EQ(view_test_h(7), static_cast<value_type>(9));
    EXPECT_EQ(view_test_h(8), static_cast<value_type>(8));
    EXPECT_EQ(view_test_h(9), static_cast<value_type>(0));
    EXPECT_EQ(view_test_h(10), static_cast<value_type>(0));
    EXPECT_EQ(view_test_h(11), static_cast<value_type>(0));
    EXPECT_EQ(view_test_h(12), static_cast<value_type>(0));
  }

  else if (name == "medium" || name == "large") {
    std::vector<value_type> tmp(view_test_h.extent(0));
    std::fill(tmp.begin(), tmp.end(), static_cast<value_type>(0));

    auto std_r = my_unique_copy(KE::cbegin(view_from_h), KE::cend(view_from_h),
                                tmp.begin(), args...);
    (void)std_r;

    for (std::size_t i = 0; i < view_from_h.extent(0); ++i) {
      EXPECT_EQ(view_test_h(i), tmp[i]);
    }
  }

  else {
    throw std::runtime_error("invalid choice");
  }
}

std::string value_type_to_string(int) { return "int"; }
std::string value_type_to_string(double) { return "double"; }

template <class Tag, class ValueType>
void print_scenario_details(const std::string& name) {
  std::cout << "unique_copy: default predicate: " << name << ", "
            << view_tag_to_string(Tag{}) << " "
            << value_type_to_string(ValueType()) << '\n';
}

template <class Tag, class ValueType, class Predicate>
void print_scenario_details(const std::string& name, Predicate pred) {
  (void)pred;
  std::cout << "unique_copy: custom  predicate: " << name << ", "
            << view_tag_to_string(Tag{}) << " "
            << value_type_to_string(ValueType()) << '\n';
}

template <class Tag, class ValueType, class InfoType, class... Args>
void run_single_scenario(const InfoType& scenario_info, Args... args) {
  const auto name            = std::get<0>(scenario_info);
  const std::size_t view_ext = std::get<1>(scenario_info);
  // print_scenario_details<Tag, ValueType>(name, args...);

  auto view_from = create_view<ValueType>(Tag{}, view_ext, "unique_copy_from");
  auto n         = fill_view(view_from, name);

  {
    auto view_dest =
        create_view<ValueType>(Tag{}, view_ext, "unique_copy_dest");
    auto rit =
        KE::unique_copy(exespace(), KE::cbegin(view_from), KE::cend(view_from),
                        KE::begin(view_dest), args...);
    verify_data(name, view_from, view_dest, args...);
    EXPECT_EQ(rit, (KE::begin(view_dest) + n));
  }

  {
    auto view_dest =
        create_view<ValueType>(Tag{}, view_ext, "unique_copy_dest");
    auto rit =
        KE::unique_copy("label", exespace(), KE::cbegin(view_from),
                        KE::cend(view_from), KE::begin(view_dest), args...);
    verify_data(name, view_from, view_dest, args...);
    EXPECT_EQ(rit, (KE::begin(view_dest) + n));
  }

  {
    auto view_dest =
        create_view<ValueType>(Tag{}, view_ext, "unique_copy_dest");
    auto rit = KE::unique_copy(exespace(), view_from, view_dest, args...);
    verify_data(name, view_from, view_dest, args...);
    EXPECT_EQ(rit, (KE::begin(view_dest) + n));
  }

  {
    auto view_dest =
        create_view<ValueType>(Tag{}, view_ext, "unique_copy_dest");
    auto rit =
        KE::unique_copy("label", exespace(), view_from, view_dest, args...);
    verify_data(name, view_from, view_dest, args...);
    EXPECT_EQ(rit, (KE::begin(view_dest) + n));
  }

  Kokkos::fence();
}

template <class Tag, class ValueType>
void run_all_scenarios() {
  const std::map<std::string, std::size_t> scenarios = {
      {"empty", 0},          {"one-element-a", 1},  {"one-element-b", 1},
      {"two-elements-a", 2}, {"two-elements-b", 2}, {"small-a", 11},
      {"small-b", 13},       {"medium", 21103},     {"large", 101513}};

  for (const auto& it : scenarios) {
    run_single_scenario<Tag, ValueType>(it);

    using func_t = IsEqualFunctor<ValueType>;
    run_single_scenario<Tag, ValueType>(it, func_t());
  }
}

TEST(std_algorithms_mod_seq_ops, unique_copy) {
  run_all_scenarios<DynamicTag, int>();
  run_all_scenarios<StridedThreeTag, int>();
}

}  // namespace UniqueCopy
}  // namespace stdalgos
}  // namespace Test
