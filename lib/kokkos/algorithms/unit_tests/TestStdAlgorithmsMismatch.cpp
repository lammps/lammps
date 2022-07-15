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
#include <iterator>
#include <std_algorithms/Kokkos_BeginEnd.hpp>
#include <std_algorithms/Kokkos_NonModifyingSequenceOperations.hpp>
#include <algorithm>
#include <numeric>

namespace Test {
namespace stdalgos {
namespace Mismatch {

namespace KE = Kokkos::Experimental;

std::string value_type_to_string(int) { return "int"; }
std::string value_type_to_string(double) { return "double"; }

template <class Tag, class ValueType>
void print_scenario_details(std::size_t ext1, std::size_t ext2,
                            const std::string& flag) {
  std::cout << "mismatch: "
            << "ext1 = " << ext1 << ", "
            << "ext2 = " << ext2 << ", " << flag << ", "
            << view_tag_to_string(Tag{}) << ", "
            << value_type_to_string(ValueType()) << std::endl;
}

template <class Tag, class ViewType, class... Args>
void run_single_scenario(ViewType view1, ViewType view2,
                         const std::string& flag, Args... args) {
  using value_type = typename ViewType::value_type;
  using exe_space  = typename ViewType::execution_space;
  using aux_view_t = Kokkos::View<value_type*, exe_space>;

  const std::size_t ext1 = view1.extent(0);
  const std::size_t ext2 = view2.extent(0);
  // print_scenario_details<Tag, value_type>(ext1, ext2, flag);

  aux_view_t aux_view1("aux_view1", ext1);
  auto v1_h = create_mirror_view(Kokkos::HostSpace(), aux_view1);
  aux_view_t aux_view2("aux_view2", ext2);
  auto v2_h = create_mirror_view(Kokkos::HostSpace(), aux_view2);

  // note that the checks ext1>0 and ext2>0 are there
  // otherwise we get an error for CUDA NVCC DEBUG CI

  // view is is always filled with 8's
  if (ext1 > 0) {
    for (std::size_t i = 0; i < ext1; ++i) {
      v1_h(i) = static_cast<value_type>(8);
    }
  }

  if (flag == "fill-to-match") {
    if (ext2 > 0) {
      for (std::size_t i = 0; i < ext2; ++i) {
        v2_h(i) = static_cast<value_type>(8);
      }
    }
  }

  else if (flag == "fill-to-mismatch") {
    // need to make them mismatch, so we fill
    // with same value and only modifify the
    // second view arbitrarily at middle point

    if (ext2 > 0) {
      for (std::size_t i = 0; i < ext2; ++i) {
        v2_h(i) = static_cast<value_type>(8);
      }

      // make them mismatch at middle
      v2_h(ext2 / 2) = -5;
    }
  } else {
    throw std::runtime_error("Kokkos: stdalgo: test: mismatch: Invalid string");
  }

  Kokkos::deep_copy(aux_view1, v1_h);
  CopyFunctor<aux_view_t, ViewType> F1(aux_view1, view1);
  Kokkos::parallel_for("copy1", view1.extent(0), F1);

  Kokkos::deep_copy(aux_view2, v2_h);
  CopyFunctor<aux_view_t, ViewType> F2(aux_view2, view2);
  Kokkos::parallel_for("copy2", view2.extent(0), F2);

  // run the std::mismatch on a host copy of the data
  auto view1_h         = create_host_space_copy(view1);
  auto view2_h         = create_host_space_copy(view2);
  auto f1_h            = KE::cbegin(view1_h);
  auto l1_h            = KE::cend(view1_h);
  auto f2_h            = KE::cbegin(view2_h);
  auto l2_h            = KE::cend(view2_h);
  auto std_res         = std::mismatch(f1_h, l1_h, f2_h, l2_h, args...);
  const auto std_diff1 = std_res.first - f1_h;
  const auto std_diff2 = std_res.second - f2_h;

  {
    // check our overloads with iterators
    auto f1      = KE::cbegin(view1);
    auto l1      = KE::cend(view1);
    auto f2      = KE::cbegin(view2);
    auto l2      = KE::cend(view2);
    auto my_res1 = KE::mismatch(exespace(), f1, l1, f2, l2, args...);
    auto my_res2 = KE::mismatch("label", exespace(), f1, l1, f2, l2, args...);
    const auto my_diff11 = my_res1.first - f1;
    const auto my_diff12 = my_res1.second - f2;
    const auto my_diff21 = my_res2.first - f1;
    const auto my_diff22 = my_res2.second - f2;
    EXPECT_TRUE(my_diff11 == std_diff1);
    EXPECT_TRUE(my_diff12 == std_diff2);
    EXPECT_TRUE(my_diff21 == std_diff1);
    EXPECT_TRUE(my_diff22 == std_diff2);
  }

  {
    // check our overloads with views
    auto my_res1 = KE::mismatch(exespace(), view1, view2, args...);
    auto my_res2 = KE::mismatch("label", exespace(), view1, view2, args...);
    const auto my_diff11 = my_res1.first - KE::begin(view1);
    const auto my_diff12 = my_res1.second - KE::begin(view2);
    const auto my_diff21 = my_res2.first - KE::begin(view1);
    const auto my_diff22 = my_res2.second - KE::begin(view2);
    EXPECT_TRUE(my_diff11 == std_diff1);
    EXPECT_TRUE(my_diff12 == std_diff2);
    EXPECT_TRUE(my_diff21 == std_diff1);
    EXPECT_TRUE(my_diff22 == std_diff2);
  }
}

template <class Tag, class ValueType>
void run_all_scenarios() {
  using vecs_t = std::vector<std::string>;

  const std::map<std::string, std::size_t> scenarios = {
      {"empty", 0},  {"one-element", 1}, {"two-elements", 2},
      {"small", 11}, {"medium", 21103},  {"large", 101513}};

  for (const auto& scenario : scenarios) {
    {
      const std::size_t view1_ext = scenario.second;
      auto view1 = create_view<ValueType>(Tag{}, view1_ext, "mismatch_view_1");

      // for each view1 scenario, I want to test the case of a
      // second view that is smaller, equal size and greater than the view1
      const vecs_t view2cases = (scenario.first != "empty")
                                    ? vecs_t({"smaller", "equalsize", "larger"})
                                    : vecs_t({"equalsize", "larger"});

      for (auto it2 : view2cases) {
        std::size_t view2_ext = view1_ext;

        // modify extent of view2 based on what we want
        if (std::string(it2) == "smaller") {
          view2_ext -= 1;
        } else if (std::string(it2) == "larger") {
          view2_ext += 3;
        }

        auto view2 =
            create_view<ValueType>(Tag{}, view2_ext, "mismatch_view_2");

        // and now we want to test both the case view1 and view2 match,
        // as well as the case where they don't match
        for (const auto& it3 : {"fill-to-match", "fill-to-mismatch"}) {
          // run to use default predicate
          run_single_scenario<Tag>(view1, view2, it3);

          // run using an arbitrary predicate
          using predicate_type = IsEqualFunctor<ValueType>;
          run_single_scenario<Tag>(view1, view2, it3, predicate_type());
        }
      }
    }
  }
}

TEST(std_algorithms_mismatch_test, test) {
  run_all_scenarios<DynamicTag, double>();
  run_all_scenarios<StridedThreeTag, int>();
}

}  // namespace Mismatch
}  // namespace stdalgos
}  // namespace Test
