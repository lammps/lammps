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
#include <Kokkos_Random.hpp>

namespace Test {
namespace stdalgos {
namespace MoveBackward {

namespace KE = Kokkos::Experimental;

template <class Tag, class ValueType, class InfoType>
void run_single_scenario(const InfoType& scenario_info, int apiId) {
  const std::size_t view_ext = std::get<1>(scenario_info);

  auto v = create_view<ValueType>(Tag{}, view_ext, "v");

  // v might not be deep copyable so to modify it on the host
  // need to do all this
  auto v_dc   = create_deep_copyable_compatible_view_with_same_extent(v);
  auto v_dc_h = create_mirror_view(Kokkos::HostSpace(), v_dc);
  Kokkos::Random_XorShift64_Pool<Kokkos::DefaultHostExecutionSpace> pool(12371);
  Kokkos::fill_random(v_dc_h, pool, 0, 523);
  // copy to v_dc and then to v
  Kokkos::deep_copy(v_dc, v_dc_h);
  CopyFunctor<decltype(v_dc), decltype(v)> F1(v_dc, v);
  Kokkos::parallel_for("copy", v.extent(0), F1);

  // make a gold copy of v before calling the algorithm
  // since the algorithm will modify v
  auto gold = create_host_space_copy(v);

  // create another view that is bigger than v
  // because we need it to test the move_backward
  auto v2 = create_view<ValueType>(Tag{}, view_ext + 5, "v2");

  if (apiId == 0) {
    auto rit =
        KE::move_backward(exespace(), KE::begin(v), KE::end(v), KE::end(v2));
    const int dist = KE::distance(KE::begin(v2), rit);
    EXPECT_EQ(dist, 5);
  } else if (apiId == 1) {
    auto rit       = KE::move_backward("mylabel", exespace(), KE::begin(v),
                                 KE::end(v), KE::end(v2));
    const int dist = KE::distance(KE::begin(v2), rit);
    EXPECT_EQ(dist, 5);
  } else if (apiId == 2) {
    auto rit       = KE::move_backward(exespace(), v, v2);
    const int dist = KE::distance(KE::begin(v2), rit);
    EXPECT_EQ(dist, 5);
  } else if (apiId == 3) {
    auto rit       = KE::move_backward("mylabel", exespace(), v, v2);
    const int dist = KE::distance(KE::begin(v2), rit);
    EXPECT_EQ(dist, 5);
  }

  // check
  auto v2_h = create_host_space_copy(v2);
  for (std::size_t j = 0; j < v2_h.extent(1); ++j) {
    if (j < 5) {
      EXPECT_TRUE(v2_h(j) == static_cast<ValueType>(0));
    } else {
      EXPECT_TRUE(gold(j - 5) == v2_h(j));
    }
  }
}

template <class Tag, class ValueType>
void run_all_scenarios() {
  const std::map<std::string, std::size_t> scenarios = {
      {"empty", 0},          {"one-element-a", 1},  {"one-element-b", 1},
      {"two-elements-a", 2}, {"two-elements-b", 2}, {"small-a", 9},
      {"small-b", 13},       {"medium", 1103},      {"large", 101513}};

  for (const auto& it : scenarios) {
    run_single_scenario<Tag, ValueType>(it, 0);
    run_single_scenario<Tag, ValueType>(it, 1);
    run_single_scenario<Tag, ValueType>(it, 2);
    run_single_scenario<Tag, ValueType>(it, 3);
  }
}

TEST(std_algorithms_mod_seq_ops, move_backward) {
  run_all_scenarios<DynamicTag, int>();
  run_all_scenarios<DynamicTag, double>();
  run_all_scenarios<StridedThreeTag, int>();
  run_all_scenarios<StridedThreeTag, double>();
}

}  // namespace MoveBackward
}  // namespace stdalgos
}  // namespace Test
