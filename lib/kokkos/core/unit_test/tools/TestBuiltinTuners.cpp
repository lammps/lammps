// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif
#include <impl/Kokkos_Tools_Generic.hpp>

#include <gtest/gtest.h>

namespace {

using ExecSpace  = Kokkos::DefaultHostExecutionSpace;
using TeamMember = Kokkos::TeamPolicy<ExecSpace>::member_type;

struct TestTeamFunctor {
  KOKKOS_FUNCTION void operator()(TeamMember) const {}
};

struct TestMDFunctor {
  KOKKOS_FUNCTION void operator()(const int, const int) const {}
};

TEST(kokkosp, builtin_tuner) {
  Kokkos::TeamPolicy<ExecSpace> teamp(1, Kokkos::AUTO, Kokkos::AUTO);
  Kokkos::MDRangePolicy<Kokkos::Rank<2>> mdp({0, 0}, {1, 1});
  Kokkos::Tools::Experimental::TeamSizeTuner team_tune_this(
      "team_tuner", teamp, TestTeamFunctor{}, Kokkos::ParallelForTag{},
      Kokkos::Tools::Experimental::Impl::Impl::SimpleTeamSizeCalculator{});

  Kokkos::Tools::Experimental::MDRangeTuner<2> md_tune_this(
      "md_tuner", mdp, TestMDFunctor{}, Kokkos::ParallelForTag{},
      Kokkos::Tools::Experimental::Impl::Impl::SimpleTeamSizeCalculator{});

  std::vector<int> options{1, 2, 3, 4, 5};

  auto new_team_tuner = team_tune_this.combine("options", options);
  auto new_md_tuner   = md_tune_this.combine("options", options);
  using namespace Kokkos::Tools::Experimental;
  VariableInfo info;
  info.category      = StatisticalCategory::kokkos_value_categorical;
  info.valueQuantity = CandidateValueType::kokkos_value_unbounded;
  info.type          = ValueType::kokkos_value_string;
  size_t input       = declare_input_type("kernel", info);
  VariableValue team_kernel_value = make_variable_value(input, "abs");
  VariableValue md_kernel_value   = make_variable_value(input, "abs");
  size_t kernel_context           = get_new_context_id();
  begin_context(kernel_context);
  set_input_values(kernel_context, 1, &team_kernel_value);
  for (int x = 0; x < 10000; ++x) {
    auto config = new_md_tuner.begin();
    int option  = std::get<0>(config);
    (void)option;
    int tile_x = std::get<1>(config);
    int tile_y = std::get<2>(config);
    Kokkos::parallel_for("mdrange",
                         Kokkos::MDRangePolicy<Kokkos::Rank<2>>(
                             {0, 0}, {1, 1}, {tile_x, tile_y}),
                         TestMDFunctor{});
    new_md_tuner.end();
  }
  end_context(kernel_context);
  begin_context(kernel_context);
  set_input_values(kernel_context, 1, &md_kernel_value);

  /**
   * Note that 0.0 is basically a floating point index into
   * the outermost index in this, which is the options vector
   * above. The At 0.0, this will be the first element (1).
   * At 0.9 this will be the last element (5)
   */
  auto begin_point = new_team_tuner.get_point(0.0, 0.0, 0.0);
  EXPECT_EQ(std::get<0>(begin_point), 1);
  auto end_point = new_team_tuner.get_point(0.9, 0.0, 0.0);
  EXPECT_EQ(std::get<0>(end_point), 5);
  for (int x = 0; x < 10000; ++x) {
    auto config                 = new_team_tuner.begin();
    [[maybe_unused]] int option = std::get<0>(config);
    int team                    = std::get<1>(config);
    int vector                  = std::get<2>(config);
    Kokkos::parallel_for("mdrange",
                         Kokkos::TeamPolicy<ExecSpace>(1, team, vector),
                         TestTeamFunctor{});
    new_team_tuner.end();
  }
  end_context(kernel_context);
}

// Last level of recursion, validates all possible sizes for the last dimension
void validate_tile_sizes(const std::vector<int>& cont,
                         std::array<int, 6> current_tile,
                         const std::array<int, 3>& hw_tile_limits,
                         const int current_rank, const int policy_rank) {
  namespace KTE = Kokkos::Tools::Experimental;
  for (const auto& size : cont) {
    current_tile[current_rank] = size;
    EXPECT_TRUE(
        KTE::Impl::is_valid_tile(hw_tile_limits, current_tile, policy_rank));
  }
}

// Builds candidate tiles by traversing the nested map structure. Each recursive
// level sets one dimension.
template <typename Mapped>
void validate_tile_sizes(const std::map<int, Mapped>& cont,
                         std::array<int, 6> current_tile,
                         const std::array<int, 3>& hw_tile_limits,
                         const int current_rank, const int policy_rank) {
  for (const auto& [dimension_size, sub_cont] : cont) {
    current_tile[current_rank] = dimension_size;
    validate_tile_sizes(sub_cont, current_tile, hw_tile_limits,
                        current_rank + 1, policy_rank);
  }
}

// Verifies that all tiles in the configuration space satisfy hardware
// constraints.
template <typename Mapped>
void validate_tile_configurations(std::map<int, Mapped>& cont,
                                  const std::array<int, 3>& hw_tile_limits,
                                  const int policy_rank) {
  std::array<int, 6> current_tile{1, 1, 1, 1, 1, 1};
  validate_tile_sizes(cont, current_tile, hw_tile_limits, 0, policy_rank);
}

// Test for a given rank that the tile constraints are correctly applied
template <int test_policy_rank>
void test_tile_constraints() {
  namespace KTE          = Kokkos::Tools::Experimental;
  using SpaceDescription = typename KTE::Impl::n_dimensional_sparse_structure<
      int, test_policy_rank>::type;

  std::array<int, 3> hw_thread_limits{128, 128, 128};
  SpaceDescription tile_configuration_space;
  const int max_tile_product = 512;

  KTE::Impl::fill_tile(tile_configuration_space, max_tile_product);
  KTE::Impl::apply_tiles_constraints(tile_configuration_space, hw_thread_limits,
                                     test_policy_rank);
  validate_tile_configurations(tile_configuration_space, hw_thread_limits,
                               test_policy_rank);
}

TEST(kokkosp, constraint_tiles) {
  namespace KTE = Kokkos::Tools::Experimental;

  std::array<int, 3> hw_thread_limits{128, 128, 128};

  std::array<int, 6> tile_6d_within_limits{4, 2, 4, 4, 2, 4};
  std::array<int, 6> tile_6d_large_product{1, 1, 1, 1, 32, 32};
  std::array<int, 6> tile_3d_hw_max{128, 128, 128, 1, 1, 1};
  std::array<int, 6> tile_3d_with_zero{128, 0, 1, 1, 1, 1};
  std::array<int, 6> tile_3d_exceeds_limit{1, 129, 1, 1, 1, 1};

  // KTE::Impl::is_valid_tile checks each hardware dimension but not the total
  // product
  EXPECT_TRUE(
      KTE::Impl::is_valid_tile(hw_thread_limits, tile_6d_within_limits, 6));
  EXPECT_FALSE(
      KTE::Impl::is_valid_tile(hw_thread_limits, tile_6d_large_product, 6));

  EXPECT_TRUE(KTE::Impl::is_valid_tile(hw_thread_limits, tile_3d_hw_max, 3));
  EXPECT_FALSE(
      KTE::Impl::is_valid_tile(hw_thread_limits, tile_3d_with_zero, 3));
  EXPECT_FALSE(
      KTE::Impl::is_valid_tile(hw_thread_limits, tile_3d_exceeds_limit, 3));

  test_tile_constraints<2>();
  test_tile_constraints<3>();
  test_tile_constraints<4>();
  test_tile_constraints<5>();
  test_tile_constraints<6>();
}

}  // namespace
