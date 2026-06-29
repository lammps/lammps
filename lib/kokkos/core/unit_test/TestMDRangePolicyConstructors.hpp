// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <gtest/gtest.h>

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
import kokkos.core_impl;
#else
#include <Kokkos_Core.hpp>
#endif

#include <regex>

namespace {

template <class IndexType>
void construct_mdrange_policy_variable_type() {
  (void)Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>{
      Kokkos::Array<IndexType, 2>{}, Kokkos::Array<IndexType, 2>{}};

  (void)Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>{
      {{IndexType(0), IndexType(0)}}, {{IndexType(2), IndexType(2)}}};

  (void)Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>{
      {IndexType(0), IndexType(0)}, {IndexType(2), IndexType(2)}};
}

TEST(TEST_CATEGORY, md_range_policy_construction_from_arrays) {
  {
    // Check that construction from Kokkos::Array of the specified index type
    // works.
    using IndexType = unsigned long long;
    Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>,
                          Kokkos::IndexType<IndexType>>
        p1(Kokkos::Array<IndexType, 2>{{0, 1}},
           Kokkos::Array<IndexType, 2>{{2, 3}});
    Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>,
                          Kokkos::IndexType<IndexType>>
        p2(Kokkos::Array<IndexType, 2>{{0, 1}},
           Kokkos::Array<IndexType, 2>{{2, 3}});
    Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>,
                          Kokkos::IndexType<IndexType>>
        p3(Kokkos::Array<IndexType, 2>{{0, 1}},
           Kokkos::Array<IndexType, 2>{{2, 3}},
           Kokkos::Array<IndexType, 1>{{4}});
  }
  {
    // Check that construction from double-braced initializer list
    // works.
    using index_type = unsigned long long;
    Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>> p1({{0, 1}},
                                                              {{2, 3}});
    Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>,
                          Kokkos::IndexType<index_type>>
        p2({{0, 1}}, {{2, 3}});
  }
  {
    // Check that construction from Kokkos::Array of long compiles for backwards
    // compability.  This was broken in
    // https://github.com/kokkos/kokkos/pull/3527/commits/88ea8eec6567c84739d77bdd25fdbc647fae28bb#r512323639
    Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>> p1(
        Kokkos::Array<long, 2>{{0, 1}}, Kokkos::Array<long, 2>{{2, 3}});
    Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>> p2(
        Kokkos::Array<long, 2>{{0, 1}}, Kokkos::Array<long, 2>{{2, 3}});
    Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>> p3(
        Kokkos::Array<long, 2>{{0, 1}}, Kokkos::Array<long, 2>{{2, 3}},
        Kokkos::Array<long, 1>{{4}});
  }

  // Check that construction from various index types works.
  construct_mdrange_policy_variable_type<char>();
  construct_mdrange_policy_variable_type<int>();
  construct_mdrange_policy_variable_type<unsigned long>();
  construct_mdrange_policy_variable_type<std::int64_t>();
}

TEST(TEST_CATEGORY_DEATH, md_range_policy_bounds_unsafe_narrowing_conversions) {
  using Policy = Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>,
                                       Kokkos::IndexType<unsigned>>;

  std::string msg =
      "Kokkos::MDRangePolicy bound type error: an unsafe implicit conversion "
      "is "
      "performed on a bound (-1) in dimension (0), which may not preserve its "
      "original value.\n";
  std::string expected = std::regex_replace(msg, std::regex("\\(|\\)"), "\\$&");

  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  ASSERT_DEATH({ (void)Policy({-1, 0}, {2, 3}); }, expected);
}

TEST(TEST_CATEGORY_DEATH, md_range_policy_invalid_bounds) {
  using Policy = Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>;

  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  auto dim0 = (Policy::inner_direction == Kokkos::Iterate::Right) ? 1 : 0;

  std::string msg1 =
      "Kokkos::MDRangePolicy bounds error: The lower bound (100) is greater "
      "than its upper bound (90) in dimension " +
      std::to_string(dim0) + ".\n";

  // escape the parentheses in the regex to match the error message
  msg1 = std::regex_replace(msg1, std::regex("\\(|\\)"), "\\$&");
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  ASSERT_DEATH({ (void)Policy({100, 100}, {90, 90}); }, msg1);
}

// Verify that we get an error if the user requests tile dimensions too large
// for the specified LaunchBounds.
TEST(TEST_CATEGORY_DEATH, md_range_policy_tile_dims_exceed_launch_bounds) {
#if defined(KOKKOS_ENABLE_CUDA)
  if constexpr (!std::is_same_v<TEST_EXECSPACE, Kokkos::Cuda>) {
#elif defined(KOKKOS_ENABLE_HIP)
  if constexpr (!std::is_same_v<TEST_EXECSPACE, Kokkos::HIP>) {
#endif
    GTEST_SKIP()
        << "LaunchBounds verification only applies to CUDA and HIP backends";
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
  }

  // Check error message when user provided tile dims exceed user specified
  // LaunchBounds.
  std::string msg =
      "Kokkos::MDRangePolicy tile dimensions error: Product of tile "
      "dimensions (256) is greater than the maximum specified via "
      "LaunchBounds (32) - choose smaller tile dims\n";

  std::string expected = std::regex_replace(msg, std::regex("\\(|\\)"), "\\$&");
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  using Policy = Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>,
                                       Kokkos::LaunchBounds<32>>;
  ASSERT_DEATH({ (void)Policy({0, 0}, {128, 128}, {64, 4}); }, expected);
#endif
}

// Test tile size recommendation
template <int Rank, Kokkos::Iterate InnerDirection>
void test_get_tile_size() {
  using Policy =
      Kokkos::MDRangePolicy<TEST_EXECSPACE,
                            Kokkos::Rank<Rank, InnerDirection, InnerDirection>>;

  using tile_type = typename Policy::tile_type;
  using point_t   = typename Policy::point_type;

  point_t lower{};
  point_t upper{};

  const int dim_length = 32;

  for (int i = 0; i < Rank; i++) {
    lower[i] = 0;
    upper[i] = dim_length;
  }

  {
    Policy policy_default(lower, upper);
    auto rec_tile_sizes      = policy_default.tile_size_recommended();
    auto internal_tile_sizes = policy_default.m_tile;

    for (std::size_t i = 0; i < Rank; ++i) {
      EXPECT_EQ(rec_tile_sizes[i], internal_tile_sizes[i])
          << " incorrect recommended tile size returned for rank " << i;
    }
  }

  {
    tile_type user_tile_sizes{};
    for (int i = 0; i < Rank; i++) {
      user_tile_sizes[i] = 2;
    }

    Policy policy(lower, upper, user_tile_sizes);

    auto rec_tile_sizes = policy.tile_size_recommended();

    int prod_rec_tile_size = 1;
    for (std::size_t i = 0; i < Rank; ++i) {
      EXPECT_GT(rec_tile_sizes[i], 0)
          << " invalid default tile size for rank " << i;
      prod_rec_tile_size *= rec_tile_sizes[i];
    }
    EXPECT_LT(prod_rec_tile_size, policy.max_total_tile_size());
  }
}

template <int... Ranks>
void test_get_tile_size_for_ranks(std::integer_sequence<int, Ranks...>) {
  (test_get_tile_size<Ranks, Kokkos::Iterate::Left>(), ...);
  (test_get_tile_size<Ranks, Kokkos::Iterate::Right>(), ...);
}

// Check that tile_size_recommended() returns valid tile sizes consistent with
// internal tile dimensions
TEST(TEST_CATEGORY, md_range_policy_get_tile_size) {
  constexpr auto ranks = std::integer_sequence<int, 2, 3, 4, 5, 6>{};
  test_get_tile_size_for_ranks(ranks);
}

template <int Rank, int MaxTperB, Kokkos::Iterate InnerDirection>
void test_default_tiles_respect_launch_bounds() {
  using policy_t =
      Kokkos::MDRangePolicy<TEST_EXECSPACE,
                            Kokkos::Rank<Rank, InnerDirection, InnerDirection>,
                            Kokkos::LaunchBounds<MaxTperB>>;
  using point_t    = typename policy_t::point_type;
  using index_type = typename policy_t::index_type;

  point_t lower{};
  point_t upper{};
  for (int i = 0; i < Rank; i++) {
    lower[i] = 0;
    upper[i] = 32;
  }

  policy_t policy_with_default_tile(lower, upper);

  EXPECT_LE(policy_with_default_tile.m_prod_tile_dims,
            static_cast<index_type>(MaxTperB))
      << " for Rank-" << Rank << " with LaunchBounds<" << MaxTperB << ">"
      << " and InnerDirection "
      << (InnerDirection == Kokkos::Iterate::Left ? "Left" : "Right");
}

// Expand ranks (2, 3, 4, 5, 6)
template <int MaxTperB, Kokkos::Iterate InnerDirection, int... Ranks>
void test_default_tiles_for_ranks(std::integer_sequence<int, Ranks...>) {
  (test_default_tiles_respect_launch_bounds<Ranks, MaxTperB, InnerDirection>(),
   ...);
}

// Expand launch bounds (256, 128, 64, 32, 16)
template <int... MaxTperBs>
void test_default_tiles_for_all_configs(
    std::integer_sequence<int, MaxTperBs...>) {
  constexpr auto ranks = std::integer_sequence<int, 2, 3, 4, 5, 6>{};
  (test_default_tiles_for_ranks<MaxTperBs, Kokkos::Iterate::Left>(ranks), ...);
  (test_default_tiles_for_ranks<MaxTperBs, Kokkos::Iterate::Right>(ranks), ...);
}

// Check that MDRangePolicy auto-computed tile sizes never exceed user-specified
// LaunchBounds
TEST(TEST_CATEGORY, md_range_policy_default_tiles_respect_launch_bounds) {
#if defined(KOKKOS_ENABLE_CUDA)
  if constexpr (!std::is_same_v<TEST_EXECSPACE, Kokkos::Cuda>) {
#elif defined(KOKKOS_ENABLE_HIP)
    if constexpr (!std::is_same_v<TEST_EXECSPACE, Kokkos::HIP>) {
#endif
    GTEST_SKIP()
        << "LaunchBounds verification only applies to CUDA and HIP backends";
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
  }
  // Verify that auto-computed tiles never exceed LaunchBounds.
  test_default_tiles_for_all_configs(
      std::integer_sequence<int, 256, 128, 64, 32, 16>{});
#endif
}

// The execution space is defaulted if not given to the constructor.
TEST(TEST_CATEGORY, md_range_policy_default_space) {
  using policy_t = Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>;

  policy_t defaulted({42, 47}, {666, 999});

  ASSERT_EQ(defaulted.space(), TEST_EXECSPACE{});
}

// The execution space instance can be updated.
TEST(TEST_CATEGORY, md_range_policy_impl_set_space) {
  using policy_t = Kokkos::MDRangePolicy<TEST_EXECSPACE, Kokkos::Rank<2>>;

  const auto [exec_old, exec_new] =
      Kokkos::Experimental::partition_space(TEST_EXECSPACE{}, 1, 1);

  const policy_t policy_old(exec_old, {42, 47}, {666, 999});
  ASSERT_EQ(policy_old.space(), exec_old);

  const policy_t policy_new(Kokkos::Impl::PolicyUpdate{}, policy_old, exec_new);
  ASSERT_EQ(policy_new.space(), exec_new);
  ASSERT_EQ(policy_new.m_lower, (typename policy_t::point_type{42, 47}));
  ASSERT_EQ(policy_new.m_upper, (typename policy_t::point_type{666, 999}));
}

}  // namespace
