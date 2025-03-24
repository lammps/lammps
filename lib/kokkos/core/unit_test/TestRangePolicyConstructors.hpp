//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

#include <regex>
#include <limits>
#include <type_traits>

namespace {

TEST(TEST_CATEGORY, range_policy_runtime_parameters) {
  using Policy     = Kokkos::RangePolicy<>;
  using Index      = Policy::index_type;
  Index work_begin = 5;
  Index work_end   = 15;
  Index chunk_size = 10;
  {
    Policy p(work_begin, work_end);
    ASSERT_EQ(p.begin(), work_begin);
    ASSERT_EQ(p.end(), work_end);
  }
  {
    Policy p(Kokkos::DefaultExecutionSpace(), work_begin, work_end);
    ASSERT_EQ(p.begin(), work_begin);
    ASSERT_EQ(p.end(), work_end);
  }
  {
    Policy p(work_begin, work_end, Kokkos::ChunkSize(chunk_size));
    ASSERT_EQ(p.begin(), work_begin);
    ASSERT_EQ(p.end(), work_end);
    ASSERT_EQ(p.chunk_size(), chunk_size);
  }
  {
    Policy p(Kokkos::DefaultExecutionSpace(), work_begin, work_end,
             Kokkos::ChunkSize(chunk_size));
    ASSERT_EQ(p.begin(), work_begin);
    ASSERT_EQ(p.end(), work_end);
    ASSERT_EQ(p.chunk_size(), chunk_size);
  }
  {
    Policy p;  // default-constructed
    ASSERT_EQ(p.begin(), Index(0));
    ASSERT_EQ(p.end(), Index(0));
    ASSERT_EQ(p.chunk_size(), Index(0));

    // copy-assigned
    p = Policy(work_begin, work_end, Kokkos::ChunkSize(chunk_size));
    ASSERT_EQ(p.begin(), work_begin);
    ASSERT_EQ(p.end(), work_end);
    ASSERT_EQ(p.chunk_size(), chunk_size);
  }
  {
    Policy p1(work_begin, work_end, Kokkos::ChunkSize(chunk_size));
    Policy p2(p1);  // copy-constructed
    ASSERT_EQ(p1.begin(), p2.begin());
    ASSERT_EQ(p1.end(), p2.end());
    ASSERT_EQ(p1.chunk_size(), p2.chunk_size());
  }
}

TEST(TEST_CATEGORY_DEATH, range_policy_invalid_bounds) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  using Policy    = Kokkos::RangePolicy<TEST_EXECSPACE>;
  using ChunkSize = Kokkos::ChunkSize;

  std::string msg =
      "Kokkos::RangePolicy bounds error: The lower bound (100) is greater than "
      "the upper bound (90).\n";
#ifndef KOKKOS_ENABLE_DEPRECATED_CODE_4
  // escape the parentheses in the regex to match the error message
  msg = std::regex_replace(msg, std::regex("\\(|\\)"), "\\$&");
  ASSERT_DEATH({ (void)Policy(100, 90); }, msg);

  ASSERT_DEATH({ (void)Policy(TEST_EXECSPACE(), 100, 90, ChunkSize(10)); },
               msg);
#else

  if (!Kokkos::show_warnings()) {
    GTEST_SKIP() << "Kokkos warning messages are disabled";
  }

  {
    ::testing::internal::CaptureStderr();
    Policy policy(100, 90);
    ASSERT_EQ((int)policy.begin(), 0);
    ASSERT_EQ((int)policy.end(), 0);
#ifdef KOKKOS_ENABLE_DEPRECATION_WARNINGS
    ASSERT_EQ(::testing::internal::GetCapturedStderr(), msg);
#else
    ASSERT_TRUE(::testing::internal::GetCapturedStderr().empty());
    (void)msg;
#endif
  }

  {
    ::testing::internal::CaptureStderr();
    Policy policy(TEST_EXECSPACE(), 100, 90, ChunkSize(10));
    ASSERT_EQ((int)policy.begin(), 0);
    ASSERT_EQ((int)policy.end(), 0);
#ifdef KOKKOS_ENABLE_DEPRECATION_WARNINGS
    ASSERT_EQ(::testing::internal::GetCapturedStderr(), msg);
#else
    ASSERT_TRUE(::testing::internal::GetCapturedStderr().empty());
    (void)msg;
#endif
  }

#endif
}

struct W {  // round-trip conversion check for narrowing should "fire"
  W(int const* ptr) : val_(*ptr) {}
  W(int) : val_(0) {}
  operator int() const { return val_; }

  int val_;
};

TEST(TEST_CATEGORY_DEATH, range_policy_round_trip_conversion_fires) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  using Policy = Kokkos::RangePolicy<>;

  static_assert(std::is_convertible_v<W, Policy::index_type>);
  static_assert(std::is_convertible_v<Policy::index_type, W>);

  int const n = 1;
  [[maybe_unused]] std::string msg =
      "Kokkos::RangePolicy bound type error: an unsafe implicit conversion is "
      "performed";
#ifndef KOKKOS_ENABLE_DEPRECATED_CODE_4
  ASSERT_DEATH((void)Policy(0, W(&n)), msg);
#else
  ::testing::internal::CaptureStderr();
  (void)Policy(0, W(&n));
  auto s = std::string(::testing::internal::GetCapturedStderr());
#ifdef KOKKOS_ENABLE_DEPRECATION_WARNINGS
  if (Kokkos::show_warnings()) {
    ASSERT_NE(s.find(msg), std::string::npos) << msg;
  } else
#endif
    ASSERT_TRUE(s.empty());
#endif
}

struct B {  // round-trip conversion would not compile
  B(int const* ptr) : val_(*ptr) {}
  operator int() const { return val_; }

  int val_;
};

TEST(TEST_CATEGORY, range_policy_one_way_convertible_bounds) {
  using Policy    = Kokkos::RangePolicy<>;
  using IndexType = Policy::index_type;

  static_assert(std::is_convertible_v<B, IndexType>);
  static_assert(!std::is_convertible_v<IndexType, B>);

  int const n = 1;
  Policy policy(0, B(&n));
  EXPECT_EQ(policy.begin(), static_cast<IndexType>(0));
  EXPECT_EQ(policy.end(), static_cast<IndexType>(1));
}

TEST(TEST_CATEGORY_DEATH, range_policy_check_sign_changes) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  using UInt32Policy =
      Kokkos::RangePolicy<TEST_EXECSPACE, Kokkos::IndexType<std::uint32_t>>;

  [[maybe_unused]] std::string msg =
      "Kokkos::RangePolicy bound type error: an unsafe implicit conversion is "
      "performed";
#ifndef KOKKOS_ENABLE_DEPRECATED_CODE_4
  {
    std::int64_t n = std::numeric_limits<std::int64_t>::max();
    ASSERT_DEATH((void)UInt32Policy(0, n), msg);
  }
  {
    std::int64_t n = std::numeric_limits<std::int64_t>::min();
    ASSERT_DEATH((void)UInt32Policy(n, 0), msg);
  }
#else
  {
    ::testing::internal::CaptureStderr();
    std::int64_t n = std::numeric_limits<std::int64_t>::max();
    (void)UInt32Policy(0, n);
    auto s = std::string(::testing::internal::GetCapturedStderr());
#ifdef KOKKOS_ENABLE_DEPRECATION_WARNINGS
    if (Kokkos::show_warnings()) {
      ASSERT_NE(s.find(msg), std::string::npos) << msg;
    }
#endif
  }
  {
    ::testing::internal::CaptureStderr();
    std::int64_t n = std::numeric_limits<std::int64_t>::min();
    (void)UInt32Policy(n, 0);
    auto s = std::string(::testing::internal::GetCapturedStderr());
#ifdef KOKKOS_ENABLE_DEPRECATION_WARNINGS
    if (Kokkos::show_warnings()) {
      ASSERT_NE(s.find(msg), std::string::npos) << msg;
    }
#endif
  }
#endif
}

TEST(TEST_CATEGORY_DEATH, range_policy_implicitly_converted_bounds) {
  using UIntIndexType = Kokkos::IndexType<unsigned>;
  using IntIndexType  = Kokkos::IndexType<int>;
  using UIntPolicy    = Kokkos::RangePolicy<TEST_EXECSPACE, UIntIndexType>;
  using IntPolicy     = Kokkos::RangePolicy<TEST_EXECSPACE, IntIndexType>;

  std::string msg =
      "Kokkos::RangePolicy bound type error: an unsafe implicit conversion is "
      "performed on a bound (), which may not preserve its original value.\n";

  [[maybe_unused]] auto get_error_msg = [](auto str, auto val) {
    return str.insert(str.find("(") + 1, std::to_string(val).c_str());
  };
#ifndef KOKKOS_ENABLE_DEPRECATED_CODE_4
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  std::string expected = std::regex_replace(msg, std::regex("\\(|\\)"), "\\$&");
  {
    int test_val = -1;
    ASSERT_DEATH({ (void)UIntPolicy(test_val, 10); },
                 get_error_msg(expected, test_val));
  }
  {
    unsigned test_val = std::numeric_limits<unsigned>::max();
    ASSERT_DEATH({ (void)IntPolicy(0u, test_val); },
                 get_error_msg(expected, test_val));
  }
  {
    long long test_val = std::numeric_limits<long long>::max();
    ASSERT_DEATH({ (void)IntPolicy(0LL, test_val); },
                 get_error_msg(expected, test_val));
  }
  {
    int test_val = -1;
    ASSERT_DEATH({ (void)UIntPolicy(test_val, 10, Kokkos::ChunkSize(2)); },
                 get_error_msg(expected, test_val));
  }

#else
  {
    ::testing::internal::CaptureStderr();
    int test_val = -1;
    UIntPolicy policy(test_val, 10);
    ASSERT_EQ(policy.begin(), 0u);
    ASSERT_EQ(policy.end(), 0u);
#ifdef KOKKOS_ENABLE_DEPRECATION_WARNINGS
    if (Kokkos::show_warnings()) {
      auto s = std::string(::testing::internal::GetCapturedStderr());
      ASSERT_EQ(s.substr(0, s.find("\n") + 1), get_error_msg(msg, test_val));
    }
#else
    ASSERT_TRUE(::testing::internal::GetCapturedStderr().empty());
    (void)msg;
    (void)get_error_msg;
#endif
  }
  {
    ::testing::internal::CaptureStderr();
    unsigned test_val = std::numeric_limits<unsigned>::max();
    IntPolicy policy(0u, test_val);
    ASSERT_EQ(policy.begin(), 0);
    ASSERT_EQ(policy.end(), 0);
#ifdef KOKKOS_ENABLE_DEPRECATION_WARNINGS
    if (Kokkos::show_warnings()) {
      auto s = std::string(::testing::internal::GetCapturedStderr());
      ASSERT_EQ(s.substr(0, s.find("\n") + 1), get_error_msg(msg, test_val));
    }
#else
    ASSERT_TRUE(::testing::internal::GetCapturedStderr().empty());
    (void)msg;
    (void)get_error_msg;
#endif
  }
#endif
}

constexpr bool test_chunk_size_explicit() {
  using ExecutionSpace = TEST_EXECSPACE;
  using Kokkos::ChunkSize;

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  static_assert(std::is_convertible_v<int, ChunkSize>);
  static_assert(std::is_constructible_v<ChunkSize, int>);
  // Some execution spaces were implicitly constructible from int
  // which made the constructor call ambiguous.
  static_assert(
      std::is_constructible_v<Kokkos::DefaultExecutionSpace, int> ||
      std::is_constructible_v<
          Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>, int, int, int>);
  static_assert(std::is_constructible_v<
                Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>, int, int,
                ChunkSize>);
  static_assert(std::is_constructible_v<Kokkos::RangePolicy<ExecutionSpace>,
                                        ExecutionSpace, int, int, int>);
  static_assert(std::is_constructible_v<Kokkos::RangePolicy<ExecutionSpace>,
                                        ExecutionSpace, int, int, ChunkSize>);
#else
  static_assert(!std::is_convertible_v<int, ChunkSize>);
  static_assert(std::is_constructible_v<ChunkSize, int>);
  static_assert(
      !std::is_constructible_v<
          Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>, int, int, int>);
  static_assert(std::is_constructible_v<
                Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>, int, int,
                ChunkSize>);
  static_assert(!std::is_constructible_v<Kokkos::RangePolicy<ExecutionSpace>,
                                         ExecutionSpace, int, int, int>);
  static_assert(std::is_constructible_v<Kokkos::RangePolicy<ExecutionSpace>,
                                        ExecutionSpace, int, int, ChunkSize>);
#endif
  return true;
}

static_assert(test_chunk_size_explicit());

}  // namespace
