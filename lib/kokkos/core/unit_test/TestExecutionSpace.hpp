// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <gtest/gtest.h>

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif

#include "tools/include/ToolTestingUtilities.hpp"

namespace {

template <class ExecutionSpace>
struct CheckClassWithExecutionSpaceAsDataMemberIsCopyable {
  Kokkos::DefaultExecutionSpace device;
  Kokkos::DefaultHostExecutionSpace host;

  KOKKOS_FUNCTION void operator()(int i, int& e) const { e += i; }

  CheckClassWithExecutionSpaceAsDataMemberIsCopyable() {
    int errors;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecutionSpace>(0, 1), *this,
                            errors);
    EXPECT_EQ(errors, 0);
  }
};

TEST(TEST_CATEGORY, execution_space_as_class_data_member) {
  CheckClassWithExecutionSpaceAsDataMemberIsCopyable<TEST_EXECSPACE>();
}

TEST(TEST_CATEGORY, execution_space_moved_from) {
  TEST_EXECSPACE exec;
  TEST_EXECSPACE other = std::move(exec);
  // NOLINTNEXTLINE(bugprone-use-after-move)
  ASSERT_EQ(other, exec);
  exec = std::move(other);
  // NOLINTNEXTLINE(bugprone-use-after-move)
  ASSERT_EQ(exec, other);
}

constexpr bool test_execspace_explicit_construction() {
#ifdef KOKKOS_ENABLE_SERIAL
  static_assert(!std::is_convertible_v<Kokkos::NewInstance, Kokkos::Serial>);
#endif
#ifdef KOKKOS_ENABLE_OPENMP
  static_assert(!std::is_convertible_v<int, Kokkos::OpenMP>);
#endif
#ifdef KOKKOS_ENABLE_CUDA
  static_assert(!std::is_convertible_v<cudaStream_t, Kokkos::Cuda>);
#endif
#ifdef KOKKOS_ENABLE_HIP
  static_assert(!std::is_convertible_v<hipStream_t, Kokkos::HIP>);
#endif
#ifdef KOKKOS_ENABLE_HPX
  static_assert(!std::is_convertible_v<Kokkos::Experimental::HPX::instance_mode,
                                       Kokkos::Experimental::HPX>);
  static_assert(!std::is_convertible_v<
                hpx::execution::experimental::unique_any_sender<>&&,
                Kokkos::Experimental::HPX>);
#endif
#ifdef KOKKOS_ENABLE_OPENACC
  static_assert(!std::is_convertible_v<int, Kokkos::Experimental::OpenACC>);
#endif
#ifdef KOKKOS_ENABLE_SYCL
  static_assert(!std::is_convertible_v<sycl::queue, Kokkos::SYCL>);
#endif

  return true;
}

static_assert(test_execspace_explicit_construction());

consteval bool test_execspace_nothrow_copy_and_move() {
  static_assert(std::is_nothrow_copy_constructible_v<TEST_EXECSPACE>);
  static_assert(std::is_nothrow_copy_assignable_v<TEST_EXECSPACE>);
  static_assert(std::is_nothrow_move_constructible_v<TEST_EXECSPACE>);
  static_assert(std::is_nothrow_move_assignable_v<TEST_EXECSPACE>);
  return true;
}

static_assert(test_execspace_nothrow_copy_and_move());

// We don't actually promise a tool-observable event and acknowledge that some
// backend mights not need a fence to ensure that all enqueued work has finished
// before an execution space instance is destroyed. Therefore we might want to
// revisit this test.
TEST(TEST_CATEGORY, execution_space_fence_on_destruction) {
  auto [dummy_instance] =
      Kokkos::Experimental::partition_space(TEST_EXECSPACE(), 1);
  bool created_new_instance = TEST_EXECSPACE() != dummy_instance;
  if (!created_new_instance)
    GTEST_SKIP() << "partition_space doesn't create a new instance";

  Kokkos::Test::Tools::listen_tool_events(
      Kokkos::Test::Tools::Config::DisableAll(),
      Kokkos::Test::Tools::Config::EnableFences());

  auto success = Kokkos::Test::Tools::validate_existence(
      [&]() {
        [[maybe_unused]] auto [new_instance] =
            Kokkos::Experimental::partition_space(TEST_EXECSPACE(), 1);
      },
      [&](Kokkos::Test::Tools::BeginFenceEvent event) {
        return Kokkos::Test::Tools::MatchDiagnostic{
            event.descriptor().find("fence on destruction") !=
            std::string::npos};
      });
  ASSERT_TRUE(success);
  Kokkos::Test::Tools::listen_tool_events(
      Kokkos::Test::Tools::Config::DisableAll());
}

}  // namespace
