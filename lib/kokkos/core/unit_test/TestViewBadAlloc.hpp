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

#include <Kokkos_Core.hpp>

#include <gtest/gtest.h>

namespace {

template <class MemorySpace>
void test_view_bad_alloc() {
  bool did_throw    = false;
  auto too_large    = std::numeric_limits<size_t>::max() - 42;
  std::string label = "my_label";
  try {
    auto should_always_fail =
        Kokkos::View<double *, MemorySpace>(label, too_large);
  } catch (std::runtime_error const &error) {
    std::string msg = error.what();
    ASSERT_PRED_FORMAT2(
        ::testing::IsSubstring,
        std::string(MemorySpace::name()) + " memory space failed to allocate",
        msg)
        << "memory space name is missing";
    ASSERT_PRED_FORMAT2(::testing::IsSubstring,
                        std::string("(label=\"") + label + "\")", msg)
        << "label is missing";
    did_throw = true;
  }
  ASSERT_TRUE(did_throw);
}

TEST(TEST_CATEGORY, view_bad_alloc) {
  using ExecutionSpace = TEST_EXECSPACE;
  using MemorySpace    = ExecutionSpace::memory_space;
#if defined(__has_feature)
#if __has_feature(address_sanitizer)
  if (std::is_same_v<MemorySpace, Kokkos::HostSpace>) {
    GTEST_SKIP() << "AddressSanitizer detects allocating too much memory "
                    "preventing our checks to run";
  }
#endif
#endif
#if ((HIP_VERSION_MAJOR == 5) && (HIP_VERSION_MINOR == 3))
  if (std::is_same_v<ExecutionSpace, Kokkos::HIP>) {
    GTEST_SKIP()
        << "ROCm 5.3 segfaults when trying to allocate too much memory";
  }
#endif
#if defined(KOKKOS_ENABLE_OPENACC)  // FIXME_OPENACC
  if (std::is_same_v<ExecutionSpace, Kokkos::Experimental::OpenACC>) {
    GTEST_SKIP() << "acc_malloc() not properly returning nullptr";
  }
#endif

  test_view_bad_alloc<MemorySpace>();

  constexpr bool execution_space_is_device =
      std::is_same_v<ExecutionSpace, Kokkos::DefaultExecutionSpace> &&
      !std::is_same_v<Kokkos::DefaultExecutionSpace,
                      Kokkos::DefaultHostExecutionSpace>;

  if constexpr (execution_space_is_device) {
    if constexpr (Kokkos::has_shared_space) {
      test_view_bad_alloc<Kokkos::SharedSpace>();
    }
    if constexpr (Kokkos::has_shared_host_pinned_space) {
      test_view_bad_alloc<Kokkos::SharedHostPinnedSpace>();
    }
  }
}

}  // namespace
