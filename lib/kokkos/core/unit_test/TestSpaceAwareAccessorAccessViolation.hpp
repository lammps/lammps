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

template <class MemorySpace, class ExecutionSpace>
struct TestMemoryAccessViolation {
  Kokkos::Impl::SpaceAwareAccessor<MemorySpace, Kokkos::default_accessor<int>>
      acc;

  KOKKOS_FUNCTION decltype(auto) bad_access() const {
    return acc.access(nullptr, 0);
  }

  KOKKOS_FUNCTION void operator()(int) const { ++bad_access(); }

  TestMemoryAccessViolation(ExecutionSpace const& s,
                            std::string const& matcher) {
    constexpr bool accessible_from_execution_space = Kokkos::SpaceAccessibility<
        /*AccessSpace=*/ExecutionSpace,
        /*MemorySpace=*/MemorySpace>::accessible;
    EXPECT_FALSE(accessible_from_execution_space);
    EXPECT_DEATH(
        {
          Kokkos::parallel_for(Kokkos::RangePolicy<ExecutionSpace>(s, 0, 1),
                               *this);
          Kokkos::fence();
        },
        matcher);
  }
};

template <class MemorySpace, class ExecutionSpace>
void test_memory_access_violation(ExecutionSpace const& s,
                                  std::string const& m) {
  TestMemoryAccessViolation<MemorySpace, ExecutionSpace>(s, m);
}

template <class ExecutionSpace>
void test_memory_access_violations_from_host() {
  using memory_space_t = typename ExecutionSpace::memory_space;
  using exec_space_t   = Kokkos::DefaultHostExecutionSpace;
  const exec_space_t exec_space{};
  std::string const message =
      "Kokkos::SpaceAwareAccessor ERROR: attempt to access inaccessible memory "
      "space";
  test_memory_access_violation<memory_space_t, exec_space_t>(exec_space,
                                                             message);
}

template <class ExecutionSpace>
void test_memory_access_violations_from_device() {
  using memory_space_t = Kokkos::HostSpace;
  using exec_space_t   = ExecutionSpace;
  const exec_space_t exec_space{};
  std::string const message =
      "Kokkos::SpaceAwareAccessor ERROR: attempt to access inaccessible memory "
      "space";
  test_memory_access_violation<memory_space_t, exec_space_t>(exec_space,
                                                             message);
}

// FIXME_SYCL
#if !(defined(KOKKOS_COMPILER_INTEL_LLVM) && defined(KOKKOS_ENABLE_SYCL))
TEST(TEST_CATEGORY_DEATH,
     mdspan_space_aware_accessor_invalid_access_from_host) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  using ExecutionSpace = TEST_EXECSPACE;

  if (Kokkos::SpaceAccessibility<
          /*AccessSpace=*/Kokkos::HostSpace,
          /*MemorySpace=*/typename ExecutionSpace::memory_space>::accessible) {
    GTEST_SKIP() << "skipping since no memory access violation would occur";
  }

  test_memory_access_violations_from_host<ExecutionSpace>();
}
#endif

TEST(TEST_CATEGORY_DEATH,
     mdspan_space_aware_accessor_invalid_access_from_device) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  using ExecutionSpace = TEST_EXECSPACE;

  if (Kokkos::SpaceAccessibility<
          /*AccessSpace=*/ExecutionSpace,
          /*MemorySpace=*/Kokkos::HostSpace>::accessible) {
    GTEST_SKIP() << "skipping since no memory access violation would occur";
  }

#if defined(KOKKOS_ENABLE_SYCL) && defined(NDEBUG)  // FIXME_SYCL
  if (std::is_same<ExecutionSpace, Kokkos::Experimental::SYCL>::value) {
    GTEST_SKIP() << "skipping SYCL device-side abort does not work when NDEBUG "
                    "is defined";
  }
#endif
#if defined(KOKKOS_ENABLE_OPENMPTARGET)  // FIXME_OPENMPTARGET
  if (std::is_same<ExecutionSpace, Kokkos::Experimental::OpenMPTarget>::value) {
    GTEST_SKIP() << "skipping because OpenMPTarget backend is currently not "
                    "able to abort from the device";
  }
#endif
#if defined(KOKKOS_ENABLE_OPENACC)  // FIXME_OPENACC
  if (std::is_same<ExecutionSpace, Kokkos::Experimental::OpenACC>::value) {
    GTEST_SKIP() << "skipping because OpenACC backend is currently not "
                    "able to abort from the device";
  }
#endif

  test_memory_access_violations_from_device<ExecutionSpace>();
}
