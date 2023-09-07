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

namespace {

template <typename ViewType>
struct Increment {
  ViewType view_;

  template <typename ExecutionSpace>
  explicit Increment(ExecutionSpace, ViewType view) : view_(view) {
    Kokkos::parallel_for(
        "increment",
        Kokkos::RangePolicy<ExecutionSpace, Kokkos::IndexType<size_t>>{
            0, view_.size()},
        *this);
  }

  KOKKOS_FUNCTION
  void operator()(const size_t idx) const { ++view_(idx); }
};

template <typename ViewType>
struct CheckResult {
  ViewType view_;
  int targetVal_;
  unsigned numErrors = 0;

  template <typename ExecutionSpace>
  CheckResult(ExecutionSpace, ViewType view, int targetVal)
      : view_(view), targetVal_(targetVal) {
    Kokkos::parallel_reduce(
        "check",
        Kokkos::RangePolicy<ExecutionSpace, Kokkos::IndexType<size_t>>{
            0, view_.size()},
        *this, Kokkos::Sum<unsigned>(numErrors));
  }

  KOKKOS_FUNCTION
  void operator()(const size_t idx, unsigned& errors) const {
    if (view_(idx) != targetVal_) ++errors;
  }
};

TEST(defaultdevicetype, shared_host_pinned_space) {
  ASSERT_TRUE(Kokkos::has_shared_host_pinned_space);

  if constexpr (std::is_same_v<Kokkos::DefaultExecutionSpace,
                               Kokkos::DefaultHostExecutionSpace>)
    GTEST_SKIP() << "Skipping as host and device are the same space";

  const unsigned int numDeviceHostCycles = 3;
  size_t numInts                         = 1024;

  using DeviceExecutionSpace = Kokkos::DefaultExecutionSpace;
  using HostExecutionSpace   = Kokkos::DefaultHostExecutionSpace;

  // ALLOCATION
  Kokkos::View<int*, Kokkos::SharedHostPinnedSpace> sharedData("sharedData",
                                                               numInts);
  // MAIN LOOP
  unsigned incrementCount = 0;

  for (unsigned i = 0; i < numDeviceHostCycles; ++i) {
    // INCREMENT DEVICE
    Increment incrementOnDevice(DeviceExecutionSpace{}, sharedData);
    ++incrementCount;
    Kokkos::fence();
    // CHECK RESULTS HOST
    ASSERT_EQ(
        CheckResult(HostExecutionSpace{}, sharedData, incrementCount).numErrors,
        0u)
        << "Changes to SharedHostPinnedSpace made on device not visible to "
           "host. Iteration "
        << i << " of " << numDeviceHostCycles;

    // INCREMENT HOST
    Increment incrementOnHost(HostExecutionSpace{}, sharedData);
    ++incrementCount;
    Kokkos::fence();
    // CHECK RESULTS Device
    ASSERT_EQ(CheckResult(DeviceExecutionSpace{}, sharedData, incrementCount)
                  .numErrors,
              0u)
        << "Changes to SharedHostPinnedSpace made on host not visible to "
           "device. Iteration "
        << i << " of " << numDeviceHostCycles;
  }
}
}  // namespace
