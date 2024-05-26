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

#if defined(_WIN32)
#include <windows.h>
unsigned getBytesPerPage() {
  SYSTEM_INFO si;
  GetSystemInfo(&si);
  return si.dwPageSize;
}

#else  // unix/posix system
#include <unistd.h>
unsigned getBytesPerPage() { return sysconf(_SC_PAGESIZE); }
#endif

#include <algorithm>
#include <numeric>
#include <iostream>

namespace {
void printTimings(std::ostream& out, std::vector<uint64_t> const& tr,
                  uint64_t threshold = (std::numeric_limits<uint64_t>::max)()) {
  out << "TimingResult contains " << tr.size() << " results:\n";
  for (auto it = tr.begin(); it != tr.end(); ++it) {
    out << "Duration of loop " << it - tr.begin() << " is " << *it
        << " clock cycles. ";
    if ((*it) > threshold) out << "Migration assumed.";

    out << "\n";
  }
}

template <typename T>
T computeMean(std::vector<T> const& results) {
  return std::accumulate(results.begin(), results.end(), T{}) / results.size();
}

template <typename ViewType>
class IncrementFunctor {
 private:
  using index_type = decltype(std::declval<ViewType>().size());
  ViewType view_;

 public:
  IncrementFunctor() = delete;

  explicit IncrementFunctor(ViewType view) : view_(view) {}

  KOKKOS_INLINE_FUNCTION
  void operator()(const index_type idx, uint64_t& clockTics) const {
    uint64_t start = Kokkos::Impl::clock_tic();
    ++view_(idx);
    clockTics += Kokkos::Impl::clock_tic() - start;
  }
};

// TIMING CAPTURED KERNEL
// PREMISE: This kernel should always be memory bound, as we are measuring
// memory access times. The compute load of an increment is small enough on
// current hardware but this could be different for new hardware. As we count
// the clocks in the kernel, the core frequency of the device has to be fast
// enough to guarante that the kernel stays memory bound.
template <typename ExecSpace, typename ViewType>
std::vector<uint64_t> incrementInLoop(ViewType& view,
                                      unsigned int numRepetitions) {
  using index_type = decltype(view.size());
  std::vector<uint64_t> results;

  Kokkos::fence();
  for (unsigned i = 0; i < numRepetitions; ++i) {
    uint64_t sum_clockTics;
    IncrementFunctor<ViewType> func(view);
    Kokkos::parallel_reduce(
        "increment",
        Kokkos::RangePolicy<ExecSpace, Kokkos::IndexType<index_type>>{
            0, view.size()},
        func, sum_clockTics);
    Kokkos::fence();
    results.push_back(sum_clockTics / view.size());
  }
  return results;
}

TEST(defaultdevicetype, shared_space) {
  ASSERT_TRUE(Kokkos::has_shared_space);

  if constexpr (std::is_same_v<Kokkos::DefaultExecutionSpace,
                               Kokkos::DefaultHostExecutionSpace>)
    GTEST_SKIP() << "Skipping as host and device are the same space";

#if defined(KOKKOS_ARCH_AMD_GPU) && defined(KOKKOS_ENABLE_HIP)
  if (!Kokkos::SharedSpace().impl_hip_driver_check_page_migration())
    GTEST_SKIP()
        << "skipping because specified arch does not support page migration";
#endif
#if defined(KOKKOS_ENABLE_SYCL) && !defined(KOKKOS_ARCH_INTEL_GPU)
  GTEST_SKIP()
      << "skipping because clock_tic is only defined for sycl+intel gpu";
#endif

  const unsigned int numRepetitions      = 10;
  const unsigned int numDeviceHostCycles = 3;
  double threshold                       = 1.5;
  unsigned int numPages                  = 100;
  size_t numBytes                        = numPages * getBytesPerPage();

  using DeviceExecutionSpace = Kokkos::DefaultExecutionSpace;
  using HostExecutionSpace   = Kokkos::DefaultHostExecutionSpace;

  // ALLOCATION
  Kokkos::View<int*, Kokkos::SharedSpace> sharedData("sharedData",
                                                     numBytes / sizeof(int));
  Kokkos::View<int*, DeviceExecutionSpace::memory_space> deviceData(
      "deviceData", numBytes / sizeof(int));
  Kokkos::View<int*, HostExecutionSpace::memory_space> hostData(
      "hostData", numBytes / sizeof(int));
  Kokkos::fence();

  // GET DEFAULT EXECSPACE LOCAL TIMINGS
  auto deviceLocalTimings =
      incrementInLoop<DeviceExecutionSpace>(deviceData, numRepetitions);

  // GET DEFAULT HOSTEXECSPACE LOCAL TIMINGS
  auto hostLocalTimings =
      incrementInLoop<HostExecutionSpace>(hostData, numRepetitions);

  // GET PAGE MIGRATING TIMINGS DATA
  std::vector<decltype(deviceLocalTimings)> deviceSharedTimings{};
  std::vector<decltype(hostLocalTimings)> hostSharedTimings{};
  for (unsigned i = 0; i < numDeviceHostCycles; ++i) {
    // GET RESULTS DEVICE
    deviceSharedTimings.push_back(
        incrementInLoop<DeviceExecutionSpace>(sharedData, numRepetitions));

    // GET RESULTS HOST
    hostSharedTimings.push_back(
        incrementInLoop<HostExecutionSpace>(sharedData, numRepetitions));
  }

  // COMPUTE STATISTICS OF HOST AND DEVICE LOCAL KERNELS
  auto deviceLocalMean = computeMean(deviceLocalTimings);
  auto hostLocalMean   = computeMean(hostLocalTimings);

  // ASSESS RESULTS
  bool fastAsLocalOnRepeatedAccess = true;

  for (unsigned cycle = 0; cycle < numDeviceHostCycles; ++cycle) {
    std::for_each(std::next(deviceSharedTimings[cycle].begin()),
                  deviceSharedTimings[cycle].end(), [&](const uint64_t timing) {
                    (timing < threshold * deviceLocalMean)
                        ? fastAsLocalOnRepeatedAccess &= true
                        : fastAsLocalOnRepeatedAccess &= false;
                  });

    std::for_each(std::next(hostSharedTimings[cycle].begin()),
                  hostSharedTimings[cycle].end(), [&](const uint64_t timing) {
                    (timing < threshold * hostLocalMean)
                        ? fastAsLocalOnRepeatedAccess &= true
                        : fastAsLocalOnRepeatedAccess &= false;
                  });
  }

  // CHECK IF PASSED
  bool passed = (fastAsLocalOnRepeatedAccess);

  // PRINT IF NOT PASSED
  if (!passed) {
    std::cout << "Page size as reported by os: " << getBytesPerPage()
              << " bytes \n";
    std::cout << "Allocating " << numPages
              << " pages of memory in SharedSpace.\n";

    std::cout << "Behavior found: \n";
    std::cout << "SharedSpace is as fast as local space on repeated access: "
              << fastAsLocalOnRepeatedAccess << ", we expect true \n\n";

    std::cout
        << "Please look at the following timings. The first access in a "
           "different ExecutionSpace is not evaluated for the test. As we "
           "expect the memory to migrate during the first access it might have "
           "a higher cycle count than subsequent accesses, depending on your "
           "hardware. If the cycles are more than "
        << threshold
        << " times the cycles for pure local memory access, we assume a page "
           "migration happened.\n\n";

    std::cout << "################SHARED SPACE####################\n";
    for (unsigned cycle = 0; cycle < numDeviceHostCycles; ++cycle) {
      std::cout << "DeviceExecutionSpace timings of run " << cycle << ":\n";
      printTimings(std::cout, deviceSharedTimings[cycle],
                   threshold * deviceLocalMean);
      std::cout << "HostExecutionSpace timings of run " << cycle << ":\n";
      printTimings(std::cout, hostSharedTimings[cycle],
                   threshold * hostLocalMean);
    }
    std::cout << "################LOCAL SPACE####################\n";
    printTimings(std::cout, deviceLocalTimings);
    printTimings(std::cout, hostLocalTimings);
  }
  ASSERT_TRUE(passed);
}
}  // namespace
