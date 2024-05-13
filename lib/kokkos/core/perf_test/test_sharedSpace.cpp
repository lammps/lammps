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

#if defined _WIN32
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
void printTimings(std::ostream& out, std::vector<double> const& tr,
                  size_t numBytes,
                  double threshold = (std::numeric_limits<double>::max)()) {
  out << "TimingResult contains " << tr.size() << " results:\n";
  for (auto it = tr.begin(); it != tr.end(); ++it) {
    out << "Duration of loop " << it - tr.begin() << " is " << *it
        << " seconds.";
    if ((*it) > threshold) {
      out << " Marked as page migation.";
    }
    out << " The transfer rate is "
        << (double)numBytes / std::pow(1000.0, 3) / (*it) *
               2.0  // as we read and write
        << " GB/s \n";
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
  void operator()(const index_type idx) const { ++view_(idx); }
};

// TIMING CAPTURED KERNEL
// PREMISE: This kernel should always be memory bound, as we are measuring
// memory access times. The compute load of an increment is small enough on
// current hardware but this could be different for new hardware. As we count
// the wall-clock time in the kernel, the core frequency of the device has to be
// at the maximum to guarantee accurate masurements.
template <typename ExecSpace, typename ViewType>
std::vector<double> incrementInLoop(ViewType& view,
                                    unsigned int numRepetitions) {
  using index_type = decltype(view.size());
  Kokkos::Timer timer;
  std::vector<double> results;

  Kokkos::fence();
  for (unsigned i = 0; i < numRepetitions; ++i) {
    IncrementFunctor<ViewType> func(view);
    timer.reset();
    Kokkos::parallel_for(
        "increment",
        Kokkos::RangePolicy<ExecSpace, Kokkos::IndexType<index_type>>{
            0, view.size()},
        func);
    Kokkos::fence();
    results.push_back(timer.seconds());
  }
  return results;
}

size_t getDeviceMemorySize() {
#if defined KOKKOS_ENABLE_CUDA
  return Kokkos::Cuda{}.cuda_device_prop().totalGlobalMem;
#elif defined KOKKOS_ENABLE_HIP
  return Kokkos::HIP{}.hip_device_prop().totalGlobalMem;
#elif defined KOKKOS_ENABLE_SYCL
  auto device = Kokkos::Experimental::SYCL{}.sycl_queue().get_device();
  return device.get_info<sycl::info::device::global_mem_size>();
#else
#error \
    "The sharedMemory test is only defined for Kokkos::Cuda, Kokkos::HIP, and Kokkos::SYCL"
  return 0;
#endif
}

struct Arguments {
  unsigned int numRepetitions       = 10;
  unsigned int numWarmupRepetitions = 100;
  unsigned int numDeviceHostCycles  = 3;
  double fractionOfDeviceMemory     = 0.4;
  double threshold                  = 2.0;
};

void test_sharedSpace(Arguments args) {
  const unsigned int numRepetitions       = args.numRepetitions;
  const unsigned int numWarmupRepetitions = args.numWarmupRepetitions;
  const unsigned int numDeviceHostCycles  = args.numDeviceHostCycles;
  double fractionOfDeviceMemory           = args.fractionOfDeviceMemory;
  double threshold                        = args.threshold;
  size_t numBytes = fractionOfDeviceMemory * getDeviceMemorySize();
  size_t numPages = numBytes / getBytesPerPage();

  // ALLOCATION
  Kokkos::View<int*, Kokkos::SharedSpace> migratableData(
      "migratableData", numPages * getBytesPerPage() / sizeof(int));
  Kokkos::View<int*, Kokkos::DefaultExecutionSpace::memory_space> deviceData(
      "deviceData", numPages * getBytesPerPage() / sizeof(int));
  Kokkos::View<int*, Kokkos::DefaultHostExecutionSpace::memory_space> hostData(
      "hostData", numPages * getBytesPerPage() / sizeof(int));
  Kokkos::fence();

  // WARMUP GPU
  incrementInLoop<Kokkos::DefaultExecutionSpace>(
      deviceData,
      numWarmupRepetitions);  // warming up gpu

  // GET DEVICE LOCAL TIMINGS
  auto deviceLocalResults = incrementInLoop<Kokkos::DefaultExecutionSpace>(
      deviceData, numRepetitions);

  // WARMUP HOST
  incrementInLoop<Kokkos::DefaultHostExecutionSpace>(
      hostData,
      numWarmupRepetitions);  // warming up host
  // GET HOST LOCAL TIMINGS
  auto hostLocalResults = incrementInLoop<Kokkos::DefaultHostExecutionSpace>(
      hostData, numRepetitions);

  // GET PAGE MIGRATING TIMINGS DATA
  std::vector<decltype(deviceLocalResults)> deviceResults{};
  std::vector<decltype(hostLocalResults)> hostResults{};
  for (unsigned i = 0; i < numDeviceHostCycles; ++i) {
    // WARMUP GPU
    incrementInLoop<Kokkos::DefaultExecutionSpace>(
        deviceData,
        numWarmupRepetitions);  // warming up gpu without touching the
                                // migratableData to get measurements of initial
                                // position
    // GET RESULTS DEVICE
    deviceResults.push_back(incrementInLoop<Kokkos::DefaultExecutionSpace>(
        migratableData, numRepetitions));

    // WARMUP HOST
    incrementInLoop<Kokkos::DefaultHostExecutionSpace>(
        hostData,
        numWarmupRepetitions);  // warming up host without touching the
                                // migratableData to get measurements of initial
                                // position
    // GET RESULTS HOST
    hostResults.push_back(incrementInLoop<Kokkos::DefaultHostExecutionSpace>(
        migratableData, numRepetitions));
  }

  // COMPUTE STATISTICS OF HOST AND DEVICE LOCAL KERNELS
  auto hostLocalMean   = computeMean(hostLocalResults);
  auto deviceLocalMean = computeMean(deviceLocalResults);

  // ASSESS PAGE MIGRATIONS
  bool initialPlacementOnDevice   = false;
  bool migratesOnEverySpaceAccess = true;
  bool migratesOnlyOncePerAccess  = true;

  for (unsigned cycle = 0; cycle < numDeviceHostCycles; ++cycle) {
    unsigned int indicatedPageMigrationsDevice = std::count_if(
        deviceResults[cycle].begin(), deviceResults[cycle].end(),
        [&](auto const& val) { return val > (threshold * deviceLocalMean); });

    if (cycle == 0 && indicatedPageMigrationsDevice == 0)
      initialPlacementOnDevice = true;
    else {
      if (indicatedPageMigrationsDevice != 1) migratesOnlyOncePerAccess = false;
    }

    unsigned int indicatedPageMigrationsHost = std::count_if(
        hostResults[cycle].begin(), hostResults[cycle].end(),
        [&](auto const& val) { return val > (threshold * hostLocalMean); });

    if (indicatedPageMigrationsHost != 1) migratesOnlyOncePerAccess = false;

    if (cycle != 0 && indicatedPageMigrationsDevice != 1 &&
        indicatedPageMigrationsHost != 1)
      migratesOnEverySpaceAccess = false;
  }

  std::cout << "Page size as reported by os: " << getBytesPerPage()
            << " bytes \n";
  std::cout << "Allocating " << numPages
            << " pages of memory in pageMigratingMemorySpace.\n"
            << "This corresponds to " << fractionOfDeviceMemory * 100
            << " % of the device memory.\n"
            << "The view size is " << migratableData.size() << "\n";

  std::cout << "Behavior found: \n";
  std::cout << "Initial placement on device is " << initialPlacementOnDevice
            << "\n";
  std::cout << "Memory migrates on every space access is "
            << migratesOnEverySpaceAccess << "\n";
  std::cout << "Memory migrates only once per access "
            << migratesOnlyOncePerAccess << "\n\n";

  std::cout << "Please look at the following timings. A migration was "
               "marked detected if the time was larger than "
            << threshold * hostLocalMean << " for the host and "
            << threshold * deviceLocalMean << " for the device\n\n";

  std::cout << "#############TIMINGS WITH SHAREDSPACE##################\n";

  for (unsigned cycle = 0; cycle < numDeviceHostCycles; ++cycle) {
    std::cout << "device timings of run " << cycle << ":\n";
    printTimings(std::cout, deviceResults[cycle], numBytes,
                 threshold * deviceLocalMean);
    std::cout << "host timings of run " << cycle << ":\n";
    printTimings(std::cout, hostResults[cycle], numBytes,
                 threshold * hostLocalMean);
  }
  std::cout << "\n#############TIMINGS WITH LOCALSPACE##################\n";
  std::cout << "Device local memory timings for comparison:\n";
  printTimings(std::cout, deviceLocalResults, numBytes);
  std::cout << "Host local memory timings for comparison:\n";
  printTimings(std::cout, hostLocalResults, numBytes);
}
}  // namespace

int main(int argc, char* argv[]) {
  static const char help_flag[]                   = "--help";
  static const char numRepetitions_flag[]         = "--numRepetitions=";
  static const char numWarmupRepetitions_flag[]   = "--numWarmupRepetitions=";
  static const char numDeviceHostCycles_flag[]    = "--numDeviceHostCycles=";
  static const char fractionOfDeviceMemory_flag[] = "--fractionOfDeviceMemory=";
  static const char threshold_flag[]              = "--threshold=";

  int ask_help = 0;
  Arguments args;

  for (int i = 1; i < argc; i++) {
    const char* const a = argv[i];

    if (!strncmp(a, help_flag, strlen(help_flag))) ask_help = 1;

    if (!strncmp(a, numRepetitions_flag, strlen(numRepetitions_flag)))
      args.numRepetitions = std::stoi(a + strlen(numRepetitions_flag));

    if (!strncmp(a, numWarmupRepetitions_flag,
                 strlen(numWarmupRepetitions_flag)))
      args.numWarmupRepetitions =
          std::stoi(a + strlen(numWarmupRepetitions_flag));

    if (!strncmp(a, numDeviceHostCycles_flag, strlen(numDeviceHostCycles_flag)))
      args.numDeviceHostCycles =
          std::stoi(a + strlen(numDeviceHostCycles_flag));

    if (!strncmp(a, fractionOfDeviceMemory_flag,
                 strlen(fractionOfDeviceMemory_flag)))
      args.fractionOfDeviceMemory =
          std::stod(a + strlen(fractionOfDeviceMemory_flag));

    if (!strncmp(a, threshold_flag, strlen(threshold_flag)))
      args.threshold = std::stod(a + strlen(threshold_flag));
  }

  if (ask_help) {
    std::cout << "command line options:"
              << " " << help_flag << " " << numRepetitions_flag << "##"
              << " " << numWarmupRepetitions_flag << "##"
              << " " << numDeviceHostCycles_flag << "##"
              << " " << fractionOfDeviceMemory_flag << "##"
              << " " << threshold_flag << "##"
              << " any given Kokkos args are passed to Kokkos::initialize ##"
              << std::endl;
    return 0;
  }

  Kokkos::initialize(argc, argv);
  if constexpr (Kokkos::has_shared_space)
    test_sharedSpace(args);
  else
    std::cout
        << "The used Kokkos configuration does not support SharedSpace \n";
  Kokkos::finalize();

  return 0;
}
