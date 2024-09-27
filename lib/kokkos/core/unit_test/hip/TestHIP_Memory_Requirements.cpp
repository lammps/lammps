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
#include <TestHIP_Category.hpp>

namespace {

template <class HIPMemoryContainer>
bool checkMemoryCoarseGrainedness(HIPMemoryContainer const& container) {
  auto size           = container.size();
  auto allocationSize = HIPMemoryContainer::required_allocation_size(size);
  hipMemRangeCoherencyMode memInfo;

  KOKKOS_IMPL_HIP_SAFE_CALL(hipMemRangeGetAttribute(
      &memInfo, sizeof(hipMemRangeCoherencyMode),
      hipMemRangeAttributeCoherencyMode, container.data(), allocationSize));

  return (hipMemRangeCoherencyModeCoarseGrain == memInfo);
}

#define KOKKOS_TEST_MEMORY_COARSEGRAINEDNESS(MEMORY_SPACE, DATATYPE, SIZE)    \
  {                                                                           \
    Kokkos::View<DATATYPE*, MEMORY_SPACE> view(#MEMORY_SPACE, SIZE);          \
    ASSERT_TRUE(view.is_allocated())                                          \
        << "View in " << #MEMORY_SPACE << " with size " << SIZE               \
        << " was not allocated. This prevents checks of the grainedness.";    \
    ASSERT_TRUE(checkMemoryCoarseGrainedness(view))                           \
        << "The memory in views in " << #MEMORY_SPACE                         \
        << " is not coarse-grained. Kokkos relies on all user facing memory " \
           "being coarse-grained.";                                           \
  }

TEST(hip, memory_requirements) {
  // we want all user-facing memory in hip to be coarse grained. As of
  // today(07.01.22) the documentation is not reliable/correct, we test the
  // memory on the device and host
  KOKKOS_TEST_MEMORY_COARSEGRAINEDNESS(Kokkos::HIPSpace, int, 10);
  KOKKOS_TEST_MEMORY_COARSEGRAINEDNESS(Kokkos::HIPHostPinnedSpace, int, 10);
  KOKKOS_TEST_MEMORY_COARSEGRAINEDNESS(Kokkos::HIPManagedSpace, int, 10);
}
}  // namespace
