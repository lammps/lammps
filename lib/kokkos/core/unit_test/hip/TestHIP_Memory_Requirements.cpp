
/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

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
  KOKKOS_TEST_MEMORY_COARSEGRAINEDNESS(Kokkos::Experimental::HIPSpace, int, 10);
  KOKKOS_TEST_MEMORY_COARSEGRAINEDNESS(Kokkos::Experimental::HIPHostPinnedSpace,
                                       int, 10);
  KOKKOS_TEST_MEMORY_COARSEGRAINEDNESS(Kokkos::Experimental::HIPManagedSpace,
                                       int, 10);
}
}  // namespace
