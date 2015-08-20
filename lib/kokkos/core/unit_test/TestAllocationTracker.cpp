/*
//@HEADER
// ************************************************************************
// 
//                        Kokkos v. 2.0
//              Copyright (2014) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
// 
// ************************************************************************
//@HEADER
*/

#include <gtest/gtest.h>

#include <iostream>
#include <vector>

#include <Kokkos_Core.hpp>

#include <impl/Kokkos_AllocationTracker.hpp>
#include <impl/Kokkos_BasicAllocators.hpp>

namespace Test {

class alocation_tracker : public ::testing::Test {
protected:
  static void SetUpTestCase()
  {
    Kokkos::initialize();
  }

  static void TearDownTestCase()
  {
    Kokkos::finalize();
  }
};

TEST_F( alocation_tracker, simple)
{
  using namespace Kokkos::Impl;

  {
    AllocationTracker tracker;
    EXPECT_FALSE( tracker.is_valid() );
  }

  // test ref count and label
  {
    int size = 100;
    std::vector<AllocationTracker> trackers(size);

    trackers[0] = AllocationTracker( MallocAllocator(), 128,"Test");

    for (int i=0; i<size; ++i) {
      trackers[i] = trackers[0];
    }

    EXPECT_EQ(100u, trackers[0].ref_count());
    EXPECT_EQ(std::string("Test"), std::string(trackers[0].label()));
  }


  // test circular list
  {
    int num_allocs = 3000;
    unsigned ref_count = 100;

    std::vector<AllocationTracker> trackers(num_allocs);

    for (int i=0; i<num_allocs; ++i) {
      trackers[i] = AllocationTracker( MallocAllocator(), 128, "Test");
      std::vector<AllocationTracker> ref_trackers(ref_count);
      for (unsigned j=0; j<ref_count; ++j) {
        ref_trackers[j] = trackers[i];
      }
      EXPECT_EQ( ref_count + 1u, trackers[i].ref_count() );
    }

    for (int i=0; i<num_allocs; ++i) {
      EXPECT_EQ( 1u, trackers[i].ref_count() );
    }
  }
}

TEST_F( alocation_tracker, force_leaks)
{
// uncomment to force memory leaks
#if 0
  using namespace Kokkos::Impl;
  Kokkos::kokkos_malloc("Forced Leak", 4096*10);
  Kokkos::kokkos_malloc<Kokkos::HostSpace>("Forced Leak", 4096*10);
#endif
}

TEST_F( alocation_tracker, disable_reference_counting)
{
  using namespace Kokkos::Impl;
  // test ref count and label
  {
    int size = 100;
    std::vector<AllocationTracker> trackers(size);

    trackers[0] = AllocationTracker( MallocAllocator(), 128,"Test");

    for (int i=1; i<size; ++i) {
      trackers[i] = CopyWithoutTracking::apply(trackers[0]);
    }

    EXPECT_EQ(1u, trackers[0].ref_count());
    EXPECT_EQ(std::string("Test"), std::string(trackers[0].label()));
  }
}

} // namespace Test
