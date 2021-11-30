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
#include <iostream>
#include <gtest/gtest.h>
#include "Kokkos_Core.hpp"

#include <impl/Kokkos_Stacktrace.hpp>
#include <vector>
#include <algorithm>
namespace Kokkos {
class Serial;
class OpenMP;
class Cuda;
class Threads;
namespace Experimental {
class SYCL;
class HIP;
class OpenMPTarget;
class HPX;
}  // namespace Experimental
}  // namespace Kokkos
namespace Test {
struct FencePayload {
  std::string name;
  enum distinguishable_devices { yes, no };
  distinguishable_devices distinguishable;
  uint32_t dev_id;
};

std::vector<FencePayload> found_payloads;
template <typename Lambda>
void expect_fence_events(std::vector<FencePayload>& expected, Lambda lam) {
  found_payloads = {};
  Kokkos::Tools::Experimental::set_begin_fence_callback(
      [](const char* name, const uint32_t dev_id, uint64_t*) {
        found_payloads.push_back(
            FencePayload{std::string(name),
                         FencePayload::distinguishable_devices::no, dev_id});
      });
  Kokkos::Tools::Experimental::set_begin_parallel_for_callback(
      [](const char* name, const uint32_t dev_id, uint64_t*) {
        found_payloads.push_back(
            FencePayload{std::string(name),
                         FencePayload::distinguishable_devices::no, dev_id});
      });
  lam();
  for (auto& entry : expected) {
    std::cout << "Ref: " << entry.dev_id << std::endl;
    std::cout << "Ref: " << entry.name << std::endl;
    auto search = std::find_if(
        found_payloads.begin(), found_payloads.end(),
        [&](const auto& found_entry) {
          auto name_match =
              (found_entry.name.find(entry.name) != std::string::npos);
          auto id_match = (entry.dev_id == found_entry.dev_id);
          std::cout << found_entry.dev_id << std::endl;
          std::cout << found_entry.name << std::endl;
          if (!name_match) {
            std::cout << "Miss on name\n";
          }
          if (!id_match) {
            std::cout << "Miss on id\n";
          }
          return (name_match && id_match);
        });
    auto found = (search != found_payloads.end());
    ASSERT_TRUE(found);
  }
  Kokkos::Tools::Experimental::set_begin_fence_callback(
      [](const char*, const uint32_t, uint64_t*) {});
  Kokkos::Tools::Experimental::set_begin_parallel_for_callback(
      [](const char*, const uint32_t, uint64_t*) {});
}

template <class>
struct increment {
  constexpr static const int size = 0;
};
int num_instances = 1;
struct TestFunctor {
  KOKKOS_FUNCTION void operator()(const int) const {}
};
template <typename Lambda>
void test_wrapper(const Lambda& lambda) {
  if (!std::is_same<Kokkos::DefaultExecutionSpace, Kokkos::Serial>::value) {
    lambda();
  }
}
/**
 * Test that fencing an instance with a name yields a fence
 * event of that name, and the correct device ID
 */
TEST(defaultdevicetype, test_named_instance_fence) {
  test_wrapper([&]() {
    auto root = Kokkos::Tools::Experimental::device_id_root<
        Kokkos::DefaultExecutionSpace>();
    std::vector<FencePayload> expected{

        {"named_instance", FencePayload::distinguishable_devices::no,
         root + num_instances}};
    expect_fence_events(expected, [=]() {
      Kokkos::DefaultExecutionSpace ex;
      ex.fence("named_instance");
    });
    num_instances += increment<Kokkos::DefaultExecutionSpace>::size;
  });
}
/**
 * Test that fencing an instance without a name yields a fence
 * event of a correct name, and the correct device ID
 */
TEST(defaultdevicetype, test_unnamed_instance_fence) {
  test_wrapper([&]() {
    auto root = Kokkos::Tools::Experimental::device_id_root<
        Kokkos::DefaultExecutionSpace>();
    std::vector<FencePayload> expected{

        {"Unnamed Instance Fence", FencePayload::distinguishable_devices::no,
         root + num_instances}};
    expect_fence_events(expected, [=]() {
      Kokkos::DefaultExecutionSpace ex;
      ex.fence();
    });
    num_instances += increment<Kokkos::DefaultExecutionSpace>::size;
  });
}

/**
 * Test that invoking a global fence with a name yields a fence
 * event of a correct name, and fences the root of the default device
 */
TEST(defaultdevicetype, test_named_global_fence) {
  test_wrapper([&]() {
    auto root = Kokkos::Tools::Experimental::device_id_root<
        Kokkos::DefaultExecutionSpace>();

    std::vector<FencePayload> expected{

        {"test global fence", FencePayload::distinguishable_devices::no, root}};
    expect_fence_events(expected,
                        [=]() { Kokkos::fence("test global fence"); });
  });
}

/**
 * Test that invoking a global fence with no name yields a fence
 * event of a correct name, and fences the root of the default device
 */
TEST(defaultdevicetype, test_unnamed_global_fence) {
  test_wrapper([&]() {
    auto root = Kokkos::Tools::Experimental::device_id_root<
        Kokkos::DefaultExecutionSpace>();

    std::vector<FencePayload> expected{

        {"Unnamed Global Fence", FencePayload::distinguishable_devices::no,
         root}};
    expect_fence_events(expected, [=]() { Kokkos::fence(); });
    num_instances += increment<Kokkos::DefaultExecutionSpace>::size;
  });
}
/**
 * Test that creating two default instances and fencing both yields
 * fence on the same device ID, as these should yield the same instance
 */
TEST(defaultdevicetype, test_multiple_default_instances) {
  test_wrapper([&]() {
    std::vector<FencePayload> expected{};
    expect_fence_events(expected, [=]() {
      Kokkos::DefaultExecutionSpace ex1;
      Kokkos::DefaultExecutionSpace ex2;
      ex1.fence("named_instance_fence_one");
      ex2.fence("named_instance_fence_two");
    });
    ASSERT_TRUE(found_payloads[0].dev_id == found_payloads[1].dev_id);
  });
}

/**
 * Test that fencing and kernels yield events on the correct device ID's
 */
TEST(defaultdevicetype, test_kernel_sequence) {
  test_wrapper([&]() {
    auto root = Kokkos::Tools::Experimental::device_id_root<
        Kokkos::DefaultExecutionSpace>();
    std::vector<FencePayload> expected{

        {"named_instance", FencePayload::distinguishable_devices::no,
         root + num_instances},
        {"test_kernel", FencePayload::distinguishable_devices::no,
         root + num_instances}

    };
    expect_fence_events(expected, [=]() {
      Kokkos::DefaultExecutionSpace ex;
      TestFunctor tf;
      ex.fence("named_instance");
      Kokkos::parallel_for(
          "test_kernel",
          Kokkos::RangePolicy<Kokkos::DefaultExecutionSpace>(ex, 0, 1), tf);
    });
    num_instances += increment<Kokkos::DefaultExecutionSpace>::size;
  });
}
#ifdef KOKKOS_ENABLE_CUDA
/**
 * CUDA ONLY: test that creating instances from streams leads to events
 * on different device ID's
 */
TEST(defaultdevicetype, test_streams) {
  test_wrapper([&]() {
    // auto root = Kokkos::Tools::Experimental::device_id_root<
    //    Kokkos::DefaultExecutionSpace>();
    std::vector<FencePayload> expected{};
    expect_fence_events(expected, [=]() {
      cudaStream_t s1, s2;
      cudaStreamCreate(&s1);
      cudaStreamCreate(&s2);
      Kokkos::Cuda default_space;
      Kokkos::Cuda space_s1(s1);
      Kokkos::Cuda space_s2(s2);
      default_space.fence();
      space_s1.fence();
      space_s2.fence();
    });
    num_instances += increment<Kokkos::DefaultExecutionSpace>::size;
    found_payloads.erase(
        std::remove_if(found_payloads.begin(), found_payloads.end(),
                       [&](const auto& entry) {
                         return (
                             entry.name.find("Fence on space initialization") !=
                             std::string::npos);
                       }),
        found_payloads.end());
    ASSERT_TRUE(found_payloads[0].dev_id != found_payloads[1].dev_id);
    ASSERT_TRUE(found_payloads[2].dev_id != found_payloads[1].dev_id);
    ASSERT_TRUE(found_payloads[2].dev_id != found_payloads[0].dev_id);
  });
}

#endif

}  // namespace Test
