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

namespace Test {

void debug_print(const Kokkos_Profiling_SpaceHandle hand, const char* name,
                 const void* ptr, const size_t size) {
  std::cout << "Alloc: " << hand.name << ", [" << name << "," << ptr << "] "
            << size << std::endl;
}
void debug_dealloc(const Kokkos_Profiling_SpaceHandle hand, const char* name,
                   const void* ptr, const size_t size) {
  std::cout << "Dealloc: " << hand.name << ", [" << name << "," << ptr << "] "
            << size << std::endl;
}

void fail_on_event(const Kokkos::Profiling::SpaceHandle, const char*,
                   const void*, const uint64_t) {
  ASSERT_TRUE(false) << "Unexpected memory event";
}

void expect_no_events() {
  Kokkos::Tools::Experimental::set_allocate_data_callback(&fail_on_event);
  Kokkos::Tools::Experimental::set_deallocate_data_callback(&fail_on_event);
}

std::string expected_view_name;
std::string expected_space_name;
std::string error_message;
void expect_allocation_event(const std::string evn, const std::string esn,
                             const std::string em) {
  expected_view_name  = evn;
  expected_space_name = esn;
  error_message       = em;
  Kokkos::Tools::Experimental::set_allocate_data_callback(
      [](const Kokkos_Profiling_SpaceHandle hand, const char* name, const void*,
         const uint64_t) {
        ASSERT_EQ(std::string(hand.name), expected_space_name)
            << error_message << " (bad handle)";
        ASSERT_EQ(std::string(name), expected_view_name)
            << error_message << " (bad view name)";
        expect_no_events();
      });
}
void expect_deallocation_event(const std::string& evn, const std::string& esn,
                               const std::string em) {
  expected_view_name  = evn;
  expected_space_name = esn;
  error_message       = em;
  Kokkos::Tools::Experimental::set_deallocate_data_callback(
      [](const Kokkos_Profiling_SpaceHandle hand, const char* name, const void*,
         const uint64_t) {
        ASSERT_EQ(std::string(hand.name), expected_space_name)
            << error_message << " (bad handle)";
        ASSERT_EQ(std::string(name), expected_view_name)
            << error_message << " (bad view name)";
        expect_no_events();
      });
}

struct TestSpaceNamer {
  static constexpr const char* get_name() { return "TestSpace"; }
};
struct TestSpaceNamerTwo {
  static constexpr const char* get_name() { return "YoDawg"; }
};
struct TestSpaceNamerThree {
  static constexpr const char* get_name() { return "CustomAccessSpace"; }
};
using fake_memory_space = Kokkos::Experimental::LogicalMemorySpace<
    Kokkos::HostSpace, Kokkos::DefaultHostExecutionSpace, TestSpaceNamer,
    Kokkos::Experimental::LogicalSpaceSharesAccess::shared_access>;

void test_view_construct() {
  {
    expect_allocation_event("puppy_view", "TestSpace", "View allocation");
    Kokkos::View<double*, fake_memory_space> pup_view("puppy_view", 1000);
    expect_deallocation_event("puppy_view", "TestSpace", "View free");
  }
  Kokkos::Tools::Experimental::pause_tools();
}
void test_malloc_free() {
  expect_allocation_event("does_malloc_work", "TestSpace",
                          "Error in malloc event");
  auto* temp =
      Kokkos::kokkos_malloc<fake_memory_space>("does_malloc_work", 1000);
  expect_deallocation_event("does_malloc_work", "TestSpace", "Error in free");
  Kokkos::kokkos_free(temp);
  Kokkos::Tools::Experimental::pause_tools();
}
void test_chained_spaces() {
  using doubly_fake_memory_space = Kokkos::Experimental::LogicalMemorySpace<
      fake_memory_space, Kokkos::DefaultHostExecutionSpace, TestSpaceNamerTwo,
      Kokkos::Experimental::LogicalSpaceSharesAccess::shared_access>;
  {
    expect_allocation_event("xzibit_dot_jpeg", "YoDawg",
                            "Chained space view allocation");
    Kokkos::View<double*, doubly_fake_memory_space> pup_view("xzibit_dot_jpeg",
                                                             1000);
    expect_deallocation_event("xzibit_dot_jpeg", "YoDawg",
                              "Chained space free");
  }
  Kokkos::Tools::Experimental::pause_tools();
}
void test_space_allocations() {
  fake_memory_space debug_space;
  expect_allocation_event("allocation_from_space", "TestSpace",
                          "Space allocation");
  auto* temp = debug_space.allocate("allocation_from_space", 1000);
  expect_deallocation_event("allocation_from_space", "TestSpace",
                            "Space deallocation");
  debug_space.deallocate("allocation_from_space", temp, 1000);
  Kokkos::Tools::Experimental::pause_tools();
}
template <typename Space>
struct AccessCheckKernel {
  Kokkos::View<double*, Space> data;
  KOKKOS_FUNCTION void operator()(const int i) const { data[i] = i; }
};

template <typename Space>
void test_allowed_access() {
  constexpr const int data_size = 1000;
  Kokkos::View<double*, Space> test_view("test_view", data_size);
  AccessCheckKernel<Space> functor{test_view};
  Kokkos::parallel_for(
      "access_allowed",
      Kokkos::RangePolicy<typename Space::execution_space>(0, data_size),
      functor);
  Kokkos::fence();
}

using semantically_independent_logical_space =
    Kokkos::Experimental::LogicalMemorySpace<
        Kokkos::HostSpace, Kokkos::DefaultHostExecutionSpace,
        TestSpaceNamerThree,
        Kokkos::Experimental::LogicalSpaceSharesAccess::no_shared_access>;

TEST(defaultdevicetype, logical_space_views) { test_view_construct(); }
TEST(defaultdevicetype, logical_space_malloc) { test_malloc_free(); }
TEST(defaultdevicetype, logical_space_alloc) { test_space_allocations(); }
TEST(defaultdevicetype, chained_logical_spaces) { test_chained_spaces(); }
TEST(defaultdevicetype, access_allowed) {
  test_allowed_access<fake_memory_space>();
}
// FIXME_SYCL
#if !(defined(KOKKOS_COMPILER_INTEL) && defined(KOKKOS_ENABLE_SYCL))
TEST(defaultdevicetype_DeathTest, access_forbidden) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";
  ASSERT_DEATH(
      { test_allowed_access<semantically_independent_logical_space>(); },
      "Kokkos::View ERROR: attempt to access inaccessible memory space");
}
#endif

}  // namespace Test
