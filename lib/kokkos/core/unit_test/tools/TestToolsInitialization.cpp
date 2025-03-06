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

// This file calls most of the basic Kokkos primitives. When combined with a
// testing library this tests that our shared-library loading based profiling
// mechanisms work

#include <Kokkos_Core.hpp>
#include <gtest/gtest.h>

bool init_callback                               = false;
bool finalize_callback                           = false;
bool begin_parallel_for_callback                 = false;
bool end_parallel_for_callback                   = false;
bool begin_parallel_reduce_callback              = false;
bool end_parallel_reduce_callback                = false;
bool begin_parallel_scan_callback                = false;
bool end_parallel_scan_callback                  = false;
bool push_region_callback                        = false;
bool pop_region_callback                         = false;
bool allocate_data_callback                      = false;
bool deallocate_data_callback                    = false;
bool create_profile_section_callback             = false;
bool start_profile_section_callback              = false;
bool stop_profile_section_callback               = false;
bool destroy_profile_section_callback            = false;
bool profile_event_callback                      = false;
bool begin_deep_copy_callback                    = false;
bool end_deep_copy_callback                      = false;
bool begin_fence_callback                        = false;
bool end_fence_callback                          = false;
bool declare_metadata_callback                   = false;
bool request_tool_settings_callback              = false;
bool provide_tool_programming_interface_callback = false;

void test_tools_initialization_with_callbacks() {
  Kokkos::Tools::Experimental::set_init_callback(
      [](const int /*loadseq*/, const uint64_t /*version*/,
         const uint32_t /*num_infos*/,
         Kokkos::Profiling::KokkosPDeviceInfo* /*infos*/) {
        init_callback = true;
      });
  Kokkos::Tools::Experimental::set_finalize_callback(
      []() { finalize_callback = true; });
  Kokkos::Tools::Experimental::set_begin_parallel_for_callback(
      [](const char* /*n*/, const uint32_t /*d*/, uint64_t* /*k*/) {
        begin_parallel_for_callback = true;
      });
  Kokkos::Tools::Experimental::set_end_parallel_for_callback(
      [](const uint64_t /*k*/) { end_parallel_for_callback = true; });
  Kokkos::Tools::Experimental::set_begin_parallel_reduce_callback(
      [](const char* /*n*/, const uint32_t /*d*/, uint64_t* /*k*/) {
        begin_parallel_reduce_callback = true;
      });
  Kokkos::Tools::Experimental::set_end_parallel_reduce_callback(
      [](const uint64_t /*k*/) { end_parallel_reduce_callback = true; });
  Kokkos::Tools::Experimental::set_begin_parallel_scan_callback(
      [](const char* /*n*/, const uint32_t /*d*/, uint64_t* /*k*/) {
        begin_parallel_scan_callback = true;
      });
  Kokkos::Tools::Experimental::set_end_parallel_scan_callback(
      [](const uint64_t /*k*/) { end_parallel_scan_callback = true; });
  Kokkos::Tools::Experimental::set_push_region_callback(
      [](const char* /*name*/) { push_region_callback = true; });
  Kokkos::Tools::Experimental::set_pop_region_callback(
      []() { pop_region_callback = true; });
  Kokkos::Tools::Experimental::set_allocate_data_callback(
      [](Kokkos::Tools::SpaceHandle /*handle*/, const char* /*name*/,
         const void* /*ptr*/,
         const uint64_t /*size*/) { allocate_data_callback = true; });
  Kokkos::Tools::Experimental::set_deallocate_data_callback(
      [](Kokkos::Tools::SpaceHandle /*handle*/, const char* /*name*/,
         const void* /*ptr*/,
         const uint64_t /*size*/) { deallocate_data_callback = true; });
  Kokkos::Tools::Experimental::set_create_profile_section_callback(
      [](const char* /*name*/, uint32_t* /*id*/) {
        create_profile_section_callback = true;
      });
  Kokkos::Tools::Experimental::set_destroy_profile_section_callback(
      [](const uint32_t /*id*/) { destroy_profile_section_callback = true; });
  Kokkos::Tools::Experimental::set_start_profile_section_callback(
      [](uint32_t /*id*/) { start_profile_section_callback = true; });
  Kokkos::Tools::Experimental::set_stop_profile_section_callback(
      [](uint32_t /*id*/) { stop_profile_section_callback = true; });
  Kokkos::Tools::Experimental::set_profile_event_callback(
      [](const char* /*name*/) { profile_event_callback = true; });
  Kokkos::Tools::Experimental::set_begin_deep_copy_callback(
      [](Kokkos::Tools::SpaceHandle /*dst_handle*/, const char* /*dst_name*/,
         const void* /*dst_ptr*/, Kokkos::Tools::SpaceHandle /*src_handle*/,
         const char* /*src_name*/, const void* /*src_ptr*/,
         uint64_t /*size*/) { begin_deep_copy_callback = true; });
  Kokkos::Tools::Experimental::set_end_deep_copy_callback(
      []() { end_deep_copy_callback = true; });
  Kokkos::Tools::Experimental::set_begin_fence_callback(
      [](const char* /*n*/, const uint32_t /*d*/, uint64_t* /*k*/) {
        begin_fence_callback = true;
      });
  Kokkos::Tools::Experimental::set_end_fence_callback(
      [](const uint64_t /*k*/) { end_fence_callback = true; });
  Kokkos::Tools::Experimental::set_declare_metadata_callback(
      [](const char* /*key*/, const char* /*value*/) {
        declare_metadata_callback = true;
      });
  Kokkos::Tools::Experimental::set_request_tool_settings_callback(
      [](const uint32_t /*num_settings*/,
         Kokkos::Tools::Experimental::ToolSettings* /*settings*/) {
        request_tool_settings_callback = true;
      });
  Kokkos::Tools::Experimental::set_provide_tool_programming_interface_callback(
      [](const uint32_t /*num_functions*/,
         Kokkos::Tools::Experimental::ToolProgrammingInterface /*interface*/) {
        provide_tool_programming_interface_callback = true;
      });

  Kokkos::initialize();
  {
    ASSERT_TRUE(init_callback);
    ASSERT_FALSE(finalize_callback);
    ASSERT_FALSE(begin_parallel_for_callback);
    ASSERT_FALSE(end_parallel_for_callback);
    ASSERT_FALSE(begin_parallel_reduce_callback);
    ASSERT_FALSE(end_parallel_reduce_callback);
    ASSERT_FALSE(begin_parallel_scan_callback);
    ASSERT_FALSE(end_parallel_scan_callback);
    ASSERT_FALSE(push_region_callback);
    ASSERT_FALSE(pop_region_callback);
    ASSERT_FALSE(allocate_data_callback);
    ASSERT_FALSE(deallocate_data_callback);
    ASSERT_FALSE(create_profile_section_callback);
    ASSERT_FALSE(start_profile_section_callback);
    ASSERT_FALSE(stop_profile_section_callback);
    ASSERT_FALSE(destroy_profile_section_callback);
    ASSERT_FALSE(profile_event_callback);
    ASSERT_FALSE(begin_deep_copy_callback);
    ASSERT_FALSE(end_deep_copy_callback);
    ASSERT_FALSE(begin_fence_callback);
    ASSERT_FALSE(end_fence_callback);
    ASSERT_TRUE(declare_metadata_callback);
    ASSERT_TRUE(request_tool_settings_callback);
    ASSERT_TRUE(provide_tool_programming_interface_callback);
  }
  {
    using execution_space = Kokkos::DefaultExecutionSpace;
    using memory_space    = typename execution_space::memory_space;
    Kokkos::View<int*, memory_space> src_view("source", 10);
    Kokkos::View<int*, memory_space> dst_view("destination", 10);
    Kokkos::deep_copy(dst_view, src_view);
    Kokkos::parallel_for(
        "parallel_for", Kokkos::RangePolicy<execution_space>(0, 1),
        KOKKOS_LAMBDA(int i) { (void)i; });
    int result;
    Kokkos::parallel_reduce(
        "parallel_reduce", Kokkos::RangePolicy<execution_space>(0, 1),
        KOKKOS_LAMBDA(int i, int& hold_result) { hold_result += i; }, result);
    Kokkos::parallel_scan(
        "parallel_scan", Kokkos::RangePolicy<execution_space>(0, 1),
        KOKKOS_LAMBDA(const int i, int& hold_result, const bool final) {
          if (final) {
            hold_result += i;
          }
        });
    Kokkos::Profiling::pushRegion("push_region");
    Kokkos::Profiling::popRegion();
    uint32_t sectionId;
    Kokkos::Profiling::createProfileSection("created_section", &sectionId);
    Kokkos::Profiling::startSection(sectionId);
    Kokkos::Profiling::stopSection(sectionId);
    Kokkos::Profiling::destroyProfileSection(sectionId);
    Kokkos::Profiling::markEvent("profiling_event");
    Kokkos::Tools::declareMetadata("dogs", "good");
  }
  Kokkos::finalize();
  {
    ASSERT_TRUE(init_callback);
    ASSERT_TRUE(finalize_callback);
    ASSERT_TRUE(begin_parallel_for_callback);
    ASSERT_TRUE(end_parallel_for_callback);
    ASSERT_TRUE(begin_parallel_reduce_callback);
    ASSERT_TRUE(end_parallel_reduce_callback);
    ASSERT_TRUE(begin_parallel_scan_callback);
    ASSERT_TRUE(end_parallel_scan_callback);
    ASSERT_TRUE(push_region_callback);
    ASSERT_TRUE(pop_region_callback);
    ASSERT_TRUE(allocate_data_callback);
    ASSERT_TRUE(deallocate_data_callback);
    ASSERT_TRUE(create_profile_section_callback);
    ASSERT_TRUE(start_profile_section_callback);
    ASSERT_TRUE(stop_profile_section_callback);
    ASSERT_TRUE(destroy_profile_section_callback);
    ASSERT_TRUE(profile_event_callback);
    ASSERT_TRUE(begin_deep_copy_callback);
    ASSERT_TRUE(end_deep_copy_callback);
    ASSERT_TRUE(begin_fence_callback);
    ASSERT_TRUE(end_fence_callback);
    ASSERT_TRUE(declare_metadata_callback);
    ASSERT_TRUE(request_tool_settings_callback);
    ASSERT_TRUE(provide_tool_programming_interface_callback);
  }
}

TEST(tools, initialization_with_callbacks) {
  test_tools_initialization_with_callbacks();
}
