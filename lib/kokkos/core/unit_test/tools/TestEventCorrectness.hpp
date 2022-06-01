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
#include "Kokkos_Core_fwd.hpp"
#include "include/ToolTestingUtilities.hpp"
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
using index_type  = Kokkos::RangePolicy<>::index_type;
struct TestFunctor {
  KOKKOS_FUNCTION void operator()(const index_type) const {}
};
struct TestReduceFunctor {
  using value_type = int;
  KOKKOS_FUNCTION void operator()(const index_type, value_type&) const {}
};
struct TestScanFunctor {
  using value_type = int;
  KOKKOS_FUNCTION void operator()(const index_type, value_type&, bool) const {}
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
TEST(kokkosp, test_named_instance_fence) {
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
TEST(kokkosp, test_unnamed_instance_fence) {
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
TEST(kokkosp, test_named_global_fence) {
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
TEST(kokkosp, test_unnamed_global_fence) {
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
TEST(kokkosp, test_multiple_default_instances) {
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
 * Test that device_id() and identifier_from_devid(id) are reciprocal
 * operations
 */
TEST(kokkosp, test_id_gen) {
  using namespace Kokkos::Tools::Experimental;
  using Kokkos::Tools::Experimental::DeviceTypeTraits;
  test_wrapper([&]() {
    Kokkos::DefaultExecutionSpace ex;
    auto id      = device_id(ex);
    auto id_ref  = identifier_from_devid(id);
    auto success = (id_ref.instance_id == ex.impl_instance_id()) &&
                   (id_ref.device_id ==
                    static_cast<uint32_t>(
                        DeviceTypeTraits<Kokkos::DefaultExecutionSpace>::id));
    ASSERT_TRUE(success);
  });
}

/**
 * Test that fencing and kernels yield events on the correct device ID's
 */
TEST(kokkosp, test_kernel_sequence) {
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
TEST(kokkosp, test_streams) {
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
                         return (entry.name.find("Unnamed Instance Fence") ==
                                 std::string::npos);
                       }),
        found_payloads.end());
    ASSERT_NE(found_payloads[0].dev_id, found_payloads[1].dev_id);
    ASSERT_NE(found_payloads[2].dev_id, found_payloads[1].dev_id);
    ASSERT_NE(found_payloads[2].dev_id, found_payloads[0].dev_id);
  });
}

#endif
/** FIXME: OpenMPTarget currently has unexpected fences */
#ifndef KOKKOS_ENABLE_OPENMPTARGET
TEST(kokkosp, async_deep_copy) {
  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableFences());
  Kokkos::View<float*> left("left", 5), right("right", 5);

  auto success = validate_absence(
      [&]() {
        Kokkos::deep_copy(Kokkos::DefaultExecutionSpace(), left, right);
      },
      [&](BeginFenceEvent begin) {
        if (begin.deviceID !=
            Kokkos::DefaultExecutionSpace().impl_instance_id()) {
          std::stringstream error_message;
          error_message
              << "Fence encountered outside of the default instance, default: "
              << Kokkos::DefaultExecutionSpace().impl_instance_id()
              << ", encountered " << begin.deviceID << " , fence name "
              << begin.name;
          return MatchDiagnostic{true, {error_message.str()}};
        }
        return MatchDiagnostic{false};
      });
  ASSERT_TRUE(success);
}
#endif
TEST(kokkosp, parallel_for) {
  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableKernels());
  auto success = validate_event_set(
      [=]() {
        TestFunctor tf;
        Kokkos::parallel_for("dogs", Kokkos::RangePolicy<>(0, 1), tf);
      },
      [=](BeginParallelForEvent begin_event, EndParallelForEvent end_event) {
        if (begin_event.name != "dogs") {
          return MatchDiagnostic{false, {"No match on BeginParallelFor name"}};
        }
        if (end_event.kID != ((begin_event.kID))) {
          return MatchDiagnostic{false, {"No match on kID's"}};
        }
        return MatchDiagnostic{true};
      });
  ASSERT_TRUE(success);
}

TEST(kokkosp, parallel_reduce) {
  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableKernels());
  auto success = validate_event_set(
      [=]() {
        TestReduceFunctor tf;
        int result;
        Kokkos::parallel_reduce("dogs", Kokkos::RangePolicy<>(0, 1), tf,
                                Kokkos::Sum<int>(result));
      },
      [=](BeginParallelReduceEvent begin_event,
          EndParallelReduceEvent end_event) {
        if (begin_event.name != "dogs") {
          return MatchDiagnostic{false,
                                 {"No match on BeginParallelReduce name"}};
        }
        if (end_event.kID != ((begin_event.kID))) {
          return MatchDiagnostic{false, {"No match on kID's"}};
        }
        return MatchDiagnostic{true};
      });
  ASSERT_TRUE(success);
}

TEST(kokkosp, parallel_scan) {
  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableKernels());
  auto success = validate_event_set(
      [=]() {
        TestScanFunctor tf;
        Kokkos::parallel_scan("dogs", Kokkos::RangePolicy<>(0, 1), tf);
      },
      [=](BeginParallelScanEvent begin_event, EndParallelScanEvent end_event) {
        if (begin_event.name != "dogs") {
          return MatchDiagnostic{false, {"No match on BeginParallelScan name"}};
        }
        if (end_event.kID != ((begin_event.kID))) {
          return MatchDiagnostic{false, {"No match on kID's"}};
        }
        return MatchDiagnostic{true};
      });
// Currently, this test is known to fail with OpenMPTarget
#ifndef KOKKOS_ENABLE_OPENMPTARGET
  ASSERT_TRUE(success);
#else
  (void)success;
#endif
}

TEST(kokkosp, regions) {
  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableRegions());
  auto success = validate_event_set(
      [=]() {
        Kokkos::Tools::pushRegion("dogs");
        Kokkos::Tools::popRegion();
      },
      [=](PushRegionEvent push_event, PopRegionEvent) {
        if (push_event.name != "dogs") {
          return MatchDiagnostic{false, {"No match on PushRegion name"}};
        }

        return MatchDiagnostic{true};
      });
  ASSERT_TRUE(success);
}

TEST(kokkosp, fences) {
  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableFences());
  auto success = validate_event_set(
      [=]() { Kokkos::DefaultExecutionSpace().fence("dogs"); },
      [=](BeginFenceEvent begin_event, EndFenceEvent end_event) {
        if (begin_event.name != "dogs") {
          return MatchDiagnostic{false, {"No match on BeginFence name"}};
        }
        if (end_event.kID != ((begin_event.kID))) {
          return MatchDiagnostic{false, {"No match on kID's"}};
        }
        return MatchDiagnostic{true};
      });
  ASSERT_TRUE(success);
}

TEST(kokkosp, raw_allocation) {
  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableAllocs());
  auto success = validate_event_set(
      [=]() {
        void* foo =
            Kokkos::kokkos_malloc<Kokkos::DefaultExecutionSpace::memory_space>(
                "dogs", 1000);
        Kokkos::kokkos_free(foo);
      },
      [=](AllocateDataEvent alloc, DeallocateDataEvent free) {
        if (alloc.name != "dogs") {
          return MatchDiagnostic{false, {"No match on alloc name"}};
        }
        if (alloc.size != 1000) {
          return MatchDiagnostic{false, {"No match on alloc size"}};
        }
        if (alloc.ptr != free.ptr) {
          return MatchDiagnostic{false, {"No match on pointers"}};
        }
        if (free.name != "dogs") {
          return MatchDiagnostic{false, {"No match on free name"}};
        }
        if (free.size != 1000) {
          return MatchDiagnostic{false, {"No match on free size"}};
        }
        return MatchDiagnostic{true};
      });
// Currently, this test is known to fail with OpenMPTarget
#ifndef KOKKOS_ENABLE_OPENMPTARGET
  ASSERT_TRUE(success);
#else
  (void)success;
#endif
}

TEST(kokkosp, view) {
  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableAllocs());
  auto success = validate_event_set(
      [=]() { Kokkos::View<float*> dogs("dogs", 1000); },
      [=](AllocateDataEvent alloc, DeallocateDataEvent free) {
        if (alloc.name != "dogs") {
          return MatchDiagnostic{false, {"No match on alloc name"}};
        }
        if (alloc.size != 1000 * sizeof(float)) {
          return MatchDiagnostic{false, {"No match on alloc size"}};
        }
        if (alloc.ptr != free.ptr) {
          return MatchDiagnostic{false, {"No match on pointers"}};
        }
        if (free.name != "dogs") {
          return MatchDiagnostic{false, {"No match on free name"}};
        }
        if (free.size != 1000 * sizeof(float)) {
          return MatchDiagnostic{false, {"No match on free size"}};
        }
        return MatchDiagnostic{true};
      });
// Currently, this test is known to fail with OpenMPTarget
#ifndef KOKKOS_ENABLE_OPENMPTARGET
  ASSERT_TRUE(success);
#else
  (void)success;
#endif
}

TEST(kokkosp, sections) {
  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableSections());
  auto success = validate_event_set(
      [=]() {
        uint32_t section_id;
        Kokkos::Tools::createProfileSection("dogs", &section_id);
        Kokkos::Tools::startSection(section_id);
        Kokkos::Tools::stopSection(section_id);
        Kokkos::Tools::destroyProfileSection(section_id);
      },
      [=](CreateProfileSectionEvent create, StartProfileSectionEvent start,
          StopProfileSectionEvent stop, DestroyProfileSectionEvent destroy) {
        if (create.name != "dogs") {
          return MatchDiagnostic{false, {"No match on section name"}};
        }
        if ((create.id != start.id) || (stop.id != start.id) ||
            (stop.id != destroy.id)) {
          return MatchDiagnostic{false, {"No match on section IDs"}};
        }
        return MatchDiagnostic{true};
      });
  ASSERT_TRUE(success);
}

TEST(kokkosp, metadata) {
  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableMetadata());
  auto success = validate_event_set(
      [=]() {
        /** Attempts to decrease the value of dog_goodness will be rejected on
         * review */
        Kokkos::Tools::declareMetadata("dog_goodness", "infinity");
      },
      [=](DeclareMetadataEvent meta) {
        if (meta.key != "dog_goodness") {
          return MatchDiagnostic{false, {"No match on metadata key"}};
        }
        if (meta.value != "infinity") {
          return MatchDiagnostic{false, {"No match on metadata value"}};
        }
        return MatchDiagnostic{true};
      });
  ASSERT_TRUE(success);
}

TEST(kokkosp, profile_events) {
  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableProfileEvents());
  auto success = validate_event_set(
      [=]() { Kokkos::Tools::markEvent("dog_goodness>=infinity"); },
      [=](ProfileEvent event) {
        if (event.name != "dog_goodness>=infinity") {
          return MatchDiagnostic{false, {"No match on profiled event name"}};
        }
        return MatchDiagnostic{true};
      });
  ASSERT_TRUE(success);
}
#if defined(KOKKOS_ENABLE_TUNING)
TEST(kokkosp, tuning_sequence) {
  using namespace Kokkos::Test::Tools;
  listen_tool_events(Config::DisableAll(), Config::EnableTuning());
  size_t input_id, output_id;
  Kokkos::Tools::Experimental::VariableInfo input_info;
  input_info.type = Kokkos::Tools::Experimental::ValueType::kokkos_value_int64;
  input_info.category = Kokkos::Tools::Experimental::StatisticalCategory::
      kokkos_value_categorical;
  input_info.valueQuantity =
      Kokkos::Tools::Experimental::CandidateValueType::kokkos_value_unbounded;
  Kokkos::Tools::Experimental::VariableInfo output_info = input_info;
  output_info.valueQuantity =
      Kokkos::Tools::Experimental::CandidateValueType::kokkos_value_set;
  std::vector<int64_t> values{1, 2, 3, 4, 5};
  output_info.candidates = Kokkos::Tools::Experimental::make_candidate_set(
      values.size(), values.data());

  auto success = validate_event_set(
      [&]() {
        input_id = Kokkos::Tools::Experimental::declare_input_type("input.dogs",
                                                                   input_info);
        output_id = Kokkos::Tools::Experimental::declare_output_type(
            "output.dogs", output_info);
        auto next_context = Kokkos::Tools::Experimental::get_new_context_id();
        Kokkos::Tools::Experimental::begin_context(next_context);
        Kokkos::Tools::Experimental::VariableValue feature_value =
            Kokkos::Tools::Experimental::make_variable_value(input_id,
                                                             int64_t(0));
        Kokkos::Tools::Experimental::VariableValue tuning_value =
            Kokkos::Tools::Experimental::make_variable_value(output_id,
                                                             int64_t(1));
        Kokkos::Tools::Experimental::set_input_values(next_context, 1,
                                                      &feature_value);
        Kokkos::Tools::Experimental::request_output_values(next_context, 1,
                                                           &tuning_value);
        Kokkos::Tools::Experimental::end_context(next_context);
      },
      [&](DeclareInputTypeEvent input, DeclareOutputTypeEvent output) {
        if (input.variable_id != input_id) {
          return MatchDiagnostic{false, {"No match on input id"}};
        }
        if (output.variable_id != output_id) {
          return MatchDiagnostic{false, {"No match on output id"}};
        }
        if (output.info.candidates.set.size != 5) {
          return MatchDiagnostic{
              false, {"Candidates not properly passed through tuning system"}};
        }
        return MatchDiagnostic{true};
      },
      [=](BeginContextEvent) { return MatchDiagnostic{true}; },
      [&](RequestOutputValuesEvent value_request) {
        if (value_request.inputs[0].metadata->type != input_info.type) {
          return MatchDiagnostic{false, {"No match on input in request"}};
        }
        if (value_request.outputs[0].metadata->type != output_info.type) {
          return MatchDiagnostic{false, {"No match on output in request"}};
        }
        return MatchDiagnostic{true};
      },
      [=](EndContextEvent) { return MatchDiagnostic{true}; });
  ASSERT_TRUE(success);
}
#endif
TEST(kokkosp, no_init_kernel) {
  using namespace Kokkos::Test::Tools;

  listen_tool_events(Config::DisableAll(), Config::EnableKernels());
  auto success = validate_absence(
      [=]() {
        Kokkos::View<float*> not_inited(
            Kokkos::ViewAllocateWithoutInitializing("no_inits_here_dog"), 100);
      },
      [=](BeginParallelForEvent) {
        return MatchDiagnostic{true, {"Found begin event"}};
      },
      [=](EndParallelForEvent) {
        return MatchDiagnostic{true, {"Found end event"}};
      });
  ASSERT_TRUE(success);
}

TEST(kokkosp, get_events) {
  using namespace Kokkos::Test::Tools;
  auto event_vector = get_event_set([=]() {
    Kokkos::Tools::pushRegion("dogs");
    Kokkos::Tools::popRegion();
  });
  for (const auto& ptr : event_vector) {
    auto ptr_as_begin = std::dynamic_pointer_cast<BeginParallelForEvent>(ptr);
    ASSERT_TRUE(ptr_as_begin == nullptr);
  }
}
}  // namespace Test
