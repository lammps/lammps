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

#ifndef KOKKOS_IMPL_KOKKOS_PROFILING_HPP
#define KOKKOS_IMPL_KOKKOS_PROFILING_HPP

#include <impl/Kokkos_Profiling_Interface.hpp>
#include <memory>
#include <iosfwd>
#include <unordered_map>
#include <map>
#include <string>
#include <type_traits>
#include <mutex>
namespace Kokkos {

// forward declaration
bool tune_internals() noexcept;

namespace Tools {

struct InitArguments {
  // NOTE DZP: PossiblyUnsetOption was introduced
  // before C++17, std::optional is a better choice
  // for this long-term
  static const std::string unset_string_option;
  enum PossiblyUnsetOption { unset, off, on };
  PossiblyUnsetOption tune_internals = unset;
  PossiblyUnsetOption help           = unset;
  std::string lib                    = unset_string_option;
  std::string args                   = unset_string_option;
};

namespace Impl {

struct InitializationStatus {
  enum InitializationResult {
    success,
    failure,
    help_request,
    environment_argument_mismatch
  };
  InitializationResult result;
  std::string error_message;
};
InitializationStatus initialize_tools_subsystem(
    const Kokkos::Tools::InitArguments& args);

void parse_command_line_arguments(int& narg, char* arg[],
                                  InitArguments& arguments);
Kokkos::Tools::Impl::InitializationStatus parse_environment_variables(
    InitArguments& arguments);

}  // namespace Impl

bool profileLibraryLoaded();

void beginParallelFor(const std::string& kernelPrefix, const uint32_t devID,
                      uint64_t* kernelID);
void endParallelFor(const uint64_t kernelID);
void beginParallelScan(const std::string& kernelPrefix, const uint32_t devID,
                       uint64_t* kernelID);
void endParallelScan(const uint64_t kernelID);
void beginParallelReduce(const std::string& kernelPrefix, const uint32_t devID,
                         uint64_t* kernelID);
void endParallelReduce(const uint64_t kernelID);

void pushRegion(const std::string& kName);
void popRegion();

void createProfileSection(const std::string& sectionName, uint32_t* secID);
void startSection(const uint32_t secID);
void stopSection(const uint32_t secID);
void destroyProfileSection(const uint32_t secID);

void markEvent(const std::string& evName);

void allocateData(const SpaceHandle space, const std::string label,
                  const void* ptr, const uint64_t size);
void deallocateData(const SpaceHandle space, const std::string label,
                    const void* ptr, const uint64_t size);

void beginDeepCopy(const SpaceHandle dst_space, const std::string dst_label,
                   const void* dst_ptr, const SpaceHandle src_space,
                   const std::string src_label, const void* src_ptr,
                   const uint64_t size);
void endDeepCopy();
void beginFence(const std::string name, const uint32_t deviceId,
                uint64_t* handle);
void endFence(const uint64_t handle);

/**
 * syncDualView declares to the tool that a given DualView
 * has been synced.
 *
 * Arguments:
 *
 * label:     name of the View within the DualView
 * ptr:       that View's data ptr
 * to_device: true if the data is being synchronized to the device
 * 		false otherwise
 */
void syncDualView(const std::string& label, const void* const ptr,
                  bool to_device);
/**
 * modifyDualView declares to the tool that a given DualView
 * has been modified. Note: this means that somebody *called*
 * modify on the DualView, this doesn't get called any time
 * somebody touches the data
 *
 * Arguments:
 *
 * label:     name of the View within the DualView
 * ptr:       that View's data ptr
 * on_device: true if the data is being modified on the device
 * 		false otherwise
 */
void modifyDualView(const std::string& label, const void* const ptr,
                    bool on_device);

void declareMetadata(const std::string& key, const std::string& value);
void initialize(
    const std::string& = {});  // should rename to impl_initialize ASAP
void initialize(const Kokkos::Tools::InitArguments&);
void initialize(int argc, char* argv[]);
void finalize();
bool printHelp(const std::string&);
void parseArgs(const std::string&);

Kokkos_Profiling_SpaceHandle make_space_handle(const char* space_name);

namespace Experimental {

namespace Impl {
struct DirectFenceIDHandle {
  uint32_t value;
};
//
template <typename Space>
uint32_t idForInstance(const uintptr_t instance) {
  static std::mutex instance_mutex;
  const std::lock_guard<std::mutex> lock(instance_mutex);
  /** Needed to be a ptr due to initialization order problems*/
  using map_type = std::map<uintptr_t, uint32_t>;

  static std::shared_ptr<map_type> map;
  if (map.get() == nullptr) {
    map = std::make_shared<map_type>(map_type());
  }

  static uint32_t value = 0;
  constexpr const uint32_t offset =
      Kokkos::Tools::Experimental::NumReservedDeviceIDs;

  auto find = map->find(instance);
  if (find == map->end()) {
    auto ret         = offset + value++;
    (*map)[instance] = ret;
    return ret;
  }

  return find->second;
}

template <typename Space, typename FencingFunctor>
void profile_fence_event(const std::string& name, DirectFenceIDHandle devIDTag,
                         const FencingFunctor& func) {
  uint64_t handle = 0;
  Kokkos::Tools::beginFence(
      name,
      Kokkos::Tools::Experimental::device_id_root<Space>() + devIDTag.value,
      &handle);
  func();
  Kokkos::Tools::endFence(handle);
}

inline uint32_t int_for_synchronization_reason(
    Kokkos::Tools::Experimental::SpecialSynchronizationCases reason) {
  switch (reason) {
    case GlobalDeviceSynchronization: return 0;
    case DeepCopyResourceSynchronization: return 0x00ffffff;
  }
  return 0;
}

template <typename Space, typename FencingFunctor>
void profile_fence_event(
    const std::string& name,
    Kokkos::Tools::Experimental::SpecialSynchronizationCases reason,
    const FencingFunctor& func) {
  uint64_t handle = 0;
  Kokkos::Tools::beginFence(
      name, device_id_root<Space>() + int_for_synchronization_reason(reason),
      &handle);  // TODO: correct ID
  func();
  Kokkos::Tools::endFence(handle);
}
}  // namespace Impl
void set_init_callback(initFunction callback);
void set_finalize_callback(finalizeFunction callback);
void set_parse_args_callback(parseArgsFunction callback);
void set_print_help_callback(printHelpFunction callback);
void set_begin_parallel_for_callback(beginFunction callback);
void set_end_parallel_for_callback(endFunction callback);
void set_begin_parallel_reduce_callback(beginFunction callback);
void set_end_parallel_reduce_callback(endFunction callback);
void set_begin_parallel_scan_callback(beginFunction callback);
void set_end_parallel_scan_callback(endFunction callback);
void set_push_region_callback(pushFunction callback);
void set_pop_region_callback(popFunction callback);
void set_allocate_data_callback(allocateDataFunction callback);
void set_deallocate_data_callback(deallocateDataFunction callback);
void set_create_profile_section_callback(createProfileSectionFunction callback);
void set_start_profile_section_callback(startProfileSectionFunction callback);
void set_stop_profile_section_callback(stopProfileSectionFunction callback);
void set_destroy_profile_section_callback(
    destroyProfileSectionFunction callback);
void set_profile_event_callback(profileEventFunction callback);
void set_begin_deep_copy_callback(beginDeepCopyFunction callback);
void set_end_deep_copy_callback(endDeepCopyFunction callback);
void set_begin_fence_callback(beginFenceFunction callback);
void set_end_fence_callback(endFenceFunction callback);
void set_dual_view_sync_callback(dualViewSyncFunction callback);
void set_dual_view_modify_callback(dualViewModifyFunction callback);
void set_declare_metadata_callback(declareMetadataFunction callback);
void set_request_tool_settings_callback(requestToolSettingsFunction callback);
void set_provide_tool_programming_interface_callback(
    provideToolProgrammingInterfaceFunction callback);
void set_declare_output_type_callback(outputTypeDeclarationFunction callback);
void set_declare_input_type_callback(inputTypeDeclarationFunction callback);
void set_request_output_values_callback(requestValueFunction callback);
void set_declare_optimization_goal_callback(
    optimizationGoalDeclarationFunction callback);
void set_end_context_callback(contextEndFunction callback);
void set_begin_context_callback(contextBeginFunction callback);

void pause_tools();
void resume_tools();

EventSet get_callbacks();
void set_callbacks(EventSet new_events);
}  // namespace Experimental

namespace Experimental {
// forward declarations
size_t get_new_context_id();
size_t get_current_context_id();
}  // namespace Experimental

}  // namespace Tools
namespace Profiling {

bool profileLibraryLoaded();

void beginParallelFor(const std::string& kernelPrefix, const uint32_t devID,
                      uint64_t* kernelID);
void beginParallelReduce(const std::string& kernelPrefix, const uint32_t devID,
                         uint64_t* kernelID);
void beginParallelScan(const std::string& kernelPrefix, const uint32_t devID,
                       uint64_t* kernelID);
void endParallelFor(const uint64_t kernelID);
void endParallelReduce(const uint64_t kernelID);
void endParallelScan(const uint64_t kernelID);
void pushRegion(const std::string& kName);
void popRegion();

void createProfileSection(const std::string& sectionName, uint32_t* secID);
void destroyProfileSection(const uint32_t secID);
void startSection(const uint32_t secID);

void stopSection(const uint32_t secID);

void markEvent(const std::string& eventName);
void allocateData(const SpaceHandle handle, const std::string name,
                  const void* data, const uint64_t size);
void deallocateData(const SpaceHandle space, const std::string label,
                    const void* ptr, const uint64_t size);
void beginDeepCopy(const SpaceHandle dst_space, const std::string dst_label,
                   const void* dst_ptr, const SpaceHandle src_space,
                   const std::string src_label, const void* src_ptr,
                   const uint64_t size);
void endDeepCopy();
void finalize();
void initialize(const std::string& = {});

SpaceHandle make_space_handle(const char* space_name);

namespace Experimental {
using Kokkos::Tools::Experimental::set_allocate_data_callback;
using Kokkos::Tools::Experimental::set_begin_deep_copy_callback;
using Kokkos::Tools::Experimental::set_begin_parallel_for_callback;
using Kokkos::Tools::Experimental::set_begin_parallel_reduce_callback;
using Kokkos::Tools::Experimental::set_begin_parallel_scan_callback;
using Kokkos::Tools::Experimental::set_create_profile_section_callback;
using Kokkos::Tools::Experimental::set_deallocate_data_callback;
using Kokkos::Tools::Experimental::set_destroy_profile_section_callback;
using Kokkos::Tools::Experimental::set_end_deep_copy_callback;
using Kokkos::Tools::Experimental::set_end_parallel_for_callback;
using Kokkos::Tools::Experimental::set_end_parallel_reduce_callback;
using Kokkos::Tools::Experimental::set_end_parallel_scan_callback;
using Kokkos::Tools::Experimental::set_finalize_callback;
using Kokkos::Tools::Experimental::set_init_callback;
using Kokkos::Tools::Experimental::set_parse_args_callback;
using Kokkos::Tools::Experimental::set_pop_region_callback;
using Kokkos::Tools::Experimental::set_print_help_callback;
using Kokkos::Tools::Experimental::set_profile_event_callback;
using Kokkos::Tools::Experimental::set_push_region_callback;
using Kokkos::Tools::Experimental::set_start_profile_section_callback;
using Kokkos::Tools::Experimental::set_stop_profile_section_callback;

using Kokkos::Tools::Experimental::EventSet;

using Kokkos::Tools::Experimental::pause_tools;
using Kokkos::Tools::Experimental::resume_tools;

using Kokkos::Tools::Experimental::get_callbacks;
using Kokkos::Tools::Experimental::set_callbacks;

}  // namespace Experimental
}  // namespace Profiling

namespace Tools {
namespace Experimental {

VariableValue make_variable_value(size_t id, int64_t val);
VariableValue make_variable_value(size_t id, double val);
VariableValue make_variable_value(size_t id, const std::string& val);

SetOrRange make_candidate_set(size_t size, std::string* data);
SetOrRange make_candidate_set(size_t size, int64_t* data);
SetOrRange make_candidate_set(size_t size, double* data);
SetOrRange make_candidate_range(double lower, double upper, double step,
                                bool openLower, bool openUpper);

SetOrRange make_candidate_range(int64_t lower, int64_t upper, int64_t step,
                                bool openLower, bool openUpper);

void declare_optimization_goal(const size_t context,
                               const OptimizationGoal& goal);

size_t declare_output_type(const std::string& typeName, VariableInfo info);

size_t declare_input_type(const std::string& typeName, VariableInfo info);

void set_input_values(size_t contextId, size_t count, VariableValue* values);

void end_context(size_t contextId);
void begin_context(size_t contextId);

void request_output_values(size_t contextId, size_t count,
                           VariableValue* values);

bool have_tuning_tool();

size_t get_new_context_id();
size_t get_current_context_id();

size_t get_new_variable_id();
}  // namespace Experimental
}  // namespace Tools

}  // namespace Kokkos

#endif
