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

#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_ExecPolicy.hpp>
#include <Kokkos_Macros.hpp>
#include <Kokkos_Tuners.hpp>
#include <impl/Kokkos_Profiling_Interface.hpp>
#include <map>
#include <string>
#include <type_traits>
namespace Kokkos {

// forward declaration
bool tune_internals() noexcept;

namespace Tools {

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
void initialize(const std::string& = {});
void finalize();
bool printHelp(const std::string&);
void parseArgs(const std::string&);

Kokkos_Profiling_SpaceHandle make_space_handle(const char* space_name);

namespace Experimental {

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

namespace Impl {

static std::map<std::string, Kokkos::Tools::Experimental::TeamSizeTuner>
    team_tuners;

template <int Rank>
using MDRangeTuningMap =
    std::map<std::string, Kokkos::Tools::Experimental::MDRangeTuner<Rank>>;

template <int Rank>
static MDRangeTuningMap<Rank> mdrange_tuners;

// For any policies without a tuning implementation, with a reducer
template <class ReducerType, class ExecPolicy, class Functor, typename TagType>
void tune_policy(const size_t, const std::string&, ExecPolicy&, const Functor&,
                 TagType) {}

// For any policies without a tuning implementation, without a reducer
template <class ExecPolicy, class Functor, typename TagType>
void tune_policy(const size_t, const std::string&, ExecPolicy&, const Functor&,
                 const TagType&) {}

/**
 * Tuning for parallel_fors and parallel_scans is a fairly simple process.
 *
 * Tuning for a parallel_reduce turns out to be a little more complicated.
 *
 * If you're tuning a reducer, it might be a complex or a simple reducer
 * (an example of simple would be one where the join is just "+".
 *
 * Unfortunately these two paths are very different in terms of which classes
 * get instantiated. Thankfully, all of this complexity is encoded in the
 * ReducerType. If it's a "simple" reducer, this will be Kokkos::InvalidType,
 * otherwise it'll be something else.
 *
 * If the type is complex, for the code to be generally right you _must_
 * pass an instance of that ReducerType to functions that determine
 * eligible team sizes. If the type is simple, you can't construct one,
 * you use the simpler 2-arg formulation of team_size_recommended/max.
 */

namespace Impl {

struct SimpleTeamSizeCalculator {
  template <typename Policy, typename Functor, typename Tag>
  int get_max_team_size(const Policy& policy, const Functor& functor,
                        const Tag tag) {
    auto max = policy.team_size_max(functor, tag);
    return max;
  }
  template <typename Policy, typename Functor, typename Tag>
  int get_recommended_team_size(const Policy& policy, const Functor& functor,
                                const Tag tag) {
    auto max = policy.team_size_recommended(functor, tag);
    return max;
  }
  template <typename Policy, typename Functor>
  int get_mdrange_max_tile_size_product(const Policy& policy,
                                        const Functor& functor,
                                        const Kokkos::ParallelForTag&) {
    using exec_space = typename Policy::execution_space;
    using driver     = Kokkos::Impl::ParallelFor<Functor, Policy, exec_space>;
    return driver::max_tile_size_product(policy, functor);
  }
  template <typename Policy, typename Functor>
  int get_mdrange_max_tile_size_product(const Policy& policy,
                                        const Functor& functor,
                                        const Kokkos::ParallelReduceTag&) {
    using exec_space = typename Policy::execution_space;
    using driver =
        Kokkos::Impl::ParallelReduce<Functor, Policy, Kokkos::InvalidType,
                                     exec_space>;
    return driver::max_tile_size_product(policy, functor);
  }
};

// when we have a complex reducer, we need to pass an
// instance to team_size_recommended/max. Reducers
// aren't default constructible, but they are
// constructible from a reference to an
// instance of their value_type so we construct
// a value_type and temporary reducer here
template <typename ReducerType>
struct ComplexReducerSizeCalculator {
  template <typename Policy, typename Functor, typename Tag>
  int get_max_team_size(const Policy& policy, const Functor& functor,
                        const Tag tag) {
    using value_type = typename ReducerType::value_type;
    value_type value;
    ReducerType reducer_example = ReducerType(value);
    return policy.team_size_max(functor, reducer_example, tag);
  }
  template <typename Policy, typename Functor, typename Tag>
  int get_recommended_team_size(const Policy& policy, const Functor& functor,
                                const Tag tag) {
    using value_type = typename ReducerType::value_type;
    value_type value;
    ReducerType reducer_example = ReducerType(value);
    return policy.team_size_recommended(functor, reducer_example, tag);
  }
  template <typename Policy, typename Functor>
  int get_mdrange_max_tile_size_product(const Policy& policy,
                                        const Functor& functor,
                                        const Kokkos::ParallelReduceTag&) {
    using exec_space = typename Policy::execution_space;
    using driver =
        Kokkos::Impl::ParallelReduce<Functor, Policy, ReducerType, exec_space>;
    return driver::max_tile_size_product(policy, functor);
  }
};

}  // namespace Impl

template <class Tuner, class Functor, class TagType,
          class TuningPermissionFunctor, class Map, class Policy>
void generic_tune_policy(const std::string& label_in, Map& map, Policy& policy,
                         const Functor& functor, const TagType& tag,
                         const TuningPermissionFunctor& should_tune) {
  if (should_tune(policy)) {
    std::string label = label_in;
    if (label_in.empty()) {
      using policy_type =
          typename std::remove_reference<decltype(policy)>::type;
      using work_tag = typename policy_type::work_tag;
      Kokkos::Impl::ParallelConstructName<Functor, work_tag> name(label);
      label = name.get();
    }
    auto tuner_iter = [&]() {
      auto my_tuner = map.find(label);
      if (my_tuner == map.end()) {
        return (map.emplace(label, Tuner(label, policy, functor, tag,
                                         Impl::SimpleTeamSizeCalculator{}))
                    .first);
      }
      return my_tuner;
    }();
    tuner_iter->second.tune(policy);
  }
}
template <class Tuner, class ReducerType, class Functor, class TagType,
          class TuningPermissionFunctor, class Map, class Policy>
void generic_tune_policy(const std::string& label_in, Map& map, Policy& policy,
                         const Functor& functor, const TagType& tag,
                         const TuningPermissionFunctor& should_tune) {
  if (should_tune(policy)) {
    std::string label = label_in;
    if (label_in.empty()) {
      using policy_type =
          typename std::remove_reference<decltype(policy)>::type;
      using work_tag = typename policy_type::work_tag;
      Kokkos::Impl::ParallelConstructName<Functor, work_tag> name(label);
      label = name.get();
    }
    auto tuner_iter = [&]() {
      auto my_tuner = map.find(label);
      if (my_tuner == map.end()) {
        return (map.emplace(
                       label,
                       Tuner(label, policy, functor, tag,
                             Impl::ComplexReducerSizeCalculator<ReducerType>{}))
                    .first);
      }
      return my_tuner;
    }();
    tuner_iter->second.tune(policy);
  }
}

// tune a TeamPolicy, without reducer
template <class Functor, class TagType, class... Properties>
void tune_policy(const size_t /**tuning_context*/, const std::string& label_in,
                 Kokkos::TeamPolicy<Properties...>& policy,
                 const Functor& functor, const TagType& tag) {
  generic_tune_policy<Experimental::TeamSizeTuner>(
      label_in, team_tuners, policy, functor, tag,
      [](const Kokkos::TeamPolicy<Properties...>& candidate_policy) {
        return (candidate_policy.impl_auto_team_size() ||
                candidate_policy.impl_auto_vector_length());
      });
}

// tune a TeamPolicy, with reducer
template <class ReducerType, class Functor, class TagType, class... Properties>
void tune_policy(const size_t /**tuning_context*/, const std::string& label_in,
                 Kokkos::TeamPolicy<Properties...>& policy,
                 const Functor& functor, const TagType& tag) {
  generic_tune_policy<Experimental::TeamSizeTuner, ReducerType>(
      label_in, team_tuners, policy, functor, tag,
      [](const Kokkos::TeamPolicy<Properties...>& candidate_policy) {
        return (candidate_policy.impl_auto_team_size() ||
                candidate_policy.impl_auto_vector_length());
      });
}

// tune a MDRangePolicy, without reducer
template <class Functor, class TagType, class... Properties>
void tune_policy(const size_t /**tuning_context*/, const std::string& label_in,
                 Kokkos::MDRangePolicy<Properties...>& policy,
                 const Functor& functor, const TagType& tag) {
  using Policy              = Kokkos::MDRangePolicy<Properties...>;
  static constexpr int rank = Policy::rank;
  generic_tune_policy<Experimental::MDRangeTuner<rank>>(
      label_in, mdrange_tuners<rank>, policy, functor, tag,
      [](const Policy& candidate_policy) {
        return candidate_policy.impl_tune_tile_size();
      });
}

// tune a MDRangePolicy, with reducer
template <class ReducerType, class Functor, class TagType, class... Properties>
void tune_policy(const size_t /**tuning_context*/, const std::string& label_in,
                 Kokkos::MDRangePolicy<Properties...>& policy,
                 const Functor& functor, const TagType& tag) {
  using Policy              = Kokkos::MDRangePolicy<Properties...>;
  static constexpr int rank = Policy::rank;
  generic_tune_policy<Experimental::MDRangeTuner<rank>, ReducerType>(
      label_in, mdrange_tuners<rank>, policy, functor, tag,
      [](const Policy& candidate_policy) {
        return candidate_policy.impl_tune_tile_size();
      });
}

template <class ReducerType>
struct ReductionSwitcher {
  template <class Functor, class TagType, class ExecPolicy>
  static void tune(const size_t tuning_context, const std::string& label,
                   ExecPolicy& policy, const Functor& functor,
                   const TagType& tag) {
    if (Kokkos::tune_internals()) {
      tune_policy<ReducerType>(tuning_context, label, policy, functor, tag);
    }
  }
};

template <>
struct ReductionSwitcher<Kokkos::InvalidType> {
  template <class Functor, class TagType, class ExecPolicy>
  static void tune(const size_t tuning_context, const std::string& label,
                   ExecPolicy& policy, const Functor& functor,
                   const TagType& tag) {
    if (Kokkos::tune_internals()) {
      tune_policy(tuning_context, label, policy, functor, tag);
    }
  }
};

template <class Tuner, class Functor, class TagType,
          class TuningPermissionFunctor, class Map, class Policy>
void generic_report_results(const std::string& label_in, Map& map,
                            Policy& policy, const Functor&, const TagType&,
                            const TuningPermissionFunctor& should_tune) {
  if (should_tune(policy)) {
    std::string label = label_in;
    if (label_in.empty()) {
      using policy_type =
          typename std::remove_reference<decltype(policy)>::type;
      using work_tag = typename policy_type::work_tag;
      Kokkos::Impl::ParallelConstructName<Functor, work_tag> name(label);
      label = name.get();
    }
    auto tuner_iter = map[label];
    tuner_iter.end();
  }
}

// report results for a policy type we don't tune (do nothing)
template <class ExecPolicy, class Functor, typename TagType>
void report_policy_results(const size_t, const std::string&, ExecPolicy&,
                           const Functor&, const TagType&) {}

// report results for a TeamPolicy
template <class Functor, class TagType, class... Properties>
void report_policy_results(const size_t /**tuning_context*/,
                           const std::string& label_in,
                           Kokkos::TeamPolicy<Properties...>& policy,
                           const Functor& functor, const TagType& tag) {
  generic_report_results<Experimental::TeamSizeTuner>(
      label_in, team_tuners, policy, functor, tag,
      [](const Kokkos::TeamPolicy<Properties...>& candidate_policy) {
        return (candidate_policy.impl_auto_team_size() ||
                candidate_policy.impl_auto_vector_length());
      });
}

// report results for an MDRangePolicy
template <class Functor, class TagType, class... Properties>
void report_policy_results(const size_t /**tuning_context*/,
                           const std::string& label_in,
                           Kokkos::MDRangePolicy<Properties...>& policy,
                           const Functor& functor, const TagType& tag) {
  using Policy              = Kokkos::MDRangePolicy<Properties...>;
  static constexpr int rank = Policy::rank;
  generic_report_results<Experimental::MDRangeTuner<rank>>(
      label_in, mdrange_tuners<rank>, policy, functor, tag,
      [](const Policy& candidate_policy) {
        return candidate_policy.impl_tune_tile_size();
      });
}

template <class ExecPolicy, class FunctorType>
void begin_parallel_for(ExecPolicy& policy, FunctorType& functor,
                        const std::string& label, uint64_t& kpID) {
  if (Kokkos::Tools::profileLibraryLoaded()) {
    Kokkos::Impl::ParallelConstructName<FunctorType,
                                        typename ExecPolicy::work_tag>
        name(label);
    Kokkos::Tools::beginParallelFor(
        name.get(), Kokkos::Profiling::Experimental::device_id(policy.space()),
        &kpID);
  }
#ifdef KOKKOS_ENABLE_TUNING
  size_t context_id = Kokkos::Tools::Experimental::get_new_context_id();
  if (Kokkos::tune_internals()) {
    tune_policy(context_id, label, policy, functor, Kokkos::ParallelForTag{});
  }
#else
  (void)functor;
#endif
}

template <class ExecPolicy, class FunctorType>
void end_parallel_for(ExecPolicy& policy, FunctorType& functor,
                      const std::string& label, uint64_t& kpID) {
  if (Kokkos::Tools::profileLibraryLoaded()) {
    Kokkos::Tools::endParallelFor(kpID);
  }
#ifdef KOKKOS_ENABLE_TUNING
  size_t context_id = Kokkos::Tools::Experimental::get_current_context_id();
  if (Kokkos::tune_internals()) {
    report_policy_results(context_id, label, policy, functor,
                          Kokkos::ParallelForTag{});
  }
#else
  (void)policy;
  (void)functor;
  (void)label;
#endif
}

template <class ExecPolicy, class FunctorType>
void begin_parallel_scan(ExecPolicy& policy, FunctorType& functor,
                         const std::string& label, uint64_t& kpID) {
  if (Kokkos::Tools::profileLibraryLoaded()) {
    Kokkos::Impl::ParallelConstructName<FunctorType,
                                        typename ExecPolicy::work_tag>
        name(label);
    Kokkos::Tools::beginParallelScan(
        name.get(), Kokkos::Profiling::Experimental::device_id(policy.space()),
        &kpID);
  }
#ifdef KOKKOS_ENABLE_TUNING
  size_t context_id = Kokkos::Tools::Experimental::get_new_context_id();
  if (Kokkos::tune_internals()) {
    tune_policy(context_id, label, policy, functor, Kokkos::ParallelScanTag{});
  }
#else
  (void)functor;
#endif
}

template <class ExecPolicy, class FunctorType>
void end_parallel_scan(ExecPolicy& policy, FunctorType& functor,
                       const std::string& label, uint64_t& kpID) {
  if (Kokkos::Tools::profileLibraryLoaded()) {
    Kokkos::Tools::endParallelScan(kpID);
  }
#ifdef KOKKOS_ENABLE_TUNING
  size_t context_id = Kokkos::Tools::Experimental::get_current_context_id();
  if (Kokkos::tune_internals()) {
    report_policy_results(context_id, label, policy, functor,
                          Kokkos::ParallelScanTag{});
  }
#else
  (void)policy;
  (void)functor;
  (void)label;
#endif
}

template <class ReducerType, class ExecPolicy, class FunctorType>
void begin_parallel_reduce(ExecPolicy& policy, FunctorType& functor,
                           const std::string& label, uint64_t& kpID) {
  if (Kokkos::Tools::profileLibraryLoaded()) {
    Kokkos::Impl::ParallelConstructName<FunctorType,
                                        typename ExecPolicy::work_tag>
        name(label);
    Kokkos::Tools::beginParallelReduce(
        name.get(), Kokkos::Profiling::Experimental::device_id(policy.space()),
        &kpID);
  }
#ifdef KOKKOS_ENABLE_TUNING
  size_t context_id = Kokkos::Tools::Experimental::get_new_context_id();
  ReductionSwitcher<ReducerType>::tune(context_id, label, policy, functor,
                                       Kokkos::ParallelReduceTag{});
#else
  (void)functor;
#endif
}

template <class ReducerType, class ExecPolicy, class FunctorType>
void end_parallel_reduce(ExecPolicy& policy, FunctorType& functor,
                         const std::string& label, uint64_t& kpID) {
  if (Kokkos::Tools::profileLibraryLoaded()) {
    Kokkos::Tools::endParallelReduce(kpID);
  }
#ifdef KOKKOS_ENABLE_TUNING
  size_t context_id = Kokkos::Tools::Experimental::get_current_context_id();
  if (Kokkos::tune_internals()) {
    report_policy_results(context_id, label, policy, functor,
                          Kokkos::ParallelReduceTag{});
  }
#else
  (void)policy;
  (void)functor;
  (void)label;
#endif
}

}  // namespace Impl

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
