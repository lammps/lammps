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

#ifndef KOKKOS_IMPL_KOKKOS_TOOLS_GENERIC_HPP
#define KOKKOS_IMPL_KOKKOS_TOOLS_GENERIC_HPP

#include <impl/Kokkos_Profiling.hpp>
#include <impl/Kokkos_FunctorAnalysis.hpp>

#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_ExecPolicy.hpp>
#include <Kokkos_Macros.hpp>
#include <Kokkos_Tuners.hpp>

namespace Kokkos {

namespace Tools {

namespace Experimental {

namespace Impl {

static std::map<std::string, Kokkos::Tools::Experimental::TeamSizeTuner>
    team_tuners;

static std::map<std::string,
                Kokkos::Tools::Experimental::RangePolicyOccupancyTuner>
    range_policy_tuners;

template <int Rank>
using MDRangeTuningMap =
    std::map<std::string, Kokkos::Tools::Experimental::MDRangeTuner<Rank>>;

template <int Rank>
static MDRangeTuningMap<Rank> mdrange_tuners;

// For any policies without a tuning implementation, with a reducer
template <class ReducerType, class ExecPolicy, class Functor, typename TagType>
auto tune_policy(const size_t, const std::string&, const ExecPolicy& policy,
                 const Functor&, TagType) {
  return policy;
}

// For any policies without a tuning implementation, without a reducer
template <class ExecPolicy, class Functor, typename TagType>
auto tune_policy(const size_t, const std::string&, const ExecPolicy& policy,
                 const Functor&, const TagType&) {
  return policy;
}

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
  template <typename Policy, typename FunctorReducer>
  int get_max_team_size(const Policy& policy,
                        const FunctorReducer& functor_reducer,
                        const Kokkos::ParallelReduceTag tag) {
    auto max = policy.team_size_max(functor_reducer.get_functor(),
                                    functor_reducer.get_reducer(), tag);
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
  template <typename Policy, typename FunctorReducer>
  int get_mdrange_max_tile_size_product(const Policy& policy,
                                        const FunctorReducer& functor_reducer,
                                        const Kokkos::ParallelReduceTag&) {
    using exec_space = typename Policy::execution_space;
    using driver =
        Kokkos::Impl::ParallelReduce<FunctorReducer, Policy, exec_space>;
    return driver::max_tile_size_product(policy, functor_reducer.get_functor());
  }
};

// when we have a complex reducer, we need to pass an
// instance to team_size_recommended/max. Reducers
// aren't default constructible, but they are
// constructible from a reference to an
// instance of their value_type so we construct
// a value_type and temporary reducer here
struct ComplexReducerSizeCalculator {
  template <typename Policy, typename FunctorReducer, typename Tag>
  int get_max_team_size(const Policy& policy,
                        const FunctorReducer& functor_reducer, const Tag tag) {
    return policy.team_size_max(functor_reducer.get_functor(),
                                functor_reducer.get_reducer(), tag);
  }
  template <typename Policy, typename FunctorReducer, typename Tag>
  int get_recommended_team_size(const Policy& policy,
                                const FunctorReducer& functor_reducer,
                                const Tag tag) {
    return policy.team_size_recommended(functor_reducer.get_functor(),
                                        functor_reducer.get_reducer(), tag);
  }
  template <typename Policy, typename FunctorReducer>
  int get_mdrange_max_tile_size_product(const Policy& policy,
                                        const FunctorReducer& functor_reducer,
                                        const Kokkos::ParallelReduceTag&) {
    using exec_space = typename Policy::execution_space;
    using driver =
        Kokkos::Impl::ParallelReduce<FunctorReducer, Policy, exec_space>;
    return driver::max_tile_size_product(policy, functor_reducer.get_functor());
  }
};

template <typename Policy>
auto default_tuned_version_of(const Policy& policy) {
  return policy;
}

}  // namespace Impl

template <class Tuner, class Functor, class TagType,
          class TuningPermissionFunctor, class Map, class Policy>
auto generic_tune_policy(const std::string& label_in, Map& map,
                         const Policy& policy, const Functor& functor,
                         const TagType& tag,
                         const TuningPermissionFunctor& should_tune) {
  if (should_tune(policy)) {
    std::string label = label_in;
    if (label_in.empty()) {
      using policy_type = std::remove_reference_t<decltype(policy)>;
      using work_tag    = typename policy_type::work_tag;
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
    return tuner_iter->second.tune(policy);
  }
  return Impl::default_tuned_version_of(policy);
}
template <class Tuner, class ReducerType, class Functor, class TagType,
          class TuningPermissionFunctor, class Map, class Policy>
auto generic_tune_policy(const std::string& label_in, Map& map,
                         const Policy& policy, const Functor& functor,
                         const TagType& tag,
                         const TuningPermissionFunctor& should_tune) {
  if (should_tune(policy)) {
    std::string label = label_in;
    if (label_in.empty()) {
      using policy_type = std::remove_reference_t<decltype(policy)>;
      using work_tag    = typename policy_type::work_tag;
      Kokkos::Impl::ParallelConstructName<Functor, work_tag> name(label);
      label = name.get();
    }
    auto tuner_iter = [&]() {
      auto my_tuner = map.find(label);
      if (my_tuner == map.end()) {
        return (map.emplace(label, Tuner(label, policy, functor, tag,
                                         Impl::ComplexReducerSizeCalculator{}))
                    .first);
      }
      return my_tuner;
    }();
    return tuner_iter->second.tune(policy);
  }
  return Impl::default_tuned_version_of(policy);
}

// tune a TeamPolicy, without reducer
template <class Functor, class TagType, class... Properties>
auto tune_policy(const size_t /**tuning_context*/, const std::string& label_in,
                 const Kokkos::TeamPolicy<Properties...>& policy,
                 const Functor& functor, const TagType& tag) {
  return generic_tune_policy<Experimental::TeamSizeTuner>(
      label_in, team_tuners, policy, functor, tag,
      [](const Kokkos::TeamPolicy<Properties...>& candidate_policy) {
        return (candidate_policy.impl_auto_team_size() ||
                candidate_policy.impl_auto_vector_length());
      });
}

// tune a TeamPolicy, with reducer
template <class ReducerType, class Functor, class TagType, class... Properties>
auto tune_policy(const size_t /**tuning_context*/, const std::string& label_in,
                 const Kokkos::TeamPolicy<Properties...>& policy,
                 const Functor& functor, const TagType& tag) {
  return generic_tune_policy<Experimental::TeamSizeTuner, ReducerType>(
      label_in, team_tuners, policy, functor, tag,
      [](const Kokkos::TeamPolicy<Properties...>& candidate_policy) {
        return (candidate_policy.impl_auto_team_size() ||
                candidate_policy.impl_auto_vector_length());
      });
}

template <class Functor, class TagType, class... Properties>
auto tune_occupancy_controlled_policy(
    const size_t /**tuning_context*/, const std::string& label_in,
    const Kokkos::RangePolicy<Properties...>& policy, const Functor& functor,
    const TagType& tag) {
  return generic_tune_policy<Experimental::RangePolicyOccupancyTuner>(
      label_in, range_policy_tuners, policy, functor, tag,
      [](const Kokkos::RangePolicy<Properties...>& candidate_policy) {
        return candidate_policy.impl_get_occupancy_control().should_tune();
      });
}
template <class Functor, class TagType, class... Properties>
auto tune_range_policy(const size_t tuning_context, const std::string& label_in,
                       const Kokkos::RangePolicy<Properties...>& policy,
                       const Functor& functor, const TagType& tag,
                       std::true_type) {
  return tune_occupancy_controlled_policy(tuning_context, label_in, policy,
                                          functor, tag);
}
template <class Functor, class TagType, class... Properties>
auto tune_range_policy(const size_t /**tuning_context*/,
                       const std::string& /*label_in*/,
                       const Kokkos::RangePolicy<Properties...>& policy,
                       const Functor& /**functor*/, const TagType& /**tag*/,
                       std::false_type) {
  return policy;
}

// Reducer versions
template <class RT, class Functor, class TagType, class... Properties>
auto tune_occupancy_controlled_policy(
    const size_t /**tuning_context*/, const std::string& label_in,
    const Kokkos::RangePolicy<Properties...>& policy, const Functor& functor,
    const TagType& tag) {
  return generic_tune_policy<Experimental::RangePolicyOccupancyTuner>(
      label_in, range_policy_tuners, policy, functor, tag,
      [](const Kokkos::RangePolicy<Properties...>& candidate_policy) {
        return candidate_policy.impl_get_occupancy_control().should_tune();
      });
}
template <class RT, class Functor, class TagType, class... Properties>
auto tune_range_policy(const size_t tuning_context, const std::string& label_in,
                       const Kokkos::RangePolicy<Properties...>& policy,
                       const Functor& functor, const TagType& tag,
                       std::true_type) {
  return tune_occupancy_controlled_policy<RT>(tuning_context, label_in, policy,
                                              functor, tag);
}
template <class ReducerType, class Functor, class TagType, class... Properties>
auto tune_range_policy(const size_t /**tuning_context*/,
                       const std::string& /**label_in*/,
                       const Kokkos::RangePolicy<Properties...>& policy,
                       const Functor& /**functor*/, const TagType& /**tag*/,
                       std::false_type) {
  return policy;
}

// tune a RangePolicy, without reducer
template <class Functor, class TagType, class... Properties>
auto tune_policy(const size_t tuning_context, const std::string& label_in,
                 const Kokkos::RangePolicy<Properties...>& policy,
                 const Functor& functor, const TagType& tag) {
  using policy_t = Kokkos::RangePolicy<Properties...>;
  using has_desired_occupancy =
      typename std::is_same<typename policy_t::occupancy_control,
                            Kokkos::Experimental::DesiredOccupancy>::type;
  return tune_range_policy(tuning_context, label_in, policy, functor, tag,
                           has_desired_occupancy{});
}

// tune a RangePolicy, with reducer
template <class ReducerType, class Functor, class TagType, class... Properties>
auto tune_policy(const size_t tuning_context, const std::string& label_in,
                 const Kokkos::RangePolicy<Properties...>& policy,
                 const Functor& functor, const TagType& tag) {
  using policy_t = Kokkos::RangePolicy<Properties...>;
  using has_desired_occupancy =
      typename std::is_same<typename policy_t::occupancy_control,
                            Kokkos::Experimental::DesiredOccupancy>::type;
  return tune_range_policy<ReducerType>(tuning_context, label_in, policy,
                                        functor, tag, has_desired_occupancy{});
}

// tune a MDRangePolicy, without reducer
template <class Functor, class TagType, class... Properties>
auto tune_policy(const size_t /**tuning_context*/, const std::string& label_in,
                 const Kokkos::MDRangePolicy<Properties...>& policy,
                 const Functor& functor, const TagType& tag) {
  using Policy              = Kokkos::MDRangePolicy<Properties...>;
  static constexpr int rank = Policy::rank;
  return generic_tune_policy<Experimental::MDRangeTuner<rank>>(
      label_in, mdrange_tuners<rank>, policy, functor, tag,
      [](const Policy& candidate_policy) {
        return candidate_policy.impl_tune_tile_size();
      });
}

// tune a MDRangePolicy, with reducer
template <class ReducerType, class Functor, class TagType, class... Properties>
auto tune_policy(const size_t /**tuning_context*/, const std::string& label_in,
                 const Kokkos::MDRangePolicy<Properties...>& policy,
                 const Functor& functor, const TagType& tag) {
  using Policy              = Kokkos::MDRangePolicy<Properties...>;
  static constexpr int rank = Policy::rank;
  return generic_tune_policy<Experimental::MDRangeTuner<rank>, ReducerType>(
      label_in, mdrange_tuners<rank>, policy, functor, tag,
      [](const Policy& candidate_policy) {
        return candidate_policy.impl_tune_tile_size();
      });
}

template <class ReducerType>
struct ReductionSwitcher {
  template <class Functor, class TagType, class ExecPolicy>
  static auto tune(const size_t tuning_context, const std::string& label,
                   const ExecPolicy& policy, const Functor& functor,
                   const TagType& tag) {
    if (Kokkos::tune_internals()) {
      return tune_policy<ReducerType>(tuning_context, label, policy, functor,
                                      tag);
    }
    return Impl::default_tuned_version_of(policy);
  }
};

template <>
struct ReductionSwitcher<Kokkos::InvalidType> {
  template <class Functor, class TagType, class ExecPolicy>
  static auto tune(const size_t tuning_context, const std::string& label,
                   const ExecPolicy& policy, const Functor& functor,
                   const TagType& tag) {
    if (Kokkos::tune_internals()) {
      return tune_policy(tuning_context, label, policy, functor, tag);
    }
    return Impl::default_tuned_version_of(policy);
  }
};

template <class Tuner, class Functor, class TagType,
          class TuningPermissionFunctor, class Map, class Policy>
void generic_report_results(const std::string& label_in, Map& map,
                            const Policy& policy, const Functor&,
                            const TagType&,
                            const TuningPermissionFunctor& should_tune) {
  if (should_tune(policy)) {
    std::string label = label_in;
    if (label_in.empty()) {
      using policy_type = std::remove_reference_t<decltype(policy)>;
      using work_tag    = typename policy_type::work_tag;
      Kokkos::Impl::ParallelConstructName<Functor, work_tag> name(label);
      label = name.get();
    }
    auto tuner_iter = map[label];
    tuner_iter.end();
  }
}

// report results for a policy type we don't tune (do nothing)
template <class ExecPolicy, class Functor, typename TagType>
void report_policy_results(const size_t, const std::string&, const ExecPolicy&,
                           const Functor&, const TagType&) {}

// report results for a TeamPolicy
template <class Functor, class TagType, class... Properties>
void report_policy_results(const size_t /**tuning_context*/,
                           const std::string& label_in,
                           const Kokkos::TeamPolicy<Properties...>& policy,
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
                           const Kokkos::MDRangePolicy<Properties...>& policy,
                           const Functor& functor, const TagType& tag) {
  using Policy              = Kokkos::MDRangePolicy<Properties...>;
  static constexpr int rank = Policy::rank;
  generic_report_results<Experimental::MDRangeTuner<rank>>(
      label_in, mdrange_tuners<rank>, policy, functor, tag,
      [](const Policy& candidate_policy) {
        return candidate_policy.impl_tune_tile_size();
      });
}

// report results for an MDRangePolicy
template <class Functor, class TagType, class... Properties>
void report_policy_results(const size_t /**tuning_context*/,
                           const std::string& label_in,
                           const Kokkos::RangePolicy<Properties...>& policy,
                           const Functor& functor, const TagType& tag) {
  using Policy = Kokkos::RangePolicy<Properties...>;
  generic_report_results<Experimental::RangePolicyOccupancyTuner>(
      label_in, range_policy_tuners, policy, functor, tag, [](const Policy&) {
        return Kokkos::RangePolicy<
            Properties...>::traits::experimental_contains_desired_occupancy;
      });
}

}  // namespace Impl

}  // namespace Experimental

namespace Impl {

template <class ExecPolicy, class FunctorType>
auto begin_parallel_for(const ExecPolicy& policy, FunctorType& functor,
                        const std::string& label, uint64_t& kpID) {
  using response_type =
      Kokkos::Tools::Impl::ToolResponse<ExecPolicy, FunctorType>;
  response_type response{policy};
  if (Kokkos::Tools::profileLibraryLoaded()) {
    Kokkos::Impl::ParallelConstructName<FunctorType,
                                        typename ExecPolicy::work_tag>
        name(label);
    Kokkos::Tools::beginParallelFor(
        name.get(), Kokkos::Profiling::Experimental::device_id(policy.space()),
        &kpID);
  }
#ifdef KOKKOS_ENABLE_TUNING
  size_t context_id = Kokkos::Tools::Experimental::get_current_context_id();
  if (Kokkos::tune_internals()) {
    return response_type{Kokkos::Tools::Experimental::Impl::tune_policy(
        context_id, label, policy, functor, Kokkos::ParallelForTag{})};
  }
#else
  (void)functor;
#endif
  return response;
}

template <class ExecPolicy, class FunctorType>
void end_parallel_for(const ExecPolicy& policy, FunctorType& functor,
                      const std::string& label, uint64_t& kpID) {
  if (Kokkos::Tools::profileLibraryLoaded()) {
    Kokkos::Tools::endParallelFor(kpID);
  }
#ifdef KOKKOS_ENABLE_TUNING
  size_t context_id = Kokkos::Tools::Experimental::get_current_context_id();
  if (Kokkos::tune_internals()) {
    Experimental::Impl::report_policy_results(
        context_id, label, policy, functor, Kokkos::ParallelForTag{});
  }
#else
  (void)policy;
  (void)functor;
  (void)label;
#endif
}

template <class ExecPolicy, class FunctorType>
auto begin_parallel_scan(const ExecPolicy& policy, FunctorType& functor,
                         const std::string& label, uint64_t& kpID) {
  using response_type =
      Kokkos::Tools::Impl::ToolResponse<ExecPolicy, FunctorType>;
  response_type response{policy};
  if (Kokkos::Tools::profileLibraryLoaded()) {
    Kokkos::Impl::ParallelConstructName<FunctorType,
                                        typename ExecPolicy::work_tag>
        name(label);
    Kokkos::Tools::beginParallelScan(
        name.get(), Kokkos::Profiling::Experimental::device_id(policy.space()),
        &kpID);
  }
#ifdef KOKKOS_ENABLE_TUNING
  size_t context_id = Kokkos::Tools::Experimental::get_current_context_id();
  if (Kokkos::tune_internals()) {
    return response_type{Kokkos::Tools::Experimental::Impl::tune_policy(
        context_id, label, policy, functor, Kokkos::ParallelScanTag{})};
  }
#else
  (void)functor;
#endif
  return response;
}

template <class ExecPolicy, class FunctorType>
void end_parallel_scan(const ExecPolicy& policy, FunctorType& functor,
                       const std::string& label, uint64_t& kpID) {
  if (Kokkos::Tools::profileLibraryLoaded()) {
    Kokkos::Tools::endParallelScan(kpID);
  }
#ifdef KOKKOS_ENABLE_TUNING
  size_t context_id = Kokkos::Tools::Experimental::get_current_context_id();
  if (Kokkos::tune_internals()) {
    Experimental::Impl::report_policy_results(
        context_id, label, policy, functor, Kokkos::ParallelScanTag{});
  }
#else
  (void)policy;
  (void)functor;
  (void)label;
#endif
}

template <class ReducerType, class ExecPolicy, class FunctorType>
auto begin_parallel_reduce(const ExecPolicy& policy, FunctorType& functor,
                           const std::string& label, uint64_t& kpID) {
  using response_type = ToolResponse<ExecPolicy, FunctorType>;
  response_type response{policy};
  if (Kokkos::Tools::profileLibraryLoaded()) {
    Kokkos::Impl::ParallelConstructName<FunctorType,
                                        typename ExecPolicy::work_tag>
        name(label);
    Kokkos::Tools::beginParallelReduce(
        name.get(), Kokkos::Profiling::Experimental::device_id(policy.space()),
        &kpID);
  }
#ifdef KOKKOS_ENABLE_TUNING
  size_t context_id = Kokkos::Tools::Experimental::get_current_context_id();
  return response_type{Experimental::Impl::ReductionSwitcher<ReducerType>::tune(
      context_id, label, policy, functor, Kokkos::ParallelReduceTag{})};
#else
  (void)functor;
#endif
  return response;
}

template <class ReducerType, class ExecPolicy, class FunctorType>
void end_parallel_reduce(const ExecPolicy& policy, FunctorType& functor,
                         const std::string& label, uint64_t& kpID) {
  if (Kokkos::Tools::profileLibraryLoaded()) {
    Kokkos::Tools::endParallelReduce(kpID);
  }
#ifdef KOKKOS_ENABLE_TUNING
  size_t context_id = Kokkos::Tools::Experimental::get_current_context_id();
  if (Kokkos::tune_internals()) {
    Experimental::Impl::report_policy_results(
        context_id, label, policy, functor, Kokkos::ParallelReduceTag{});
  }
#else
  (void)policy;
  (void)functor;
  (void)label;
#endif
}

}  // namespace Impl

}  // namespace Tools

}  // namespace Kokkos

#endif  // header guard
