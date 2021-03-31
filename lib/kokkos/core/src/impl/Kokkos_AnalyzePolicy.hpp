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

#ifndef KOKKOS_IMPL_ANALYZE_POLICY_HPP
#define KOKKOS_IMPL_ANALYZE_POLICY_HPP

#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_Concepts.hpp>
#include <impl/Kokkos_Tags.hpp>
#include <impl/Kokkos_GraphImpl_fwd.hpp>
#include <impl/Kokkos_Error.hpp>
#include <impl/Kokkos_EBO.hpp>

namespace Kokkos {
namespace Experimental {
struct DesiredOccupancy {
  int m_occ = 100;
  explicit constexpr DesiredOccupancy(int occ) : m_occ(occ) {
    KOKKOS_EXPECTS(0 <= occ && occ <= 100);
  }
  explicit constexpr operator int() const { return m_occ; }
  constexpr int value() const { return m_occ; }
  explicit DesiredOccupancy() = default;
};
struct MaximizeOccupancy {
  explicit MaximizeOccupancy() = default;
};
}  // namespace Experimental

namespace Impl {
template <typename ExecutionSpace = void, typename Schedule = void,
          typename WorkTag = void, typename IndexType = void,
          typename IterationPattern = void, typename LaunchBounds = void,
          typename MyWorkItemProperty =
              Kokkos::Experimental::WorkItemProperty::None_t,
          typename IsGraphKernel    = std::false_type,
          typename OccupancyControl = Kokkos::Experimental::MaximizeOccupancy>
struct PolicyTraitsBase {
  using type =
      PolicyTraitsBase<ExecutionSpace, Schedule, WorkTag, IndexType,
                       IterationPattern, LaunchBounds, MyWorkItemProperty,
                       IsGraphKernel, OccupancyControl>;

  using execution_space    = ExecutionSpace;
  using schedule_type      = Schedule;
  using work_tag           = WorkTag;
  using index_type         = IndexType;
  using iteration_pattern  = IterationPattern;
  using launch_bounds      = LaunchBounds;
  using work_item_property = MyWorkItemProperty;
  using is_graph_kernel    = IsGraphKernel;
  using occupancy_control  = OccupancyControl;
};

template <typename PolicyBase, typename Property>
struct SetWorkItemProperty {
  static_assert(
      std::is_same<typename PolicyBase::work_item_property,
                   Kokkos::Experimental::WorkItemProperty::None_t>::value,
      "Kokkos Error: More than one work item property given");
  using type = PolicyTraitsBase<
      typename PolicyBase::execution_space, typename PolicyBase::schedule_type,
      typename PolicyBase::work_tag, typename PolicyBase::index_type,
      typename PolicyBase::iteration_pattern,
      typename PolicyBase::launch_bounds, Property,
      typename PolicyBase::is_graph_kernel,
      typename PolicyBase::occupancy_control>;
};

template <typename PolicyBase, typename ExecutionSpace>
struct SetExecutionSpace {
  static_assert(is_void<typename PolicyBase::execution_space>::value,
                "Kokkos Error: More than one execution space given");
  using type =
      PolicyTraitsBase<ExecutionSpace, typename PolicyBase::schedule_type,
                       typename PolicyBase::work_tag,
                       typename PolicyBase::index_type,
                       typename PolicyBase::iteration_pattern,
                       typename PolicyBase::launch_bounds,
                       typename PolicyBase::work_item_property,
                       typename PolicyBase::is_graph_kernel,
                       typename PolicyBase::occupancy_control>;
};

template <typename PolicyBase, typename Schedule>
struct SetSchedule {
  static_assert(is_void<typename PolicyBase::schedule_type>::value,
                "Kokkos Error: More than one schedule type given");
  using type = PolicyTraitsBase<typename PolicyBase::execution_space, Schedule,
                                typename PolicyBase::work_tag,
                                typename PolicyBase::index_type,
                                typename PolicyBase::iteration_pattern,
                                typename PolicyBase::launch_bounds,
                                typename PolicyBase::work_item_property,
                                typename PolicyBase::is_graph_kernel,
                                typename PolicyBase::occupancy_control>;
};

template <typename PolicyBase, typename WorkTag>
struct SetWorkTag {
  static_assert(is_void<typename PolicyBase::work_tag>::value,
                "Kokkos Error: More than one work tag given");
  using type = PolicyTraitsBase<typename PolicyBase::execution_space,
                                typename PolicyBase::schedule_type, WorkTag,
                                typename PolicyBase::index_type,
                                typename PolicyBase::iteration_pattern,
                                typename PolicyBase::launch_bounds,
                                typename PolicyBase::work_item_property,
                                typename PolicyBase::is_graph_kernel,
                                typename PolicyBase::occupancy_control>;
};

template <typename PolicyBase, typename IndexType>
struct SetIndexType {
  static_assert(is_void<typename PolicyBase::index_type>::value,
                "Kokkos Error: More than one index type given");
  using type = PolicyTraitsBase<typename PolicyBase::execution_space,
                                typename PolicyBase::schedule_type,
                                typename PolicyBase::work_tag, IndexType,
                                typename PolicyBase::iteration_pattern,
                                typename PolicyBase::launch_bounds,
                                typename PolicyBase::work_item_property,
                                typename PolicyBase::is_graph_kernel,
                                typename PolicyBase::occupancy_control>;
};

template <typename PolicyBase, typename IterationPattern>
struct SetIterationPattern {
  static_assert(is_void<typename PolicyBase::iteration_pattern>::value,
                "Kokkos Error: More than one iteration_pattern given");
  using type = PolicyTraitsBase<
      typename PolicyBase::execution_space, typename PolicyBase::schedule_type,
      typename PolicyBase::work_tag, typename PolicyBase::index_type,
      IterationPattern, typename PolicyBase::launch_bounds,
      typename PolicyBase::work_item_property,
      typename PolicyBase::is_graph_kernel,
      typename PolicyBase::occupancy_control>;
};

template <typename PolicyBase, typename LaunchBounds>
struct SetLaunchBounds {
  static_assert(is_void<typename PolicyBase::launch_bounds>::value,
                "Kokkos Error: More than one launch_bounds given");
  using type = PolicyTraitsBase<
      typename PolicyBase::execution_space, typename PolicyBase::schedule_type,
      typename PolicyBase::work_tag, typename PolicyBase::index_type,
      typename PolicyBase::iteration_pattern, LaunchBounds,
      typename PolicyBase::work_item_property,
      typename PolicyBase::is_graph_kernel,
      typename PolicyBase::occupancy_control>;
};

template <typename PolicyBase>
struct SetIsGraphKernel {
  using type = PolicyTraitsBase<
      typename PolicyBase::execution_space, typename PolicyBase::schedule_type,
      typename PolicyBase::work_tag, typename PolicyBase::index_type,
      typename PolicyBase::iteration_pattern,
      typename PolicyBase::launch_bounds,
      typename PolicyBase::work_item_property, std::true_type,
      typename PolicyBase::occupancy_control>;
};

template <typename PolicyBase, typename OccupancyControl>
struct SetOccupancyControl {
  using type = PolicyTraitsBase<
      typename PolicyBase::execution_space, typename PolicyBase::schedule_type,
      typename PolicyBase::work_tag, typename PolicyBase::index_type,
      typename PolicyBase::iteration_pattern,
      typename PolicyBase::launch_bounds,
      typename PolicyBase::work_item_property,
      typename PolicyBase::is_graph_kernel, OccupancyControl>;
};

template <typename Base, typename... Traits>
struct AnalyzePolicy;

// TODO DSH rewrite this to be more extensible once we have metaprogramming from
//      desul
template <typename Base, typename T, typename... Traits>
struct AnalyzePolicy<Base, T, Traits...>
    : public AnalyzePolicy<
          typename std::conditional_t<
              is_execution_space<T>::value, SetExecutionSpace<Base, T>,
              std::conditional_t<
                  is_schedule_type<T>::value, SetSchedule<Base, T>,
                  std::conditional_t<
                      is_index_type<T>::value, SetIndexType<Base, T>,
                      std::conditional_t<
                          std::is_integral<T>::value,
                          SetIndexType<Base, IndexType<T>>,
                          std::conditional_t<
                              is_iteration_pattern<T>::value,
                              SetIterationPattern<Base, T>,
                              std::conditional_t<
                                  is_launch_bounds<T>::value,
                                  SetLaunchBounds<Base, T>,
                                  std::conditional_t<
                                      Kokkos::Experimental::
                                          is_work_item_property<T>::value,
                                      SetWorkItemProperty<Base, T>,
                                      std::conditional_t<
                                          std::is_same<T,
                                                       IsGraphKernelTag>::value,
                                          SetIsGraphKernel<Base>,
                                          std::conditional_t<
                                              std::is_same<
                                                  T, Kokkos::Experimental::
                                                         DesiredOccupancy>::
                                                      value ||
                                                  std::is_same<
                                                      T,
                                                      Kokkos::Experimental::
                                                          MaximizeOccupancy>::
                                                      value,
                                              SetOccupancyControl<Base, T>,
                                              std::conditional_t<
                                                  !std::is_void<T>::value,
                                                  SetWorkTag<Base, T>,
                                                  Base>>>>>>>>>>::type,
          Traits...> {};

template <typename Base>
struct AnalyzePolicy<Base> {
  static constexpr auto execution_space_is_defaulted =
      std::is_void<typename Base::execution_space>::value;
  using execution_space =
      typename std::conditional<execution_space_is_defaulted,
                                DefaultExecutionSpace,
                                typename Base::execution_space>::type;

  using schedule_type =
      typename std::conditional<is_void<typename Base::schedule_type>::value,
                                Schedule<Static>,
                                typename Base::schedule_type>::type;

  using work_tag = typename Base::work_tag;

  using index_type =
      typename std::conditional<is_void<typename Base::index_type>::value,
                                IndexType<typename execution_space::size_type>,
                                typename Base::index_type>::type::type;
  // nasty hack to make index_type into an integral_type
  // instead of the wrapped IndexType<T> for backwards compatibility

  using iteration_pattern = typename std::conditional<
      is_void<typename Base::iteration_pattern>::value,
      void  // TODO set default iteration pattern
      ,
      typename Base::iteration_pattern>::type;

  using launch_bounds =
      typename std::conditional<is_void<typename Base::launch_bounds>::value,
                                LaunchBounds<>,
                                typename Base::launch_bounds>::type;

  using work_item_property = typename Base::work_item_property;

  using is_graph_kernel = typename Base::is_graph_kernel;

  using occupancy_control = typename Base::occupancy_control;

  using type =
      PolicyTraitsBase<execution_space, schedule_type, work_tag, index_type,
                       iteration_pattern, launch_bounds, work_item_property,
                       is_graph_kernel, occupancy_control>;
};

template <class AnalyzedPolicy>
struct PolicyDataStorage : AnalyzedPolicy,
                           NoUniqueAddressMemberEmulation<
                               typename AnalyzedPolicy::occupancy_control> {
  using occupancy_control_t = typename AnalyzedPolicy::occupancy_control;

  using occupancy_control_storage_base_t =
      NoUniqueAddressMemberEmulation<occupancy_control_t>;

  static constexpr bool experimental_contains_desired_occupancy =
      std::is_same<occupancy_control_t,
                   Kokkos::Experimental::DesiredOccupancy>::value;

  PolicyDataStorage() = default;

  // Converting constructors
  template <
      class Other,
      std::enable_if_t<
          experimental_contains_desired_occupancy &&
              PolicyDataStorage<Other>::experimental_contains_desired_occupancy,
          int> = 0>
  PolicyDataStorage(PolicyDataStorage<Other> const &other) {
    this->impl_set_desired_occupancy(other.impl_get_desired_occupancy());
  }

  template <class Other,
            std::enable_if_t<!experimental_contains_desired_occupancy ||
                                 !PolicyDataStorage<Other>::
                                     experimental_contains_desired_occupancy,
                             int> = 0>
  PolicyDataStorage(PolicyDataStorage<Other> const &) {}

  // Converting assignment operators
  template <
      class Other,
      std::enable_if_t<
          experimental_contains_desired_occupancy &&
              PolicyDataStorage<Other>::experimental_contains_desired_occupancy,
          int> = 0>
  PolicyDataStorage &operator=(PolicyDataStorage<Other> const &other) {
    this->impl_set_desired_occupancy(other.impl_get_desired_occupancy());
    return *this;
  }

  template <class Other,
            std::enable_if_t<!experimental_contains_desired_occupancy ||
                                 !PolicyDataStorage<Other>::
                                     experimental_contains_desired_occupancy,
                             int> = 0>
  PolicyDataStorage &operator=(PolicyDataStorage<Other> const &) {
    return *this;
  }

  // Access to desired occupancy (getter and setter)
  template <class Dummy = occupancy_control_t>
  std::enable_if_t<std::is_same<Dummy, occupancy_control_t>::value &&
                       experimental_contains_desired_occupancy,
                   Kokkos::Experimental::DesiredOccupancy>
  impl_get_desired_occupancy() const {
    return this
        ->occupancy_control_storage_base_t::no_unique_address_data_member();
  }

  template <class Dummy = occupancy_control_t>
  std::enable_if_t<std::is_same<Dummy, occupancy_control_t>::value &&
                   experimental_contains_desired_occupancy>
  impl_set_desired_occupancy(occupancy_control_t desired_occupancy) {
    this->occupancy_control_storage_base_t::no_unique_address_data_member() =
        desired_occupancy;
  }
};

template <typename... Traits>
struct PolicyTraits
    : PolicyDataStorage<
          typename AnalyzePolicy<PolicyTraitsBase<>, Traits...>::type> {
  using base_t = PolicyDataStorage<
      typename AnalyzePolicy<PolicyTraitsBase<>, Traits...>::type>;
  template <class... Args>
  PolicyTraits(PolicyTraits<Args...> const &p) : base_t(p) {}
  PolicyTraits() = default;
};

}  // namespace Impl
}  // namespace Kokkos

#endif  // KOKKOS_IMPL_ANALYZE_POLICY_HPP
