// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_OPENACC_SCHEDULE_TYPE_HPP
#define KOKKOS_OPENACC_SCHEDULE_TYPE_HPP

#include <Kokkos_Concepts.hpp>
#include <type_traits>

namespace Kokkos::Experimental::Impl {

template <class Policy, class ScheduleType = typename Policy::schedule_type>
struct OpenACCSchedule {
  static_assert(is_execution_policy_v<Policy>);
  static_assert(std::is_void_v<ScheduleType> ||
                std::is_same_v<ScheduleType, Schedule<Static>> ||
                std::is_same_v<ScheduleType, Schedule<Dynamic>>);
  using type =
      std::conditional_t<std::is_same_v<ScheduleType, Schedule<Static>>,
                         Schedule<Static>, Schedule<Dynamic>>;
};

template <class Policy>
using OpenACCScheduleType = typename OpenACCSchedule<Policy>::type;

}  // namespace Kokkos::Experimental::Impl

#endif
