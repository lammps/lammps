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

#ifndef KOKKOS_KOKKOS_SCHEDULETRAIT_HPP
#define KOKKOS_KOKKOS_SCHEDULETRAIT_HPP

#include <Kokkos_Macros.hpp>
#include <Kokkos_Concepts.hpp>  // is_schedule_type, Schedule
#include <traits/Kokkos_PolicyTraitAdaptor.hpp>
#include <traits/Kokkos_Traits_fwd.hpp>

namespace Kokkos {

namespace Impl {

//==============================================================================
// <editor-fold desc="trait specification"> {{{1

template <class T>
struct show_extra_schedule_type_erroneously_given_to_execution_policy;
template <>
struct show_extra_schedule_type_erroneously_given_to_execution_policy<void> {};
struct ScheduleTrait : TraitSpecificationBase<ScheduleTrait> {
  struct base_traits {
    static constexpr auto schedule_type_is_defaulted = true;

    using schedule_type = Schedule<Static>;
    KOKKOS_IMPL_MSVC_NVCC_EBO_WORKAROUND
  };
  template <class Sched, class AnalyzeNextTrait>
  struct mixin_matching_trait : AnalyzeNextTrait {
    using base_t = AnalyzeNextTrait;
    using base_t::base_t;
    using schedule_type = Sched;
    static constexpr auto show_schedule_type_error_in_compilation_message =
        show_extra_schedule_type_erroneously_given_to_execution_policy<
            std::conditional_t<base_t::schedule_type_is_defaulted, void,
                               typename base_t::schedule_type>>{};
    static_assert(base_t::schedule_type_is_defaulted,
                  "Kokkos Error: More than one schedule type given. Search "
                  "compiler output for 'show_extra_schedule_type' to see the "
                  "type of the errant tag.");
    static constexpr bool schedule_type_is_defaulted = false;
  };
};

// </editor-fold> end trait specification }}}1
//==============================================================================

//==============================================================================
// <editor-fold desc="PolicyTraitMatcher specialization"> {{{1

template <class Sched>
struct PolicyTraitMatcher<ScheduleTrait, Schedule<Sched>> : std::true_type {};

// </editor-fold> end PolicyTraitMatcher specialization }}}1
//==============================================================================

}  // end namespace Impl

namespace Experimental {

//==============================================================================
// <editor-fold desc="User interface"> {{{1

template <class Policy, class ScheduleType>
constexpr auto require(Policy const& p, Kokkos::Schedule<ScheduleType>) {
  static_assert(Kokkos::is_execution_policy<Policy>::value);
  using new_policy_t = Kokkos::Impl::ScheduleTrait::policy_with_trait<
      Policy, Kokkos::Schedule<ScheduleType>>;
  return new_policy_t{p};
}

// </editor-fold> end User interface }}}1
//==============================================================================

}  // end namespace Experimental

}  // end namespace Kokkos

#endif  // KOKKOS_KOKKOS_SCHEDULETRAIT_HPP
