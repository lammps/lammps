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

#ifndef KOKKOS_KOKKOS_WORKTAGTRAIT_HPP
#define KOKKOS_KOKKOS_WORKTAGTRAIT_HPP

#include <Kokkos_Macros.hpp>
#include <Kokkos_Concepts.hpp>  // is_execution_space
#include <traits/Kokkos_PolicyTraitAdaptor.hpp>
#include <traits/Kokkos_Traits_fwd.hpp>
#include <impl/Kokkos_Utilities.hpp>  // type_list_any, type_list_remove_first

namespace Kokkos {
namespace Impl {

//==============================================================================
// <editor-fold desc="trait specification"> {{{1

template <class T>
struct show_extra_work_tag_erroneously_given_to_execution_policy;
template <>
struct show_extra_work_tag_erroneously_given_to_execution_policy<void> {};

using _exec_policy_traits_without_work_tag = typename type_list_remove_first<
    WorkTagTrait, execution_policy_trait_specifications>::type;

template <class Trait>
struct _trait_matches_spec_predicate {
  template <class TraitSpec>
  struct apply {
    using type = typename PolicyTraitMatcher<TraitSpec, Trait>::type;
    static constexpr bool value = type::value;
  };
};

struct WorkTagTrait : TraitSpecificationBase<WorkTagTrait> {
  struct base_traits {
    using work_tag = void;
    KOKKOS_IMPL_MSVC_NVCC_EBO_WORKAROUND
  };
  template <class WorkTag, class AnalyzeNextTrait>
  struct mixin_matching_trait : AnalyzeNextTrait {
    using base_t = AnalyzeNextTrait;
    using base_t::base_t;
    using work_tag = WorkTag;
    static constexpr auto show_work_tag_error_in_compilation_message =
        show_extra_work_tag_erroneously_given_to_execution_policy<
            typename base_t::work_tag>{};
    static_assert(
        std::is_void<typename base_t::work_tag>::value,
        "Kokkos Error: More than one work tag given. Search compiler output "
        "for 'show_extra_work_tag' to see the type of the errant tag.");
  };
  // Since we don't have subsumption in pre-C++20, we need to have the work tag
  // "trait" handling code ensure that none of the other conditions are met.
  // * Compile time cost complexity note: at first glance it looks like this
  //   "rechecks" all of the other trait specs when used in the context of the
  //   full list of execution policy traits, but actually since we've already
  //   checked all of them to get to the end of the list, the compiler will
  //   have already generated those definitions, so there should be little extra
  //   cost to this. However, in the scenario where we use work tag in isolation
  //   (like if we were to add a `require()`-like thing that changes the work
  //   tag of an existing execution policy instance), we need to check all of
  //   the other traits to make sure that we're not replacing something else,
  //   given that the concept of a work tag is basically unconstrained and could
  //   be anything.  This should still be as efficient at compile time as the
  //   old code that just did a big long series of nested std::conditionals, but
  //   we should benchmark this assumption if it becomes a problem.
  template <class T>
  using trait_matches_specification = std::bool_constant<
      std::is_empty<T>::value &&
      !type_list_any<_trait_matches_spec_predicate<T>::template apply,
                     _exec_policy_traits_without_work_tag>::value>;
};

// </editor-fold> end trait specification }}}1
//==============================================================================

}  // end namespace Impl
}  // end namespace Kokkos

#endif  // KOKKOS_KOKKOS_WORKTAGTRAIT_HPP
