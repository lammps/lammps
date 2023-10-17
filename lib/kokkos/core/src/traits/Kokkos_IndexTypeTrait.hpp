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

#ifndef KOKKOS_KOKKOS_INDEXTYPETRAIT_HPP
#define KOKKOS_KOKKOS_INDEXTYPETRAIT_HPP

#include <Kokkos_Macros.hpp>
#include <Kokkos_Concepts.hpp>  // IndexType
#include <traits/Kokkos_PolicyTraitAdaptor.hpp>
#include <traits/Kokkos_Traits_fwd.hpp>

namespace Kokkos {
namespace Impl {

template <class Trait, class AnalyzeNextTrait>
struct IndexTypePolicyMixin;

//==============================================================================
// <editor-fold desc="trait specification"> {{{1

template <class T>
struct show_extra_index_type_erroneously_given_to_execution_policy;
template <>
struct show_extra_index_type_erroneously_given_to_execution_policy<void> {};
struct IndexTypeTrait : TraitSpecificationBase<IndexTypeTrait> {
  struct base_traits {
    static constexpr bool index_type_is_defaulted = true;
    using index_type = dependent_policy_trait_default;
    KOKKOS_IMPL_MSVC_NVCC_EBO_WORKAROUND
  };
  template <class IdxType, class AnalyzeNextTrait>
  using mixin_matching_trait = IndexTypePolicyMixin<IdxType, AnalyzeNextTrait>;
};

// </editor-fold> end trait specification }}}1
//==============================================================================

//==============================================================================
// <editor-fold desc="IndexTypePolicyMixin specializations"> {{{1

// Index type given as IndexType template
template <class IntegralIndexType, class AnalyzeNextTrait>
struct IndexTypePolicyMixin<Kokkos::IndexType<IntegralIndexType>,
                            AnalyzeNextTrait> : AnalyzeNextTrait {
  using base_t = AnalyzeNextTrait;
  using base_t::base_t;
  static constexpr auto show_index_type_error_in_compilation_message =
      show_extra_index_type_erroneously_given_to_execution_policy<
          std::conditional_t<base_t::index_type_is_defaulted, void,
                             typename base_t::schedule_type>>{};
  static_assert(base_t::index_type_is_defaulted,
                "Kokkos Error: More than one index type given. Search "
                "compiler output for 'show_extra_index_type' to see the "
                "type of the errant tag.");
  static constexpr bool index_type_is_defaulted = false;
  using index_type = Kokkos::IndexType<IntegralIndexType>;
};

// IndexType given as an integral type directly (the matcher already checks
// this, so we don't have specialize to re-check it here)
template <class IntegralIndexType, class AnalyzeNextTrait>
struct IndexTypePolicyMixin : AnalyzeNextTrait {
  using base_t = AnalyzeNextTrait;
  using base_t::base_t;
  static constexpr auto show_index_type_error_in_compilation_message =
      show_extra_index_type_erroneously_given_to_execution_policy<
          std::conditional_t<base_t::index_type_is_defaulted, void,
                             typename base_t::schedule_type>>{};
  static_assert(base_t::index_type_is_defaulted,
                "Kokkos Error: More than one index type given. Search "
                "compiler output for 'show_extra_index_type' to see the "
                "type of the errant tag.");
  static_assert(std::is_integral<IntegralIndexType>::value, "");
  static constexpr bool index_type_is_defaulted = false;
  using index_type = Kokkos::IndexType<IntegralIndexType>;
};

// </editor-fold> end AnalyzeExecPolicy specializations }}}1
//==============================================================================

//==============================================================================
// <editor-fold desc="PolicyTraitMatcher specialization"> {{{1

template <class IntegralIndexType>
struct PolicyTraitMatcher<IndexTypeTrait, IndexType<IntegralIndexType>>
    : std::true_type {};

template <class IntegralIndexType>
struct PolicyTraitMatcher<
    IndexTypeTrait, IntegralIndexType,
    std::enable_if_t<std::is_integral<IntegralIndexType>::value>>
    : std::true_type {};

// </editor-fold> end PolicyTraitMatcher specialization"> }}}1
//==============================================================================

}  // end namespace Impl
}  // end namespace Kokkos

#endif  // KOKKOS_KOKKOS_INDEXTYPETRAIT_HPP
