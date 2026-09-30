// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_KOKKOS_ITERATIONPATTERNTRAIT_HPP
#define KOKKOS_KOKKOS_ITERATIONPATTERNTRAIT_HPP

#include <Kokkos_Concepts.hpp>                   // is_iteration_pattern
#include <traits/Kokkos_PolicyTraitAdaptor.hpp>  // TraitSpecificationBase
#include <Kokkos_Rank.hpp>                       // Rank
#include <Kokkos_Layout.hpp>                     // Iterate
#include <type_traits>                           // is_void

namespace Kokkos {
namespace Impl {

//==============================================================================
// <editor-fold desc="trait specification"> {{{1

template <class T>
struct show_extra_iteration_pattern_erroneously_given_to_execution_policy;
template <>
struct show_extra_iteration_pattern_erroneously_given_to_execution_policy<
    void> {};

template <class IterPattern, class AnalyzeNextTrait>
struct IterPatternMixin : AnalyzeNextTrait {
  using base_t = AnalyzeNextTrait;
  using base_t::base_t;
  static constexpr auto show_iteration_pattern_error_in_compilation_message =
      show_extra_iteration_pattern_erroneously_given_to_execution_policy<
          typename base_t::iteration_pattern>{};
  static_assert(std::is_void_v<typename base_t::iteration_pattern>,
                "Kokkos Error: More than one index type given. Search "
                "compiler output for 'show_extra_iteration_pattern' to see the "
                "type of the errant tag.");
  using iteration_pattern = IterPattern;
};

struct IterationPatternTrait : TraitSpecificationBase<IterationPatternTrait> {
  struct base_traits {
    using iteration_pattern = void;  // TODO set default iteration pattern
    KOKKOS_IMPL_MSVC_NVCC_EBO_WORKAROUND
  };
  template <class IterPattern, class AnalyzeNextTrait>
  using mixin_matching_trait = IterPatternMixin<IterPattern, AnalyzeNextTrait>;
};

// </editor-fold> end trait specification }}}1
//==============================================================================

//==============================================================================
// <editor-fold desc="PolicyTraitMatcher specialization"> {{{1

template <unsigned N, Iterate OuterDir, Iterate InnerDir>
struct PolicyTraitMatcher<IterationPatternTrait, Rank<N, OuterDir, InnerDir>>
    : std::true_type {};

// </editor-fold> end  }}}1
//==============================================================================

}  // end namespace Impl
}  // end namespace Kokkos

#endif  // KOKKOS_KOKKOS_ITERATIONPATTERNTRAIT_HPP
