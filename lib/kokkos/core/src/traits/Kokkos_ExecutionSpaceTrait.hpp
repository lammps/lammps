// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_KOKKOS_EXECUTIONSPACETRAIT_HPP
#define KOKKOS_KOKKOS_EXECUTIONSPACETRAIT_HPP

#include <Kokkos_Macros.hpp>
#include <Kokkos_Concepts.hpp>  // is_execution_space
#include <traits/Kokkos_PolicyTraitAdaptor.hpp>
#include <traits/Kokkos_Traits_fwd.hpp>

namespace Kokkos {
namespace Impl {

//==============================================================================
// <editor-fold desc="trait specification"> {{{1
template <class T>
struct show_extra_execution_space_erroneously_given_to_execution_policy;
template <>
struct show_extra_execution_space_erroneously_given_to_execution_policy<void> {
};

template <class ExecSpace, class AnalyzeNextTrait>
struct ExecutionSpaceMixin : AnalyzeNextTrait {
  using base_t = AnalyzeNextTrait;
  using base_t::base_t;

  static constexpr auto show_execution_space_error_in_compilation_message =
      show_extra_execution_space_erroneously_given_to_execution_policy<
          std::conditional_t<base_t::execution_space_is_defaulted, void,
                             typename base_t::execution_space>>{};
  static_assert(base_t::execution_space_is_defaulted,
                "Kokkos Error: More than one execution space given. Search "
                "compiler output for 'show_extra_execution_space' to see the "
                "type of the errant tag.");

  static constexpr auto execution_space_is_defaulted = false;

  using execution_space = ExecSpace;
};

struct ExecutionSpaceTrait : TraitSpecificationBase<ExecutionSpaceTrait> {
  struct base_traits {
    static constexpr auto execution_space_is_defaulted = true;

    using execution_space = Kokkos::DefaultExecutionSpace;
    KOKKOS_IMPL_MSVC_NVCC_EBO_WORKAROUND
  };
  template <class T>
  using trait_matches_specification = Kokkos::is_execution_space<T>;
  template <class ExecSpace, class AnalyzeNextTrait>
  using mixin_matching_trait = ExecutionSpaceMixin<ExecSpace, AnalyzeNextTrait>;
};

// </editor-fold> end trait specification }}}1
//==============================================================================

}  // end namespace Impl
}  // end namespace Kokkos

#endif  // KOKKOS_KOKKOS_EXECUTIONSPACETRAIT_HPP
