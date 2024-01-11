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
struct ExecutionSpaceTrait : TraitSpecificationBase<ExecutionSpaceTrait> {
  struct base_traits {
    static constexpr auto execution_space_is_defaulted = true;

    using execution_space = Kokkos::DefaultExecutionSpace;
    KOKKOS_IMPL_MSVC_NVCC_EBO_WORKAROUND
  };
  template <class T>
  using trait_matches_specification = Kokkos::is_execution_space<T>;
  template <class ExecSpace, class AnalyzeNextTrait>
  struct mixin_matching_trait : AnalyzeNextTrait {
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
};

// </editor-fold> end trait specification }}}1
//==============================================================================

}  // end namespace Impl
}  // end namespace Kokkos

#endif  // KOKKOS_KOKKOS_EXECUTIONSPACETRAIT_HPP
