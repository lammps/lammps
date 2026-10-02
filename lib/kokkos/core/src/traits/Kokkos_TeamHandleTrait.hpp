// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_TEAM_HANDLE_TRAIT_HPP
#define KOKKOS_TEAM_HANDLE_TRAIT_HPP

#include <Kokkos_Macros.hpp>
#include <traits/Kokkos_PolicyTraitAdaptor.hpp>
#include <traits/Kokkos_Traits_fwd.hpp>

namespace Kokkos::Impl {

//==============================================================================
// <editor-fold desc="trait specification"> {{{1

template <class TeamHandle, class AnalyzeNextTrait>
struct TeamHandleMixin : AnalyzeNextTrait {
  using base_t = AnalyzeNextTrait;
  using base_t::base_t;

  static_assert(
      std::is_void_v<typename base_t::team_handle>,
      "Kokkos Error: More than one TeamHandleTrait specified is given.");
  static constexpr bool team_handle_is_defaulted = false;
  using team_handle                              = TeamHandle;
};

struct TeamHandleTrait : TraitSpecificationBase<TeamHandleTrait> {
  struct base_traits {
    static constexpr bool team_handle_is_defaulted = true;

    using team_handle = void;
    KOKKOS_IMPL_MSVC_NVCC_EBO_WORKAROUND
  };
  template <class TeamHandle, class AnalyzeNextTrait>
  using mixin_matching_trait = TeamHandleMixin<TeamHandle, AnalyzeNextTrait>;
  template <class T>
  using trait_matches_specification = is_team_handle<T>;
};

// </editor-fold> end trait specification }}}1
//==============================================================================

}  // end namespace Kokkos::Impl

#endif  // KOKKOS_TEAM_HANDLE_TRAIT_HPP
