// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_KOKKOS_LAUNCHBOUNDSTRAIT_HPP
#define KOKKOS_KOKKOS_LAUNCHBOUNDSTRAIT_HPP

#include <Kokkos_Macros.hpp>
#include <Kokkos_Concepts.hpp>  // LaunchBounds
#include <traits/Kokkos_PolicyTraitAdaptor.hpp>
#include <traits/Kokkos_Traits_fwd.hpp>

namespace Kokkos {
namespace Impl {

//==============================================================================
// <editor-fold desc="trait specification"> {{{1

template <class LaunchBoundParam, class AnalyzeNextTrait>
struct LaunchBoundsMixin : AnalyzeNextTrait {
  using base_t = AnalyzeNextTrait;
  using base_t::base_t;

  static constexpr bool launch_bounds_is_defaulted = false;

  static_assert(base_t::launch_bounds_is_defaulted,
                "Kokkos Error: More than one launch_bounds given");

  using launch_bounds = LaunchBoundParam;
};

struct LaunchBoundsTrait : TraitSpecificationBase<LaunchBoundsTrait> {
  struct base_traits {
    static constexpr bool launch_bounds_is_defaulted = true;

    using launch_bounds = LaunchBounds<>;
    KOKKOS_IMPL_MSVC_NVCC_EBO_WORKAROUND
  };
  template <class LaunchBoundParam, class AnalyzeNextTrait>
  using mixin_matching_trait =
      LaunchBoundsMixin<LaunchBoundParam, AnalyzeNextTrait>;
};

// </editor-fold> end trait specification }}}1
//==============================================================================

//==============================================================================
// <editor-fold desc="PolicyTraitMatcher specialization"> {{{1

template <unsigned int maxT, unsigned int minB>
struct PolicyTraitMatcher<LaunchBoundsTrait, LaunchBounds<maxT, minB>>
    : std::true_type {};

// </editor-fold> end PolicyTraitMatcher specialization }}}1
//==============================================================================

}  // end namespace Impl
}  // end namespace Kokkos

#endif  // KOKKOS_KOKKOS_LAUNCHBOUNDSTRAIT_HPP
