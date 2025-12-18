// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_STATIC_BATCH_SIZE_TRAIT_HPP
#define KOKKOS_STATIC_BATCH_SIZE_TRAIT_HPP

#include <Kokkos_Macros.hpp>
#include <traits/Kokkos_PolicyTraitAdaptor.hpp>
#include <traits/Kokkos_Traits_fwd.hpp>

namespace Kokkos::Experimental {

template <unsigned int BatchSize = 1>
struct StaticBatchSize {
  using static_batch_size = StaticBatchSize;
  using type              = StaticBatchSize<BatchSize>;
  constexpr static unsigned int batch_size{BatchSize};

  static_assert(
      BatchSize > 0,
      "Kokkos Error: StaticBatchSize factor must be greater than zero");
};

}  // end namespace Kokkos::Experimental

namespace Kokkos::Impl {

//==============================================================================
// <editor-fold desc="trait specification"> {{{1

template <class StaticBatchSizeParam, class AnalyzeNextTrait>
struct StaticBatchSizeMixin : AnalyzeNextTrait {
  using base_t = AnalyzeNextTrait;
  using base_t::base_t;

  static constexpr bool batch_size_is_defaulted = false;

  static_assert(
      base_t::batch_size_is_defaulted,
      "Kokkos Error: More than one StaticBatchSizeTrait specified is given.");

  using static_batch_size = StaticBatchSizeParam;
};

struct StaticBatchSizeTrait : TraitSpecificationBase<StaticBatchSizeTrait> {
  struct base_traits {
    static constexpr bool batch_size_is_defaulted = true;

    using static_batch_size = Kokkos::Experimental::StaticBatchSize<>;
    KOKKOS_IMPL_MSVC_NVCC_EBO_WORKAROUND
  };
  template <class StaticBatchSizeParam, class AnalyzeNextTrait>
  using mixin_matching_trait =
      StaticBatchSizeMixin<StaticBatchSizeParam, AnalyzeNextTrait>;
};
}  // end namespace Kokkos::Impl

// </editor-fold> end trait specification }}}1
//==============================================================================

//==============================================================================
// <editor-fold desc="PolicyTraitMatcher specialization"> {{{1

namespace Kokkos::Impl {

template <unsigned int BatchSize>
struct PolicyTraitMatcher<StaticBatchSizeTrait,
                          Kokkos::Experimental::StaticBatchSize<BatchSize>>
    : std::true_type {};

// </editor-fold> end PolicyTraitMatcher specialization }}}1
//==============================================================================

}  // end namespace Kokkos::Impl

#endif  // KOKKOS_STATIC_BATCH_SIZE_TRAIT_HPP
