// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_KOKKOS_GRAPHKERNELTRAIT_HPP
#define KOKKOS_KOKKOS_GRAPHKERNELTRAIT_HPP

#include <Kokkos_Macros.hpp>
#include <traits/Kokkos_PolicyTraitAdaptor.hpp>
#include <impl/Kokkos_GraphImpl_fwd.hpp>  // IsGraphKernelTag
#include <traits/Kokkos_Traits_fwd.hpp>
#include <impl/Kokkos_Utilities.hpp>

namespace Kokkos {
namespace Impl {

//==============================================================================
// <editor-fold desc="trait specification"> {{{1

template <class, class AnalyzeNextTrait>
struct GraphMixin : AnalyzeNextTrait {
  using base_t = AnalyzeNextTrait;
  using base_t::base_t;
  using is_graph_kernel = std::true_type;
};

struct GraphKernelTrait : TraitSpecificationBase<GraphKernelTrait> {
  struct base_traits {
    using is_graph_kernel = std::false_type;
    KOKKOS_IMPL_MSVC_NVCC_EBO_WORKAROUND
  };
  template <class NotUsed, class AnalyzeNextTrait>
  using mixin_matching_trait = GraphMixin<NotUsed, AnalyzeNextTrait>;
  template <class T>
  using trait_matches_specification = std::is_same<T, IsGraphKernelTag>;
};

// </editor-fold> end trait specification }}}1
//==============================================================================

}  // end namespace Impl
}  // end namespace Kokkos

#endif  // KOKKOS_KOKKOS_GRAPHKERNELTRAIT_HPP
