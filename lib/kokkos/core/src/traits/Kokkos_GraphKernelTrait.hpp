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

struct GraphKernelTrait : TraitSpecificationBase<GraphKernelTrait> {
  struct base_traits {
    using is_graph_kernel = std::false_type;
    KOKKOS_IMPL_MSVC_NVCC_EBO_WORKAROUND
  };
  template <class, class AnalyzeNextTrait>
  struct mixin_matching_trait : AnalyzeNextTrait {
    using base_t = AnalyzeNextTrait;
    using base_t::base_t;
    using is_graph_kernel = std::true_type;
  };
  template <class T>
  using trait_matches_specification = std::is_same<T, IsGraphKernelTag>;
};

// </editor-fold> end trait specification }}}1
//==============================================================================

}  // end namespace Impl
}  // end namespace Kokkos

#endif  // KOKKOS_KOKKOS_GRAPHKERNELTRAIT_HPP
