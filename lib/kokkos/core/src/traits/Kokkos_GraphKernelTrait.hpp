/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

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
