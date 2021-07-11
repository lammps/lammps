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

#ifndef KOKKOS_KOKKOS_WORKTAGTRAIT_HPP
#define KOKKOS_KOKKOS_WORKTAGTRAIT_HPP

#include <Kokkos_Macros.hpp>
#include <Kokkos_Concepts.hpp>  // is_execution_space
#include <traits/Kokkos_PolicyTraitAdaptor.hpp>
#include <traits/Kokkos_Traits_fwd.hpp>

namespace Kokkos {
namespace Impl {

//==============================================================================
// <editor-fold desc="trait specification"> {{{1

struct WorkTagTrait : TraitSpecificationBase<WorkTagTrait> {
  struct base_traits {
    using work_tag = void;
  };
};

// </editor-fold> end trait specification }}}1
//==============================================================================

//==============================================================================
// <editor-fold desc="AnalyzeExecPolicy specializations"> {{{1

// Since we don't have subsumption in pre-C++20, we need to have the work tag
// "trait" handling code be unspecialized, so we handle it instead in a class
// with a different name.
template <class... Traits>
struct AnalyzeExecPolicyHandleWorkTag : AnalyzeExecPolicy<void, Traits...> {
  using base_t = AnalyzeExecPolicy<void, Traits...>;
  using base_t::base_t;
};

template <class WorkTag, class... Traits>
struct AnalyzeExecPolicyHandleWorkTag<WorkTag, Traits...>
    : AnalyzeExecPolicy<void, Traits...> {
  using base_t = AnalyzeExecPolicy<void, Traits...>;
  using base_t::base_t;
  static_assert(std::is_void<typename base_t::work_tag>::value,
                "Kokkos Error: More than one work tag given");
  using work_tag = WorkTag;
};

// This only works if this is not a partial specialization, so we have to
// do the partial specialization elsewhere
template <class Enable, class... Traits>
struct AnalyzeExecPolicy : AnalyzeExecPolicyHandleWorkTag<Traits...> {
  using base_t = AnalyzeExecPolicyHandleWorkTag<Traits...>;
  using base_t::base_t;
};

// </editor-fold> end AnalyzeExecPolicy specializations }}}1
//==============================================================================

//==============================================================================
// <editor-fold desc="PolicyTraitMatcher specializations"> {{{1

// In order to match the work tag trait the work tag "matcher" needs to be
// unspecialized and the logic needs to be handled in a differently-named class,
// just like above.
template <class TraitSpec, class Trait>
struct PolicyTraitMatcherHandleWorkTag : std::false_type {};

template <class Trait>
struct PolicyTraitMatcherHandleWorkTag<WorkTagTrait, Trait>
    : std::integral_constant<bool, !std::is_void<Trait>::value> {};

template <class TraitSpec, class Trait, class Enable>
struct PolicyTraitMatcher /* unspecialized! */
    : PolicyTraitMatcherHandleWorkTag<TraitSpec, Trait> {};

// </editor-fold> end PolicyTraitMatcher specializations }}}1
//==============================================================================

}  // end namespace Impl
}  // end namespace Kokkos

#endif  // KOKKOS_KOKKOS_WORKTAGTRAIT_HPP
