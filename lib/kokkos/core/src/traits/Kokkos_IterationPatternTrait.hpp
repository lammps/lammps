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

#ifndef KOKKOS_KOKKOS_ITERATIONPATTERNTRAIT_HPP
#define KOKKOS_KOKKOS_ITERATIONPATTERNTRAIT_HPP

#include <Kokkos_Concepts.hpp>  // is_iteration_pattern
#include <type_traits>          // is_void

namespace Kokkos {
namespace Impl {

//==============================================================================
// <editor-fold desc="trait specification"> {{{1

struct IterationPatternTrait : TraitSpecificationBase<IterationPatternTrait> {
  struct base_traits {
    using iteration_pattern = void;  // TODO set default iteration pattern
  };
  template <class T>
  using trait_matches_specification = is_iteration_pattern<T>;
};

// </editor-fold> end trait specification }}}1
//==============================================================================

//==============================================================================
// <editor-fold desc="AnalyzeExecPolicy specializations"> {{{1

template <class IterationPattern, class... Traits>
struct AnalyzeExecPolicy<
    std::enable_if_t<is_iteration_pattern<IterationPattern>::value>,
    IterationPattern, Traits...> : AnalyzeExecPolicy<void, Traits...> {
  using base_t = AnalyzeExecPolicy<void, Traits...>;
  using base_t::base_t;
  static_assert(std::is_void<typename base_t::iteration_pattern>::value,
                "Kokkos Error: More than one iteration pattern given");
  using iteration_pattern = IterationPattern;
};

// </editor-fold> end AnalyzeExecPolicy specializations }}}1
//==============================================================================

}  // end namespace Impl
}  // end namespace Kokkos

#endif  // KOKKOS_KOKKOS_ITERATIONPATTERNTRAIT_HPP
