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

#ifndef KOKKOS_KOKKOS_INDEXTYPETRAIT_HPP
#define KOKKOS_KOKKOS_INDEXTYPETRAIT_HPP

#include <Kokkos_Macros.hpp>
#include <Kokkos_Concepts.hpp>  // IndexType, is_index_type
#include <traits/Kokkos_PolicyTraitAdaptor.hpp>
#include <traits/Kokkos_Traits_fwd.hpp>

namespace Kokkos {
namespace Impl {

//==============================================================================
// <editor-fold desc="trait specification"> {{{1

struct IndexTypeTrait : TraitSpecificationBase<IndexTypeTrait> {
  struct base_traits {
    static constexpr bool index_type_is_defaulted = true;
    using index_type = dependent_policy_trait_default;
  };
  template <class T>
  using trait_matches_specification =
      std::integral_constant<bool, std::is_integral<T>::value ||
                                       is_index_type<T>::value>;
};

// </editor-fold> end trait specification }}}1
//==============================================================================

//==============================================================================
// <editor-fold desc="AnalyzeExecPolicy specializations"> {{{1

// Index type given as IndexType template
template <class IntegralIndexType, class... Traits>
struct AnalyzeExecPolicy<void, Kokkos::IndexType<IntegralIndexType>, Traits...>
    : AnalyzeExecPolicy<void, Traits...> {
  using base_t = AnalyzeExecPolicy<void, Traits...>;
  using base_t::base_t;
  static_assert(base_t::index_type_is_defaulted,
                "Kokkos Error: More than one index type given");
  static constexpr bool index_type_is_defaulted = false;
  using index_type = Kokkos::IndexType<IntegralIndexType>;
};

// IndexType given as an integral type directly
template <class IntegralIndexType, class... Traits>
struct AnalyzeExecPolicy<
    std::enable_if_t<std::is_integral<IntegralIndexType>::value>,
    IntegralIndexType, Traits...> : AnalyzeExecPolicy<void, Traits...> {
  using base_t = AnalyzeExecPolicy<void, Traits...>;
  using base_t::base_t;
  static_assert(base_t::index_type_is_defaulted,
                "Kokkos Error: More than one index type given");
  static constexpr bool index_type_is_defaulted = false;
  using index_type = Kokkos::IndexType<IntegralIndexType>;
};

// </editor-fold> end AnalyzeExecPolicy specializations }}}1
//==============================================================================

}  // end namespace Impl
}  // end namespace Kokkos

#endif  // KOKKOS_KOKKOS_INDEXTYPETRAIT_HPP
