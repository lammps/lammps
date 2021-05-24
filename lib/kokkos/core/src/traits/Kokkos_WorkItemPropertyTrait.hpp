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

#ifndef KOKKOS_KOKKOS_WORKITEMPROPERTYTRAIT_HPP
#define KOKKOS_KOKKOS_WORKITEMPROPERTYTRAIT_HPP

#include <Kokkos_Macros.hpp>
#include <Kokkos_Concepts.hpp>  // WorkItemProperty
#include <traits/Kokkos_PolicyTraitAdaptor.hpp>
#include <traits/Kokkos_Traits_fwd.hpp>

namespace Kokkos {
namespace Impl {

//==============================================================================
// <editor-fold desc="trait specification"> {{{1

struct WorkItemPropertyTrait : TraitSpecificationBase<WorkItemPropertyTrait> {
  struct base_traits {
    using work_item_property = Kokkos::Experimental::WorkItemProperty::None_t;
  };
  template <class T>
  using trait_matches_specification =
      Kokkos::Experimental::is_work_item_property<T>;
};

// </editor-fold> end trait specification }}}1
//==============================================================================

//==============================================================================
// <editor-fold desc="AnalyzeExecPolicy specializations"> {{{1

template <class Property, class... Traits>
struct AnalyzeExecPolicy<
    std::enable_if_t<
        Kokkos::Experimental::is_work_item_property<Property>::value>,
    Property, Traits...> : AnalyzeExecPolicy<void, Traits...> {
  using base_t = AnalyzeExecPolicy<void, Traits...>;
  using base_t::base_t;
  static_assert(
      std::is_same<typename base_t::work_item_property,
                   Kokkos::Experimental::WorkItemProperty::None_t>::value,
      "Kokkos Error: More than one work item property given");
  using work_item_property = Property;
};

// </editor-fold> end AnalyzeExecPolicy specializations }}}1
//==============================================================================

}  // end namespace Impl

namespace Experimental {

//==============================================================================
// <editor-fold desc="User interface"> {{{1

template <class Policy, unsigned long Property>
constexpr auto require(const Policy p,
                       WorkItemProperty::ImplWorkItemProperty<Property>) {
  static_assert(Kokkos::is_execution_policy<Policy>::value, "");
  using new_policy_t = Kokkos::Impl::WorkItemPropertyTrait::policy_with_trait<
      Policy, WorkItemProperty::ImplWorkItemProperty<Property>>;
  return new_policy_t{p};
}

// </editor-fold> end User interface }}}1
//==============================================================================

}  // namespace Experimental

}  // end namespace Kokkos

#endif  // KOKKOS_KOKKOS_WORKITEMPROPERTYTRAIT_HPP
