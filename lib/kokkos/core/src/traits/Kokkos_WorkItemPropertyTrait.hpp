// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

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

template <class WorkItemProp, class AnalyzeNextTrait>
struct WorkItemPropMixin : AnalyzeNextTrait {
  using base_t = AnalyzeNextTrait;
  using base_t::base_t;
  using work_item_property = WorkItemProp;
};

struct WorkItemPropertyTrait : TraitSpecificationBase<WorkItemPropertyTrait> {
  struct base_traits {
    using work_item_property = Kokkos::Experimental::WorkItemProperty::None_t;
    KOKKOS_IMPL_MSVC_NVCC_EBO_WORKAROUND
  };
  template <class WorkItemProp, class AnalyzeNextTrait>
  using mixin_matching_trait =
      WorkItemPropMixin<WorkItemProp, AnalyzeNextTrait>;
  template <class T>
  using trait_matches_specification =
      Kokkos::Experimental::is_work_item_property<T>;
};

// </editor-fold> end trait specification }}}1
//==============================================================================

}  // end namespace Impl

namespace Experimental {

//==============================================================================
// <editor-fold desc="User interface"> {{{1

template <class Policy, unsigned long Property>
constexpr auto require(const Policy p,
                       WorkItemProperty::ImplWorkItemProperty<Property>) {
  static_assert(Kokkos::is_execution_policy<Policy>::value);
  using new_policy_t = Kokkos::Impl::WorkItemPropertyTrait::policy_with_trait<
      Policy, WorkItemProperty::ImplWorkItemProperty<Property>>;
  return new_policy_t{p};
}

// </editor-fold> end User interface }}}1
//==============================================================================

}  // namespace Experimental

}  // end namespace Kokkos

#endif  // KOKKOS_KOKKOS_WORKITEMPROPERTYTRAIT_HPP
