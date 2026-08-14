// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_NESTED_SORT_PUBLIC_API_HPP_
#define KOKKOS_NESTED_SORT_PUBLIC_API_HPP_

#include "impl/Kokkos_NestedSortImpl.hpp"
#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif
#include <std_algorithms/impl/Kokkos_HelperPredicates.hpp>

namespace Kokkos {
namespace Experimental {

template <class TeamMember, class ViewType>
KOKKOS_INLINE_FUNCTION void sort_team(const TeamMember& t,
                                      const ViewType& view) {
  Impl::sort_nested_impl(t, view, nullptr,
                         Experimental::Impl::StdAlgoLessThanBinaryPredicate<
                             typename ViewType::non_const_value_type>(),
                         Impl::NestedRange<true>());
}

template <class TeamMember, class ViewType, class Comparator>
KOKKOS_INLINE_FUNCTION void sort_team(const TeamMember& t, const ViewType& view,
                                      const Comparator& comp) {
  Impl::sort_nested_impl(t, view, nullptr, comp, Impl::NestedRange<true>());
}

template <class TeamMember, class KeyViewType, class ValueViewType>
KOKKOS_INLINE_FUNCTION void sort_by_key_team(const TeamMember& t,
                                             const KeyViewType& keyView,
                                             const ValueViewType& valueView) {
  Impl::sort_nested_impl(t, keyView, valueView,
                         Experimental::Impl::StdAlgoLessThanBinaryPredicate<
                             typename KeyViewType::non_const_value_type>(),
                         Impl::NestedRange<true>());
}

template <class TeamMember, class KeyViewType, class ValueViewType,
          class Comparator>
KOKKOS_INLINE_FUNCTION void sort_by_key_team(const TeamMember& t,
                                             const KeyViewType& keyView,
                                             const ValueViewType& valueView,
                                             const Comparator& comp) {
  Impl::sort_nested_impl(t, keyView, valueView, comp,
                         Impl::NestedRange<true>());
}

template <class TeamMember, class ViewType>
KOKKOS_INLINE_FUNCTION void sort_thread(const TeamMember& t,
                                        const ViewType& view) {
  Impl::sort_nested_impl(t, view, nullptr,
                         Experimental::Impl::StdAlgoLessThanBinaryPredicate<
                             typename ViewType::non_const_value_type>(),
                         Impl::NestedRange<false>());
}

template <class TeamMember, class ViewType, class Comparator>
KOKKOS_INLINE_FUNCTION void sort_thread(const TeamMember& t,
                                        const ViewType& view,
                                        const Comparator& comp) {
  Impl::sort_nested_impl(t, view, nullptr, comp, Impl::NestedRange<false>());
}

template <class TeamMember, class KeyViewType, class ValueViewType>
KOKKOS_INLINE_FUNCTION void sort_by_key_thread(const TeamMember& t,
                                               const KeyViewType& keyView,
                                               const ValueViewType& valueView) {
  Impl::sort_nested_impl(t, keyView, valueView,
                         Experimental::Impl::StdAlgoLessThanBinaryPredicate<
                             typename KeyViewType::non_const_value_type>(),
                         Impl::NestedRange<false>());
}

template <class TeamMember, class KeyViewType, class ValueViewType,
          class Comparator>
KOKKOS_INLINE_FUNCTION void sort_by_key_thread(const TeamMember& t,
                                               const KeyViewType& keyView,
                                               const ValueViewType& valueView,
                                               const Comparator& comp) {
  Impl::sort_nested_impl(t, keyView, valueView, comp,
                         Impl::NestedRange<false>());
}

}  // namespace Experimental
}  // namespace Kokkos
#endif
