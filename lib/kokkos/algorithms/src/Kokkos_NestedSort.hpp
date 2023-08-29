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

#ifndef KOKKOS_NESTEDSORT_HPP_
#define KOKKOS_NESTEDSORT_HPP_

#include <Kokkos_Core.hpp>
#include <std_algorithms/impl/Kokkos_HelperPredicates.hpp>
#include <std_algorithms/Kokkos_Swap.hpp>

namespace Kokkos {
namespace Experimental {
namespace Impl {

// true for TeamVectorRange, false for ThreadVectorRange
template <bool teamLevel>
struct NestedRange {};

// Specialization for team-level
template <>
struct NestedRange<true> {
  template <typename TeamMember, typename SizeType>
  KOKKOS_FUNCTION static auto create(const TeamMember& t, SizeType len) {
    return Kokkos::TeamVectorRange(t, len);
  }
  template <typename TeamMember>
  KOKKOS_FUNCTION static void barrier(const TeamMember& t) {
    t.team_barrier();
  }
};

// Specialization for thread-level
template <>
struct NestedRange<false> {
  template <typename TeamMember, typename SizeType>
  KOKKOS_FUNCTION static auto create(const TeamMember& t, SizeType len) {
    return Kokkos::ThreadVectorRange(t, len);
  }
  // Barrier is no-op, as vector lanes of a thread are implicitly synchronized
  // after parallel region
  template <typename TeamMember>
  KOKKOS_FUNCTION static void barrier(const TeamMember&) {}
};

// When just doing sort (not sort_by_key), use nullptr_t for ValueViewType.
// This only takes the NestedRange instance for template arg deduction.
template <class TeamMember, class KeyViewType, class ValueViewType,
          class Comparator, bool useTeamLevel>
KOKKOS_INLINE_FUNCTION void sort_nested_impl(
    const TeamMember& t, const KeyViewType& keyView,
    [[maybe_unused]] const ValueViewType& valueView, const Comparator& comp,
    const NestedRange<useTeamLevel>) {
  using SizeType  = typename KeyViewType::size_type;
  using KeyType   = typename KeyViewType::non_const_value_type;
  using Range     = NestedRange<useTeamLevel>;
  SizeType n      = keyView.extent(0);
  SizeType npot   = 1;
  SizeType levels = 0;
  // FIXME: ceiling power-of-two is a common thing to need - make it a utility
  while (npot < n) {
    levels++;
    npot <<= 1;
  }
  for (SizeType i = 0; i < levels; i++) {
    for (SizeType j = 0; j <= i; j++) {
      // n/2 pairs of items are compared in parallel
      Kokkos::parallel_for(Range::create(t, npot / 2), [=](const SizeType k) {
        // How big are the brown/pink boxes?
        // (Terminology comes from Wikipedia diagram)
        // https://commons.wikimedia.org/wiki/File:BitonicSort.svg#/media/File:BitonicSort.svg
        SizeType boxSize = SizeType(2) << (i - j);
        // Which box contains this thread?
        SizeType boxID     = k >> (i - j);          // k * 2 / boxSize;
        SizeType boxStart  = boxID << (1 + i - j);  // boxID * boxSize
        SizeType boxOffset = k - (boxStart >> 1);   // k - boxID * boxSize / 2;
        SizeType elem1     = boxStart + boxOffset;
        // In first phase (j == 0, brown box): within a box, compare with the
        // opposite value in the box.
        // In later phases (j > 0, pink box): within a box, compare with fixed
        // distance (boxSize / 2) apart.
        SizeType elem2 = (j == 0) ? (boxStart + boxSize - 1 - boxOffset)
                                  : (elem1 + boxSize / 2);
        if (elem2 < n) {
          KeyType key1 = keyView(elem1);
          KeyType key2 = keyView(elem2);
          if (comp(key2, key1)) {
            keyView(elem1) = key2;
            keyView(elem2) = key1;
            if constexpr (!std::is_same_v<ValueViewType, std::nullptr_t>) {
              Kokkos::Experimental::swap(valueView(elem1), valueView(elem2));
            }
          }
        }
      });
      Range::barrier(t);
    }
  }
}

}  // namespace Impl

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
