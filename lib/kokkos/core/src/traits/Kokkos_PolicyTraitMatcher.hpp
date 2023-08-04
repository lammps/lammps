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

#include <impl/Kokkos_Utilities.hpp>  // type_list

#include <traits/Kokkos_Traits_fwd.hpp>

#ifndef KOKKOS_KOKKOS_POLICYTRAITMATCHER_HPP
#define KOKKOS_KOKKOS_POLICYTRAITMATCHER_HPP

namespace Kokkos {
namespace Impl {

//==============================================================================
// <editor-fold desc="PolicyTraitMatcher"> {{{1

// To handle the WorkTag case, we need more than just a predicate; we need
// something that we can default to in the unspecialized case, just like we
// do for AnalyzeExecPolicy
template <class TraitSpec, class Trait, class Enable = void>
struct PolicyTraitMatcher : std::false_type {};

template <class TraitSpec, class Trait>
struct PolicyTraitMatcher<
    TraitSpec, Trait,
    std::enable_if_t<
        TraitSpec::template trait_matches_specification<Trait>::value>>
    : std::true_type {};

// </editor-fold> end PolicyTraitMatcher }}}1
//==============================================================================

}  // end namespace Impl
}  // end namespace Kokkos

#endif  // KOKKOS_KOKKOS_POLICYTRAITMATCHER_HPP
