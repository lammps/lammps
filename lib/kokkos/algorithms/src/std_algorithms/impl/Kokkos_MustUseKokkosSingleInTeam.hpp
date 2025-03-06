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

#ifndef KOKKOS_STD_ALGORITHMS_MUSTUSEKOKKOSSINGLEINTEAM_HPP
#define KOKKOS_STD_ALGORITHMS_MUSTUSEKOKKOSSINGLEINTEAM_HPP

#include <Kokkos_Core.hpp>

namespace Kokkos {
namespace Experimental {
namespace Impl {

template <typename T>
struct stdalgo_must_use_kokkos_single_for_team_scan : std::false_type {};

// the following do not support the overload for team-level scan
// accepting an "out" value to store the scan result

// FIXME_OPENACC
#if defined(KOKKOS_ENABLE_OPENACC)
template <>
struct stdalgo_must_use_kokkos_single_for_team_scan<
    Kokkos::Experimental::OpenACC> : std::true_type {};
#endif

template <typename T>
inline constexpr bool stdalgo_must_use_kokkos_single_for_team_scan_v =
    stdalgo_must_use_kokkos_single_for_team_scan<T>::value;

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos

#endif
