// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_STD_ALGORITHMS_MUSTUSEKOKKOSSINGLEINTEAM_HPP
#define KOKKOS_STD_ALGORITHMS_MUSTUSEKOKKOSSINGLEINTEAM_HPP

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif

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
