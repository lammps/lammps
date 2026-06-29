// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_HPX_MDRANGEPOLICY_HPP_
#define KOKKOS_HPX_MDRANGEPOLICY_HPP_

#include <KokkosExp_MDRangePolicy.hpp>

namespace Kokkos {
namespace Impl {

// Settings for TeamMDRangePolicy
template <typename Rank, TeamMDRangeThreadAndVector ThreadAndVector>
struct ThreadAndVectorNestLevel<Rank, Kokkos::Experimental::HPX,
                                ThreadAndVector>
    : HostBasedNestLevel<Rank, ThreadAndVector> {};

}  // namespace Impl
}  // namespace Kokkos
#endif
