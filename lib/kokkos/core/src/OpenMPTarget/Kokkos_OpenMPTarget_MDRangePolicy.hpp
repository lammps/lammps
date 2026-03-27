// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_OPENMPTARGET_MDRANGEPOLICY_HPP_
#define KOKKOS_OPENMPTARGET_MDRANGEPOLICY_HPP_

#include <KokkosExp_MDRangePolicy.hpp>

namespace Kokkos {
namespace Impl {

using OpenMPTargetIterateLeft = std::integral_constant<Iterate, Iterate::Left>;
using OpenMPTargetIterateRight =
    std::integral_constant<Iterate, Iterate::Right>;

template <typename Rank,
          ::Kokkos::Impl::TeamMDRangeThreadAndVector ThreadAndVector>
struct ThreadAndVectorNestLevel<Rank, Kokkos::Experimental::OpenMPTarget,
                                ThreadAndVector>
    : AcceleratorBasedNestLevel<Rank, ThreadAndVector> {};

}  // namespace Impl
}  // namespace Kokkos

#endif
