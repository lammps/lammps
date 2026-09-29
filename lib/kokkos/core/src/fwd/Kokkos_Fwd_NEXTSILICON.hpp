// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_NEXTSILICON_FWD_HPP_
#define KOKKOS_NEXTSILICON_FWD_HPP_

#if defined(KOKKOS_ENABLE_NEXTSILICON)
namespace Kokkos::Experimental {
class NextSilicon;             ///< NextSilicon execution space.
class NextSiliconSharedSpace;  ///< Memory migratable between Host and
                               ///< NextSilicon accelerator
}  // namespace Kokkos::Experimental
#endif
#endif
