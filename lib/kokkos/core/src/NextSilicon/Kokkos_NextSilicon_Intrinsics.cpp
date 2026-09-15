// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <NextSilicon/Kokkos_NextSilicon_Intrinsics.hpp>

extern "C" void
Kokkos::Experimental::Impl::__ns_immutable_thread_invariant_parameter_struct(
    const void*) {
  // FIXME_NEXTSILICON: this thing purposefully left empty to be filled in by NS
  // toolchain
}
