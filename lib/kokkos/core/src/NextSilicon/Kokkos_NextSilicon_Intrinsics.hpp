// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

// FIXME_NEXTSILICON: placeholder intrinsic declaration (tracked by NextSilicon
// ticket CS-694)
/// An internal-use-only intrinsic used to communicate from Kokkos C++ code to
/// the MLIR compiler stack that a struct at the address `parameter_struct` can
/// be considered immutable and thread invariant for the full duration of the
/// microtask, and that consequently all its members are read only, and that the
/// pointer to the struct and any derived pointer does not alias with anything
/// else (pointers stored in the struct can still alias). This should _only_ be
/// used to annotate the functor objects passed into Kokkos parallel constructs,
/// and should not be used _anywhere_ outside of Kokkos.
namespace Kokkos::Experimental::Impl {
extern "C" void __ns_immutable_thread_invariant_parameter_struct(
    const void* parameter_struct);
}
