// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

// Intentionally provides no Kokkos profiling interface symbols.
// Used to test that loading a library with no callbacks produces an error.

extern "C" {
void not_a_profiling_callback() {}
}
