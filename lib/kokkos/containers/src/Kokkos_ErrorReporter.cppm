// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

module;

#include <Kokkos_ErrorReporter.hpp>

export module kokkos.error_reporter;

export {
  namespace Kokkos {
  namespace Experimental {
  using ::Kokkos::Experimental::ErrorReporter;
  }
  }  // namespace Kokkos
}
