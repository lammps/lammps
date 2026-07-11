// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_OPENACC_FWD_HPP_
#define KOKKOS_OPENACC_FWD_HPP_

#if defined(KOKKOS_ENABLE_OPENACC)
namespace Kokkos {
namespace Experimental {
class OpenACC;  ///< OpenACC execution space.
class OpenACCSpace;
}  // namespace Experimental
}  // namespace Kokkos
#endif
#endif
