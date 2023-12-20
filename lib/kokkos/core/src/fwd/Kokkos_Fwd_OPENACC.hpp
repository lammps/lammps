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
