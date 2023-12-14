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

#ifndef KOKKOS_HBWSPACE_FWD_HPP_
#define KOKKOS_HBWSPACE_FWD_HPP_

#ifdef KOKKOS_ENABLE_HBWSPACE
namespace Kokkos {

namespace Experimental {
class HBWSpace;  /// Memory space for hbw_malloc from memkind (e.g. for KNL
                 /// processor)
}  // namespace Experimental
}  // namespace Kokkos
#endif
#endif
