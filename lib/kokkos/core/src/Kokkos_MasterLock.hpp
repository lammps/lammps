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

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#include <Kokkos_Macros.hpp>
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#endif
#ifndef KOKKOS_MASTER_LOCK_HPP
#define KOKKOS_MASTER_LOCK_HPP

#include <Kokkos_Macros.hpp>

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3

namespace Kokkos {
namespace Experimental {

// my be used to coordinate work between master instances
// SHOULD NOT be used within a parallel algorithm
//
// This lock should be used with with a scoped lock guard
// i.e. std::unique_lock<Lock>, std::lock_guard
//
// cannot be copied or moved
// has the following functions available
//
// Lock()
// ~Lock()
//
// void lock()
// void unlock()
// bool try_lock()
//
template <typename ExecutionSpace>
class MasterLock;

}  // namespace Experimental
}  // namespace Kokkos

#endif

#endif  // KOKKOS_MASTER_LOCK_HPP
