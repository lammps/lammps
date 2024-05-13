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

// Experimental unified task-data parallel manycore LDRD

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#include <Kokkos_Macros.hpp>
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#endif
#ifndef KOKKOS_IMPL_POINTEROWNERSHIP_HPP
#define KOKKOS_IMPL_POINTEROWNERSHIP_HPP

#include <Kokkos_Macros.hpp>

#include <Kokkos_Core_fwd.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

/// Trivial wrapper for raw pointers that express ownership.
template <class T>
using OwningRawPtr = T*;

/// Trivial wrapper for raw pointers that do not express ownership.
template <class T>
using ObservingRawPtr = T*;

}  // end namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_IMPL_POINTEROWNERSHIP_HPP */
