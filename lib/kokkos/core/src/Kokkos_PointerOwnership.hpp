// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

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

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
namespace Kokkos {

/// Trivial wrapper for raw pointers that express ownership.
template <class T>
using OwningRawPtr KOKKOS_DEPRECATED = T*;

/// Trivial wrapper for raw pointers that do not express ownership.
template <class T>
using ObservingRawPtr KOKKOS_DEPRECATED = T*;

}  // end namespace Kokkos
#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_IMPL_POINTEROWNERSHIP_HPP */
