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
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#endif

#include <Kokkos_Macros.hpp>

#include <algorithm>
#include <omp.h>

/*--------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdint.h>
#include <memory.h>

#include <iostream>
#include <sstream>
#include <cstring>

#include <OpenMPTarget/Kokkos_OpenMPTarget.hpp>
#include <OpenMPTarget/Kokkos_OpenMPTargetSpace.hpp>
#include <impl/Kokkos_Error.hpp>
#include <Kokkos_Atomic.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {
/* Default allocation mechanism */
OpenMPTargetSpace::OpenMPTargetSpace() {}

void* OpenMPTargetSpace::impl_allocate(

    const char* arg_label, const size_t arg_alloc_size,
    const size_t arg_logical_size,
    const Kokkos::Tools::SpaceHandle arg_handle) const {
  static_assert(sizeof(void*) == sizeof(uintptr_t),
                "Error sizeof(void*) != sizeof(uintptr_t)");

  void* ptr = omp_target_alloc(arg_alloc_size, omp_get_default_device());

  if (!ptr) {
    Kokkos::Impl::throw_bad_alloc(name(), arg_alloc_size, arg_label);
  }

  if (Kokkos::Profiling::profileLibraryLoaded()) {
    const size_t reported_size =
        (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;
    Kokkos::Profiling::allocateData(arg_handle, arg_label, ptr, reported_size);
  }

  return ptr;
}

void* OpenMPTargetSpace::allocate(const size_t arg_alloc_size) const {
  return allocate("[unlabeled]", arg_alloc_size);
}

void* OpenMPTargetSpace::allocate(const char* arg_label,
                                  const size_t arg_alloc_size,
                                  const size_t arg_logical_size) const {
  return impl_allocate(arg_label, arg_alloc_size, arg_logical_size);
}

void OpenMPTargetSpace::impl_deallocate(
    const char* arg_label, void* const arg_alloc_ptr,
    const size_t arg_alloc_size, const size_t arg_logical_size,
    const Kokkos::Tools::SpaceHandle arg_handle) const {
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    const size_t reported_size =
        (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;
    Kokkos::Profiling::deallocateData(arg_handle, arg_label, arg_alloc_ptr,
                                      reported_size);
  }
  if (arg_alloc_ptr) {
    omp_target_free(arg_alloc_ptr, omp_get_default_device());
  }
}

void OpenMPTargetSpace::deallocate(void* const arg_alloc_ptr,
                                   const size_t arg_alloc_size) const {
  deallocate("[unlabeled]", arg_alloc_ptr, arg_alloc_size);
}

void OpenMPTargetSpace::deallocate(const char* arg_label,
                                   void* const arg_alloc_ptr,
                                   const size_t arg_alloc_size,
                                   const size_t arg_logical_size) const

{
  impl_deallocate(arg_label, arg_alloc_ptr, arg_alloc_size, arg_logical_size);
}

}  // namespace Experimental
}  // namespace Kokkos

//==============================================================================
// <editor-fold desc="Explicit instantiations of CRTP Base classes"> {{{1

#include <impl/Kokkos_SharedAlloc_timpl.hpp>

KOKKOS_IMPL_HOST_INACCESSIBLE_SHARED_ALLOCATION_RECORD_EXPLICIT_INSTANTIATION(
    Kokkos::Experimental::OpenMPTargetSpace);

// </editor-fold> end Explicit instantiations of CRTP Base classes }}}1
//==============================================================================
