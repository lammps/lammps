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

#include <Kokkos_Atomic.hpp>
#include <Kokkos_HostSpace.hpp>
#include <impl/Kokkos_Error.hpp>
#include <impl/Kokkos_Tools.hpp>

#include <cstddef>
#include <cstdlib>
#include <cstdint>
#include <cstring>

#include <iostream>
#include <sstream>
#include <cstring>

#ifdef KOKKOS_COMPILER_INTEL
#include <aligned_new>
#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
KOKKOS_DEPRECATED HostSpace::HostSpace(const HostSpace::AllocationMechanism &)
    : HostSpace() {}
#endif

void *HostSpace::allocate(const size_t arg_alloc_size) const {
  return allocate("[unlabeled]", arg_alloc_size);
}
void *HostSpace::allocate(const char *arg_label, const size_t arg_alloc_size,
                          const size_t

                              arg_logical_size) const {
  return impl_allocate(arg_label, arg_alloc_size, arg_logical_size);
}
void *HostSpace::impl_allocate(
    const char *arg_label, const size_t arg_alloc_size,
    const size_t arg_logical_size,
    const Kokkos::Tools::SpaceHandle arg_handle) const {
  const size_t reported_size =
      (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;
  static_assert(sizeof(void *) == sizeof(uintptr_t),
                "Error sizeof(void*) != sizeof(uintptr_t)");

  static_assert(
      Kokkos::Impl::is_integral_power_of_two(Kokkos::Impl::MEMORY_ALIGNMENT),
      "Memory alignment must be power of two");

  constexpr uintptr_t alignment      = Kokkos::Impl::MEMORY_ALIGNMENT;
  constexpr uintptr_t alignment_mask = alignment - 1;

  void *ptr = nullptr;

  if (arg_alloc_size)
    ptr = operator new (arg_alloc_size, std::align_val_t(alignment),
                        std::nothrow_t{});

  if (!ptr || (reinterpret_cast<uintptr_t>(ptr) == ~uintptr_t(0)) ||
      (reinterpret_cast<uintptr_t>(ptr) & alignment_mask)) {
    Impl::throw_bad_alloc(name(), arg_alloc_size, arg_label);
  }
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    Kokkos::Profiling::allocateData(arg_handle, arg_label, ptr, reported_size);
  }
  return ptr;
}

void HostSpace::deallocate(void *const arg_alloc_ptr,
                           const size_t arg_alloc_size) const {
  deallocate("[unlabeled]", arg_alloc_ptr, arg_alloc_size);
}

void HostSpace::deallocate(const char *arg_label, void *const arg_alloc_ptr,
                           const size_t arg_alloc_size,
                           const size_t arg_logical_size) const {
  if (arg_alloc_ptr) Kokkos::fence("HostSpace::impl_deallocate before free");
  impl_deallocate(arg_label, arg_alloc_ptr, arg_alloc_size, arg_logical_size);
}
void HostSpace::impl_deallocate(
    const char *arg_label, void *const arg_alloc_ptr,
    const size_t arg_alloc_size, const size_t arg_logical_size,
    const Kokkos::Tools::SpaceHandle arg_handle) const {
  if (arg_alloc_ptr) {
    size_t reported_size =
        (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;
    if (Kokkos::Profiling::profileLibraryLoaded()) {
      Kokkos::Profiling::deallocateData(arg_handle, arg_label, arg_alloc_ptr,
                                        reported_size);
    }
    constexpr uintptr_t alignment = Kokkos::Impl::MEMORY_ALIGNMENT;
    operator delete (arg_alloc_ptr, std::align_val_t(alignment),
                     std::nothrow_t{});
  }
}

}  // namespace Kokkos

#include <impl/Kokkos_SharedAlloc_timpl.hpp>

KOKKOS_IMPL_SHARED_ALLOCATION_RECORD_EXPLICIT_INSTANTIATION(Kokkos::HostSpace);
