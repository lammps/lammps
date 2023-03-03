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

/** @file Kokkos_MemorySpace.hpp
 *
 *  Operations common to memory space instances, or at least default
 *  implementations thereof.
 */

#ifndef KOKKOS_IMPL_MEMORYSPACE_HPP
#define KOKKOS_IMPL_MEMORYSPACE_HPP

#include <Kokkos_Macros.hpp>
#include <impl/Kokkos_SharedAlloc.hpp>
#include <impl/Kokkos_Error.hpp>

#include <string>

namespace Kokkos {
namespace Impl {

// Defined in implementation file to avoid having to include iostream
void safe_throw_allocation_with_header_failure(
    std::string const &space_name, std::string const &label,
    Kokkos::Experimental::RawMemoryAllocationFailure const &failure);

template <class MemorySpace>
SharedAllocationHeader *checked_allocation_with_header(MemorySpace const &space,
                                                       std::string const &label,
                                                       size_t alloc_size) {
  try {
    return reinterpret_cast<SharedAllocationHeader *>(space.allocate(
        label.c_str(), alloc_size + sizeof(SharedAllocationHeader),
        alloc_size));
  } catch (Kokkos::Experimental::RawMemoryAllocationFailure const &failure) {
    safe_throw_allocation_with_header_failure(space.name(), label, failure);
  }
  return nullptr;  // unreachable
}

template <class ExecutionSpace, class MemorySpace>
SharedAllocationHeader *checked_allocation_with_header(
    ExecutionSpace const &exec_space, MemorySpace const &space,
    std::string const &label, size_t alloc_size) {
  try {
    return reinterpret_cast<SharedAllocationHeader *>(space.allocate(
        exec_space, label.c_str(), alloc_size + sizeof(SharedAllocationHeader),
        alloc_size));
  } catch (Kokkos::Experimental::RawMemoryAllocationFailure const &failure) {
    safe_throw_allocation_with_header_failure(space.name(), label, failure);
  }
  return nullptr;  // unreachable
}

}  // end namespace Impl
}  // end namespace Kokkos

#endif  // KOKKOS_IMPL_MEMORYSPACE_HPP
