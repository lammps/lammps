/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//              Copyright (2019) Sandia Corporation
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

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
