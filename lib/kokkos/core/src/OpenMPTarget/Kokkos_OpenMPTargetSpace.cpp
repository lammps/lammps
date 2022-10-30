/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
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

#include <Kokkos_OpenMPTarget.hpp>
#include <Kokkos_OpenMPTargetSpace.hpp>
#include <impl/Kokkos_Error.hpp>
#include <Kokkos_Atomic.hpp>
#include <impl/Kokkos_MemorySpace.hpp>

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

  void* ptr;

  ptr = omp_target_alloc(arg_alloc_size, omp_get_default_device());

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

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

#ifdef KOKKOS_ENABLE_DEBUG
SharedAllocationRecord<void, void> SharedAllocationRecord<
    Kokkos::Experimental::OpenMPTargetSpace, void>::s_root_record;
#endif

SharedAllocationRecord<Kokkos::Experimental::OpenMPTargetSpace,
                       void>::~SharedAllocationRecord() {
  auto alloc_size = SharedAllocationRecord<void, void>::m_alloc_size;
  m_space.deallocate(m_label.c_str(),
                     SharedAllocationRecord<void, void>::m_alloc_ptr,
                     alloc_size, (alloc_size - sizeof(SharedAllocationHeader)));
}

SharedAllocationRecord<Kokkos::Experimental::OpenMPTargetSpace, void>::
    SharedAllocationRecord(
        const Kokkos::Experimental::OpenMPTargetSpace& arg_space,
        const std::string& arg_label, const size_t arg_alloc_size,
        const SharedAllocationRecord<void, void>::function_type arg_dealloc)
    // Pass through allocated [ SharedAllocationHeader , user_memory ]
    // Pass through deallocation function
    : base_t(
#ifdef KOKKOS_ENABLE_DEBUG
          &SharedAllocationRecord<Kokkos::Experimental::OpenMPTargetSpace,
                                  void>::s_root_record,
#endif
          Kokkos::Impl::checked_allocation_with_header(arg_space, arg_label,
                                                       arg_alloc_size),
          sizeof(SharedAllocationHeader) + arg_alloc_size, arg_dealloc,
          arg_label),
      m_space(arg_space) {
  SharedAllocationHeader header;

  this->base_t::_fill_host_accessible_header_info(header, arg_label);

  // TODO DeepCopy
  // DeepCopy
  Kokkos::Impl::DeepCopy<Experimental::OpenMPTargetSpace, HostSpace>(
      RecordBase::m_alloc_ptr, &header, sizeof(SharedAllocationHeader));
  Kokkos::fence(
      "SharedAllocationRecord<Kokkos::Experimental::OpenMPTargetSpace, "
      "void>::SharedAllocationRecord(): fence after copying header from "
      "HostSpace");
}

//----------------------------------------------------------------------------

}  // namespace Impl
}  // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*
namespace Kokkos {
namespace {
  const unsigned HOST_SPACE_ATOMIC_MASK = 0xFFFF;
  const unsigned HOST_SPACE_ATOMIC_XOR_MASK = 0x5A39;
  static int HOST_SPACE_ATOMIC_LOCKS[HOST_SPACE_ATOMIC_MASK+1];
}

namespace Impl {
void init_lock_array_host_space() {
  static int is_initialized = 0;
  if(! is_initialized)
    for(int i = 0; i < static_cast<int> (HOST_SPACE_ATOMIC_MASK+1); i++)
      HOST_SPACE_ATOMIC_LOCKS[i] = 0;
}

bool lock_address_host_space(void* ptr) {
  return 0 == atomic_compare_exchange( &HOST_SPACE_ATOMIC_LOCKS[
      (( size_t(ptr) >> 2 ) & HOST_SPACE_ATOMIC_MASK) ^
HOST_SPACE_ATOMIC_XOR_MASK] , 0 , 1);
}

void unlock_address_host_space(void* ptr) {
   atomic_exchange( &HOST_SPACE_ATOMIC_LOCKS[
      (( size_t(ptr) >> 2 ) & HOST_SPACE_ATOMIC_MASK) ^
HOST_SPACE_ATOMIC_XOR_MASK] , 0);
}

}
}*/

//==============================================================================
// <editor-fold desc="Explicit instantiations of CRTP Base classes"> {{{1

#include <impl/Kokkos_SharedAlloc_timpl.hpp>

namespace Kokkos {
namespace Impl {

// To avoid additional compilation cost for something that's (mostly?) not
// performance sensitive, we explicity instantiate these CRTP base classes here,
// where we have access to the associated *_timpl.hpp header files.
template class HostInaccessibleSharedAllocationRecordCommon<
    Kokkos::Experimental::OpenMPTargetSpace>;
template class SharedAllocationRecordCommon<
    Kokkos::Experimental::OpenMPTargetSpace>;

}  // end namespace Impl
}  // end namespace Kokkos

// </editor-fold> end Explicit instantiations of CRTP Base classes }}}1
//==============================================================================
