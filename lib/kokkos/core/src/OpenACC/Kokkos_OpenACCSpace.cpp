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

#define KOKKOS_IMPL_PUBLIC_INCLUDE

#include <OpenACC/Kokkos_OpenACC.hpp>
#include <OpenACC/Kokkos_OpenACCSpace.hpp>
#include <impl/Kokkos_MemorySpace.hpp>
#include <impl/Kokkos_Profiling_Interface.hpp>

#include <openacc.h>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
void *Kokkos::Experimental::OpenACCSpace::allocate(
    const Kokkos::Experimental::OpenACC &exec_space,
    const size_t arg_alloc_size) const {
  return allocate(exec_space, "[unlabeled]", arg_alloc_size);
}

void *Kokkos::Experimental::OpenACCSpace::allocate(
    const size_t arg_alloc_size) const {
  return allocate("[unlabeled]", arg_alloc_size);
}

void *Kokkos::Experimental::OpenACCSpace::allocate(
    const Kokkos::Experimental::OpenACC &exec_space, const char *arg_label,
    const size_t arg_alloc_size, const size_t arg_logical_size) const {
  return impl_allocate(exec_space, arg_label, arg_alloc_size, arg_logical_size);
}

void *Kokkos::Experimental::OpenACCSpace::allocate(
    const char *arg_label, const size_t arg_alloc_size,
    const size_t arg_logical_size) const {
  return impl_allocate(arg_label, arg_alloc_size, arg_logical_size);
}

void *Kokkos::Experimental::OpenACCSpace::impl_allocate(
    const Kokkos::Experimental::OpenACC &exec_space, const char *arg_label,
    const size_t arg_alloc_size, const size_t arg_logical_size,
    const Kokkos::Tools::SpaceHandle arg_handle) const {
  static_assert(sizeof(void *) == sizeof(uintptr_t),
                "Error sizeof(void*) != sizeof(uintptr_t)");

  void *ptr = nullptr;

  // FIXME_OPENACC multiple device instances are not yet supported, and thus
  // exec_space is ignored for now.
  (void)exec_space;

  ptr = acc_malloc(arg_alloc_size);

  if (Kokkos::Profiling::profileLibraryLoaded()) {
    const size_t reported_size =
        (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;
    Kokkos::Profiling::allocateData(arg_handle, arg_label, ptr, reported_size);
  }

  return ptr;
}

void *Kokkos::Experimental::OpenACCSpace::impl_allocate(
    const char *arg_label, const size_t arg_alloc_size,
    const size_t arg_logical_size,
    const Kokkos::Tools::SpaceHandle arg_handle) const {
  static_assert(sizeof(void *) == sizeof(uintptr_t),
                "Error sizeof(void*) != sizeof(uintptr_t)");

  void *ptr = nullptr;

  //[DEBUG] Disabled due to the synchronous behavior of the current
  // implementation.
  /*
    OpenACC::impl_static_fence(
        "Kokkos::OpenACCSpace::impl_allocate: Pre OpenACC Allocation");
  */

  ptr = acc_malloc(arg_alloc_size);

  //[DEBUG] Disabled due to the synchronous behavior of the current
  // implementation.
  /*
    OpenACC::impl_static_fence(
        "Kokkos::OpenACCSpace::impl_allocate: Post OpenACC Allocation");
  */
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    const size_t reported_size =
        (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;
    Kokkos::Profiling::allocateData(arg_handle, arg_label, ptr, reported_size);
  }

  return ptr;
}

void Kokkos::Experimental::OpenACCSpace::deallocate(
    void *const arg_alloc_ptr, const size_t arg_alloc_size) const {
  deallocate("[unlabeled]", arg_alloc_ptr, arg_alloc_size);
}

void Kokkos::Experimental::OpenACCSpace::deallocate(
    const char *arg_label, void *const arg_alloc_ptr,
    const size_t arg_alloc_size, const size_t arg_logical_size) const {
  impl_deallocate(arg_label, arg_alloc_ptr, arg_alloc_size, arg_logical_size);
}

void Kokkos::Experimental::OpenACCSpace::impl_deallocate(
    const char *arg_label, void *const arg_alloc_ptr,
    const size_t arg_alloc_size, const size_t arg_logical_size,
    const Kokkos::Tools::SpaceHandle arg_handle) const {
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    const size_t reported_size =
        (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;
    Kokkos::Profiling::deallocateData(arg_handle, arg_label, arg_alloc_ptr,
                                      reported_size);
  }

  if (arg_alloc_ptr) {
    acc_free(arg_alloc_ptr);
  }
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#ifdef KOKKOS_ENABLE_DEBUG
Kokkos::Impl::SharedAllocationRecord<void, void> SharedAllocationRecord<
    Kokkos::Experimental::OpenACCSpace, void>::s_root_record;
#endif

Kokkos::Impl::SharedAllocationRecord<Kokkos::Experimental::OpenACCSpace,
                                     void>::~SharedAllocationRecord() {
  m_space.deallocate(m_label.c_str(),
                     SharedAllocationRecord<void, void>::m_alloc_ptr,
                     (SharedAllocationRecord<void, void>::m_alloc_size -
                      sizeof(SharedAllocationHeader)));
}

Kokkos::Impl::SharedAllocationRecord<Kokkos::Experimental::OpenACCSpace, void>::
    SharedAllocationRecord(
        const Kokkos::Experimental::OpenACCSpace &arg_space,
        const std::string &arg_label, const size_t arg_alloc_size,
        const SharedAllocationRecord<void, void>::function_type arg_dealloc)
    // Pass through allocated [ SharedAllocationHeader , user_memory ]
    // Pass through deallocation function
    : base_t(
#ifdef KOKKOS_ENABLE_DEBUG
          &SharedAllocationRecord<Kokkos::Experimental::OpenACCSpace,
                                  void>::s_root_record,
#endif
          Impl::checked_allocation_with_header(arg_space, arg_label,
                                               arg_alloc_size),
          sizeof(SharedAllocationHeader) + arg_alloc_size, arg_dealloc,
          arg_label),
      m_space(arg_space) {
  SharedAllocationHeader header;

  this->base_t::_fill_host_accessible_header_info(header, arg_label);

  Kokkos::Impl::DeepCopy<Experimental::OpenACCSpace, HostSpace>(
      RecordBase::m_alloc_ptr, &header, sizeof(SharedAllocationHeader));
  Kokkos::fence(
      "SharedAllocationRecord<Kokkos::Experimental::OpenACCSpace, "
      "void>::SharedAllocationRecord(): fence after copying header from "
      "HostSpace");
}

//==============================================================================
// <editor-fold desc="Explicit instantiations of CRTP Base classes"> {{{1

#include <impl/Kokkos_SharedAlloc_timpl.hpp>

// To avoid additional compilation cost for something that's (mostly?) not
// performance sensitive, we explicitly instantiate these CRTP base classes
// here, where we have access to the associated *_timpl.hpp header files.
template class Kokkos::Impl::HostInaccessibleSharedAllocationRecordCommon<
    Kokkos::Experimental::OpenACCSpace>;
template class Kokkos::Impl::SharedAllocationRecordCommon<
    Kokkos::Experimental::OpenACCSpace>;

// </editor-fold> end Explicit instantiations of CRTP Base classes }}}1
//==============================================================================
