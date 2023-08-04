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

#define KOKKOS_IMPL_PUBLIC_INCLUDE

#include <OpenACC/Kokkos_OpenACC_SharedAllocationRecord.hpp>
#include <OpenACC/Kokkos_OpenACC_DeepCopy.hpp>
#include <impl/Kokkos_MemorySpace.hpp>
#include <Kokkos_HostSpace.hpp>

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

Kokkos::Impl::SharedAllocationRecord<Kokkos::Experimental::OpenACCSpace, void>::
    SharedAllocationRecord(
        const Kokkos::Experimental::OpenACC &arg_exec_space,
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
          Impl::checked_allocation_with_header(arg_exec_space, arg_space,
                                               arg_label, arg_alloc_size),
          sizeof(SharedAllocationHeader) + arg_alloc_size, arg_dealloc,
          arg_label),
      m_space(arg_space) {
  SharedAllocationHeader header;

  this->base_t::_fill_host_accessible_header_info(header, arg_label);

  Kokkos::Impl::DeepCopy<Experimental::OpenACCSpace, HostSpace>(
      arg_exec_space, RecordBase::m_alloc_ptr, &header,
      sizeof(SharedAllocationHeader));
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
