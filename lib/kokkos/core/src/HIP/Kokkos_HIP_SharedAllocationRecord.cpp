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

#include <HIP/Kokkos_HIP_SharedAllocationRecord.hpp>
#include <HIP/Kokkos_HIP_DeepCopy.hpp>
#include <HIP/Kokkos_HIP.hpp>

namespace Kokkos {
namespace Impl {

#ifdef KOKKOS_ENABLE_DEBUG
SharedAllocationRecord<void, void>
    SharedAllocationRecord<HIPSpace, void>::s_root_record;

SharedAllocationRecord<void, void>
    SharedAllocationRecord<HIPHostPinnedSpace, void>::s_root_record;

SharedAllocationRecord<void, void>
    SharedAllocationRecord<HIPManagedSpace, void>::s_root_record;
#endif

SharedAllocationRecord<HIPSpace, void>::~SharedAllocationRecord() {
  auto alloc_size = SharedAllocationRecord<void, void>::m_alloc_size;
  m_space.deallocate(m_label.c_str(),
                     SharedAllocationRecord<void, void>::m_alloc_ptr,
                     alloc_size, (alloc_size - sizeof(SharedAllocationHeader)));
}

SharedAllocationRecord<HIPHostPinnedSpace, void>::~SharedAllocationRecord() {
  m_space.deallocate(m_label.c_str(),
                     SharedAllocationRecord<void, void>::m_alloc_ptr,
                     SharedAllocationRecord<void, void>::m_alloc_size);
}

SharedAllocationRecord<HIPManagedSpace, void>::~SharedAllocationRecord() {
  m_space.deallocate(m_label.c_str(),
                     SharedAllocationRecord<void, void>::m_alloc_ptr,
                     SharedAllocationRecord<void, void>::m_alloc_size);
}

SharedAllocationRecord<HIPSpace, void>::SharedAllocationRecord(
    const HIPSpace& arg_space, const std::string& arg_label,
    const size_t arg_alloc_size,
    const SharedAllocationRecord<void, void>::function_type arg_dealloc)
    // Pass through allocated [ SharedAllocationHeader , user_memory ]
    // Pass through deallocation function
    : base_t(
#ifdef KOKKOS_ENABLE_DEBUG
          &SharedAllocationRecord<HIPSpace, void>::s_root_record,
#endif
          Kokkos::Impl::checked_allocation_with_header(arg_space, arg_label,
                                                       arg_alloc_size),
          sizeof(SharedAllocationHeader) + arg_alloc_size, arg_dealloc,
          arg_label),
      m_space(arg_space) {

  SharedAllocationHeader header;

  this->base_t::_fill_host_accessible_header_info(header, arg_label);

  // Copy to device memory
  HIP exec;
  Kokkos::Impl::DeepCopy<HIPSpace, HostSpace>(
      exec, RecordBase::m_alloc_ptr, &header, sizeof(SharedAllocationHeader));
  exec.fence(
      "SharedAllocationRecord<Kokkos::HIPSpace, "
      "void>::SharedAllocationRecord(): fence after copying header from "
      "HostSpace");
}

SharedAllocationRecord<HIPSpace, void>::SharedAllocationRecord(
    const HIP& arg_exec_space, const HIPSpace& arg_space,
    const std::string& arg_label, const size_t arg_alloc_size,
    const SharedAllocationRecord<void, void>::function_type arg_dealloc)
    // Pass through allocated [ SharedAllocationHeader , user_memory ]
    // Pass through deallocation function
    : base_t(
#ifdef KOKKOS_ENABLE_DEBUG
          &SharedAllocationRecord<HIPSpace, void>::s_root_record,
#endif
          Kokkos::Impl::checked_allocation_with_header(arg_space, arg_label,
                                                       arg_alloc_size),
          sizeof(SharedAllocationHeader) + arg_alloc_size, arg_dealloc,
          arg_label),
      m_space(arg_space) {

  SharedAllocationHeader header;

  this->base_t::_fill_host_accessible_header_info(header, arg_label);

  // Copy to device memory
  Kokkos::Impl::DeepCopy<HIPSpace, HostSpace>(arg_exec_space,
                                              RecordBase::m_alloc_ptr, &header,
                                              sizeof(SharedAllocationHeader));
}

SharedAllocationRecord<HIPHostPinnedSpace, void>::SharedAllocationRecord(
    const HIPHostPinnedSpace& arg_space, const std::string& arg_label,
    const size_t arg_alloc_size,
    const SharedAllocationRecord<void, void>::function_type arg_dealloc)
    // Pass through allocated [ SharedAllocationHeader , user_memory ]
    // Pass through deallocation function
    : base_t(
#ifdef KOKKOS_ENABLE_DEBUG
          &SharedAllocationRecord<HIPHostPinnedSpace, void>::s_root_record,
#endif
          Kokkos::Impl::checked_allocation_with_header(arg_space, arg_label,
                                                       arg_alloc_size),
          sizeof(SharedAllocationHeader) + arg_alloc_size, arg_dealloc,
          arg_label),
      m_space(arg_space) {
  // Fill in the Header information, directly accessible via host pinned memory
  this->base_t::_fill_host_accessible_header_info(*RecordBase::m_alloc_ptr,
                                                  arg_label);
}

SharedAllocationRecord<HIPManagedSpace, void>::SharedAllocationRecord(
    const HIPManagedSpace& arg_space, const std::string& arg_label,
    const size_t arg_alloc_size,
    const SharedAllocationRecord<void, void>::function_type arg_dealloc)
    // Pass through allocated [ SharedAllocationHeader , user_memory ]
    // Pass through deallocation function
    : base_t(
#ifdef KOKKOS_ENABLE_DEBUG
          &SharedAllocationRecord<HIPManagedSpace, void>::s_root_record,
#endif
          Kokkos::Impl::checked_allocation_with_header(arg_space, arg_label,
                                                       arg_alloc_size),
          sizeof(SharedAllocationHeader) + arg_alloc_size, arg_dealloc,
          arg_label),
      m_space(arg_space) {
  // Fill in the Header information, directly accessible via managed memory
  this->base_t::_fill_host_accessible_header_info(*RecordBase::m_alloc_ptr,
                                                  arg_label);
}

}  // namespace Impl
}  // namespace Kokkos
