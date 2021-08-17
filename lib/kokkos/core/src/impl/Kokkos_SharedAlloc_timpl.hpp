/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (12/8/20) National Technology & Engineering
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

#ifndef KOKKOS_IMPL_SHAREDALLOC_TIMPL_HPP
#define KOKKOS_IMPL_SHAREDALLOC_TIMPL_HPP

#include <Kokkos_Macros.hpp>
#include <Kokkos_Core_fwd.hpp>

#include <impl/Kokkos_SharedAlloc.hpp>

#include <Kokkos_HostSpace.hpp>  // used with HostInaccessible specializations

#include <string>    // std::string
#include <cstring>   // strncpy
#include <iostream>  // ostream

namespace Kokkos {
namespace Impl {

template <class MemorySpace>
auto SharedAllocationRecordCommon<MemorySpace>::allocate(
    MemorySpace const& arg_space, std::string const& arg_label,
    size_t arg_alloc_size) -> derived_t* {
  return new derived_t(arg_space, arg_label, arg_alloc_size);
}

template <class MemorySpace>
void* SharedAllocationRecordCommon<MemorySpace>::allocate_tracked(
    const MemorySpace& arg_space, const std::string& arg_alloc_label,
    size_t arg_alloc_size) {
  if (!arg_alloc_size) return nullptr;

  SharedAllocationRecord* const r =
      allocate(arg_space, arg_alloc_label, arg_alloc_size);

  record_base_t::increment(r);

  return r->data();
}

template <class MemorySpace>
void SharedAllocationRecordCommon<MemorySpace>::deallocate(
    SharedAllocationRecordCommon::record_base_t* arg_rec) {
  delete static_cast<derived_t*>(arg_rec);
}

template <class MemorySpace>
void SharedAllocationRecordCommon<MemorySpace>::deallocate_tracked(
    void* arg_alloc_ptr) {
  if (arg_alloc_ptr != nullptr) {
    SharedAllocationRecord* const r = derived_t::get_record(arg_alloc_ptr);
    record_base_t::decrement(r);
  }
}

template <class MemorySpace>
void* SharedAllocationRecordCommon<MemorySpace>::reallocate_tracked(
    void* arg_alloc_ptr, size_t arg_alloc_size) {
  derived_t* const r_old = derived_t::get_record(arg_alloc_ptr);
  derived_t* const r_new =
      allocate(r_old->m_space, r_old->get_label(), arg_alloc_size);

  Kokkos::Impl::DeepCopy<MemorySpace, MemorySpace>(
      r_new->data(), r_old->data(), std::min(r_old->size(), r_new->size()));

  record_base_t::increment(r_new);
  record_base_t::decrement(r_old);

  return r_new->data();
}

template <class MemorySpace>
auto SharedAllocationRecordCommon<MemorySpace>::get_record(void* alloc_ptr)
    -> derived_t* {
  using Header = SharedAllocationHeader;

  Header const* const h = alloc_ptr ? Header::get_header(alloc_ptr) : nullptr;

  if (!alloc_ptr || h->m_record->m_alloc_ptr != h) {
    Kokkos::Impl::throw_runtime_exception(
        std::string("Kokkos::Impl::SharedAllocationRecordCommon<") +
        std::string(MemorySpace::name()) +
        std::string(">::get_record() ERROR"));
  }

  return static_cast<derived_t*>(h->m_record);
}

template <class MemorySpace>
std::string SharedAllocationRecordCommon<MemorySpace>::get_label() const {
  return std::string(record_base_t::head()->m_label);
}

template <class MemorySpace>
void SharedAllocationRecordCommon<MemorySpace>::
    _fill_host_accessible_header_info(SharedAllocationHeader& arg_header,
                                      std::string const& arg_label) {
  // Fill in the Header information, directly accessible on the host

  arg_header.m_record = &self();

  strncpy(arg_header.m_label, arg_label.c_str(),
          SharedAllocationHeader::maximum_label_length);
  // Set last element zero, in case c_str is too long
  arg_header.m_label[SharedAllocationHeader::maximum_label_length - 1] = '\0';
}

template <class MemorySpace>
void SharedAllocationRecordCommon<MemorySpace>::print_records(
    std::ostream& s, const MemorySpace&, bool detail) {
  (void)s;
  (void)detail;
#ifdef KOKKOS_ENABLE_DEBUG
  SharedAllocationRecord<void, void>::print_host_accessible_records(
      s, MemorySpace::name(), &derived_t::s_root_record, detail);
#else
  Kokkos::Impl::throw_runtime_exception(
      std::string("SharedAllocationHeader<") +
      std::string(MemorySpace::name()) +
      std::string(
          ">::print_records only works with KOKKOS_ENABLE_DEBUG enabled"));
#endif
}

template <class MemorySpace>
void HostInaccessibleSharedAllocationRecordCommon<MemorySpace>::print_records(
    std::ostream& s, const MemorySpace&, bool detail) {
  (void)s;
  (void)detail;
#ifdef KOKKOS_ENABLE_DEBUG
  SharedAllocationRecord<void, void>* r = &derived_t::s_root_record;

  char buffer[256];

  SharedAllocationHeader head;

  if (detail) {
    do {
      if (r->m_alloc_ptr) {
        Kokkos::Impl::DeepCopy<HostSpace, MemorySpace>(
            &head, r->m_alloc_ptr, sizeof(SharedAllocationHeader));
      } else {
        head.m_label[0] = 0;
      }

      // Formatting dependent on sizeof(uintptr_t)
      const char* format_string;

      if (sizeof(uintptr_t) == sizeof(unsigned long)) {
        format_string =
            "%s addr( 0x%.12lx ) list( 0x%.12lx 0x%.12lx ) extent[ 0x%.12lx "
            "+ %.8ld ] count(%d) dealloc(0x%.12lx) %s\n";
      } else if (sizeof(uintptr_t) == sizeof(unsigned long long)) {
        format_string =
            "%s addr( 0x%.12llx ) list( 0x%.12llx 0x%.12llx ) extent[ "
            "0x%.12llx + %.8ld ] count(%d) dealloc(0x%.12llx) %s\n";
      }

      snprintf(buffer, 256, format_string, MemorySpace::execution_space::name(),
               reinterpret_cast<uintptr_t>(r),
               reinterpret_cast<uintptr_t>(r->m_prev),
               reinterpret_cast<uintptr_t>(r->m_next),
               reinterpret_cast<uintptr_t>(r->m_alloc_ptr), r->m_alloc_size,
               r->m_count, reinterpret_cast<uintptr_t>(r->m_dealloc),
               head.m_label);
      s << buffer;
      r = r->m_next;
    } while (r != &derived_t::s_root_record);
  } else {
    do {
      if (r->m_alloc_ptr) {
        Kokkos::Impl::DeepCopy<HostSpace, MemorySpace>(
            &head, r->m_alloc_ptr, sizeof(SharedAllocationHeader));

        // Formatting dependent on sizeof(uintptr_t)
        const char* format_string;

        if (sizeof(uintptr_t) == sizeof(unsigned long)) {
          format_string = "%s [ 0x%.12lx + %ld ] %s\n";
        } else if (sizeof(uintptr_t) == sizeof(unsigned long long)) {
          format_string = "%s [ 0x%.12llx + %ld ] %s\n";
        }

        snprintf(
            buffer, 256, format_string, MemorySpace::execution_space::name(),
            reinterpret_cast<uintptr_t>(r->data()), r->size(), head.m_label);
      } else {
        snprintf(buffer, 256, "%s [ 0 + 0 ]\n",
                 MemorySpace::execution_space::name());
      }
      s << buffer;
      r = r->m_next;
    } while (r != &derived_t::s_root_record);
  }
#else
  Kokkos::Impl::throw_runtime_exception(
      std::string("SharedAllocationHeader<") +
      std::string(MemorySpace::name()) +
      std::string(
          ">::print_records only works with KOKKOS_ENABLE_DEBUG enabled"));
#endif
}

template <class MemorySpace>
auto HostInaccessibleSharedAllocationRecordCommon<MemorySpace>::get_record(
    void* alloc_ptr) -> derived_t* {
  // Copy the header from the allocation
  SharedAllocationHeader head;

  SharedAllocationHeader const* const head_cuda =
      alloc_ptr ? SharedAllocationHeader::get_header(alloc_ptr) : nullptr;

  if (alloc_ptr) {
    Kokkos::Impl::DeepCopy<HostSpace, MemorySpace>(
        &head, head_cuda, sizeof(SharedAllocationHeader));
  }

  derived_t* const record =
      alloc_ptr ? static_cast<derived_t*>(head.m_record) : nullptr;

  if (!alloc_ptr || record->m_alloc_ptr != head_cuda) {
    Kokkos::Impl::throw_runtime_exception(
        std::string("Kokkos::Impl::SharedAllocationRecord<") +
        std::string(MemorySpace::name()) +
        std::string(", void>::get_record ERROR"));
  }

  return record;
}

template <class MemorySpace>
std::string
HostInaccessibleSharedAllocationRecordCommon<MemorySpace>::get_label() const {
  SharedAllocationHeader header;

  Kokkos::Impl::DeepCopy<Kokkos::HostSpace, MemorySpace>(
      &header, this->record_base_t::head(), sizeof(SharedAllocationHeader));

  return std::string(header.m_label);
}

}  // end namespace Impl
}  // end namespace Kokkos

#endif  // KOKKOS_IMPL_SHAREDALLOC_TIMPL_HPP
