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

#include <Kokkos_HostSpace.hpp>
#include <Kokkos_SYCL.hpp>
#include <Kokkos_SYCL_Space.hpp>
#include <SYCL/Kokkos_SYCL_DeepCopy.hpp>
#include <SYCL/Kokkos_SYCL_Instance.hpp>
#include <impl/Kokkos_MemorySpace.hpp>
#include <impl/Kokkos_Profiling.hpp>

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
namespace Kokkos {
namespace Impl {
namespace {
auto USM_memcpy(cl::sycl::queue& q, void* dst, const void* src, size_t n) {
  return q.memcpy(dst, src, n);
}

void USM_memcpy(Kokkos::Experimental::Impl::SYCLInternal& space, void* dst,
                const void* src, size_t n) {
  (void)USM_memcpy(*space.m_queue, dst, src, n);
}

void USM_memcpy(void* dst, const void* src, size_t n) {
  Kokkos::Experimental::Impl::SYCLInternal::singleton().m_queue->wait();
  USM_memcpy(*Kokkos::Experimental::Impl::SYCLInternal::singleton().m_queue,
             dst, src, n)
      .wait();
}
}  // namespace

DeepCopy<Kokkos::Experimental::SYCLDeviceUSMSpace,
         Kokkos::Experimental::SYCLDeviceUSMSpace, Kokkos::Experimental::SYCL>::
    DeepCopy(const Kokkos::Experimental::SYCL& instance, void* dst,
             const void* src, size_t n) {
  USM_memcpy(*instance.impl_internal_space_instance(), dst, src, n);
}

DeepCopy<Kokkos::Experimental::SYCLDeviceUSMSpace,
         Kokkos::Experimental::SYCLDeviceUSMSpace,
         Kokkos::Experimental::SYCL>::DeepCopy(void* dst, const void* src,
                                               size_t n) {
  USM_memcpy(dst, src, n);
}

DeepCopy<Kokkos::HostSpace, Kokkos::Experimental::SYCLDeviceUSMSpace,
         Kokkos::Experimental::SYCL>::DeepCopy(const Kokkos::Experimental::SYCL&
                                                   instance,
                                               void* dst, const void* src,
                                               size_t n) {
  USM_memcpy(*instance.impl_internal_space_instance(), dst, src, n);
}

DeepCopy<Kokkos::HostSpace, Kokkos::Experimental::SYCLDeviceUSMSpace,
         Kokkos::Experimental::SYCL>::DeepCopy(void* dst, const void* src,
                                               size_t n) {
  USM_memcpy(dst, src, n);
}

DeepCopy<Kokkos::Experimental::SYCLDeviceUSMSpace, Kokkos::HostSpace,
         Kokkos::Experimental::SYCL>::DeepCopy(const Kokkos::Experimental::SYCL&
                                                   instance,
                                               void* dst, const void* src,
                                               size_t n) {
  USM_memcpy(*instance.impl_internal_space_instance(), dst, src, n);
}

DeepCopy<Kokkos::Experimental::SYCLDeviceUSMSpace, Kokkos::HostSpace,
         Kokkos::Experimental::SYCL>::DeepCopy(void* dst, const void* src,
                                               size_t n) {
  USM_memcpy(dst, src, n);
}

}  // namespace Impl
}  // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Experimental {

SYCLDeviceUSMSpace::SYCLDeviceUSMSpace() : m_device(SYCL().sycl_device()) {}

void* SYCLDeviceUSMSpace::allocate(const size_t arg_alloc_size) const {
  return allocate("[unlabeled]", arg_alloc_size);
}
void* SYCLDeviceUSMSpace::allocate(const char* arg_label,
                                   const size_t arg_alloc_size,
                                   const size_t arg_logical_size) const {
  return impl_allocate(arg_label, arg_alloc_size, arg_logical_size);
}

void* SYCLDeviceUSMSpace::impl_allocate(
    const char* arg_label, const size_t arg_alloc_size,
    const size_t arg_logical_size,
    const Kokkos::Tools::SpaceHandle arg_handle) const {
  const cl::sycl::queue& queue =
      *SYCL().impl_internal_space_instance()->m_queue;
  void* const hostPtr = cl::sycl::malloc_device(arg_alloc_size, queue);

  if (hostPtr == nullptr)
    throw RawMemoryAllocationFailure(
        arg_alloc_size, 1, RawMemoryAllocationFailure::FailureMode::Unknown,
        RawMemoryAllocationFailure::AllocationMechanism::SYCLMalloc);

  if (Kokkos::Profiling::profileLibraryLoaded()) {
    const size_t reported_size =
        (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;
    Kokkos::Profiling::allocateData(arg_handle, arg_label, hostPtr,
                                    reported_size);
  }

  return hostPtr;
}

void SYCLDeviceUSMSpace::deallocate(void* const arg_alloc_ptr,
                                    const size_t arg_alloc_size) const {
  deallocate("[unlabeled]", arg_alloc_ptr, arg_alloc_size);
}
void SYCLDeviceUSMSpace::deallocate(const char* arg_label,
                                    void* const arg_alloc_ptr,
                                    const size_t arg_alloc_size,
                                    const size_t arg_logical_size) const {
  impl_deallocate(arg_label, arg_alloc_ptr, arg_alloc_size, arg_logical_size);
}
void SYCLDeviceUSMSpace::impl_deallocate(
    const char* arg_label, void* const arg_alloc_ptr,
    const size_t arg_alloc_size, const size_t arg_logical_size,
    const Kokkos::Tools::SpaceHandle arg_handle) const {
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    const size_t reported_size =
        (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;
    Kokkos::Profiling::deallocateData(arg_handle, arg_label, arg_alloc_ptr,
                                      reported_size);
  }
  const cl::sycl::queue& queue =
      *SYCL().impl_internal_space_instance()->m_queue;
  cl::sycl::free(arg_alloc_ptr, queue);
}

}  // namespace Experimental
}  // namespace Kokkos

namespace Kokkos {
namespace Impl {

#ifdef KOKKOS_ENABLE_DEBUG
SharedAllocationRecord<void, void> SharedAllocationRecord<
    Kokkos::Experimental::SYCLDeviceUSMSpace, void>::s_root_record;
#endif

SharedAllocationRecord<Kokkos::Experimental::SYCLDeviceUSMSpace, void>::
    SharedAllocationRecord(
        const Kokkos::Experimental::SYCLDeviceUSMSpace& space,
        const std::string& label, const size_t size,
        const SharedAllocationRecord<void, void>::function_type dealloc)
    // Pass through allocated [ SharedAllocationHeader , user_memory ]
    // Pass through deallocation function
    : SharedAllocationRecord<void, void>(
#ifdef KOKKOS_ENABLE_DEBUG
          &SharedAllocationRecord<Kokkos::Experimental::SYCLDeviceUSMSpace,
                                  void>::s_root_record,
#endif
          Kokkos::Impl::checked_allocation_with_header(space, label, size),
          sizeof(SharedAllocationHeader) + size, dealloc),
      m_space(space) {
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    Kokkos::Profiling::allocateData(
        Kokkos::Profiling::make_space_handle(space.name()), label, data(),
        size);
  }

  SharedAllocationHeader header;

  // Fill in the Header information
  header.m_record = static_cast<SharedAllocationRecord<void, void>*>(this);

  strncpy(header.m_label, label.c_str(),
          SharedAllocationHeader::maximum_label_length);
  // Set last element zero, in case c_str is too long
  header.m_label[SharedAllocationHeader::maximum_label_length - 1] = (char)0;

  // Copy to device memory
  Kokkos::Impl::DeepCopy<Kokkos::Experimental::SYCLDeviceUSMSpace, HostSpace>(
      RecordBase::m_alloc_ptr, &header, sizeof(SharedAllocationHeader));
}

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

std::string SharedAllocationRecord<Kokkos::Experimental::SYCLDeviceUSMSpace,
                                   void>::get_label() const {
  SharedAllocationHeader header;

  Kokkos::Impl::DeepCopy<Kokkos::HostSpace,
                         Kokkos::Experimental::SYCLDeviceUSMSpace>(
      &header, RecordBase::head(), sizeof(SharedAllocationHeader));

  return std::string(header.m_label);
}

SharedAllocationRecord<Kokkos::Experimental::SYCLDeviceUSMSpace, void>*
SharedAllocationRecord<Kokkos::Experimental::SYCLDeviceUSMSpace, void>::
    allocate(const Kokkos::Experimental::SYCLDeviceUSMSpace& arg_space,
             const std::string& arg_label, const size_t arg_alloc_size) {
  return new SharedAllocationRecord(arg_space, arg_label, arg_alloc_size);
}

void SharedAllocationRecord<Kokkos::Experimental::SYCLDeviceUSMSpace, void>::
    deallocate(SharedAllocationRecord<void, void>* arg_rec) {
  delete static_cast<SharedAllocationRecord*>(arg_rec);
}

SharedAllocationRecord<Kokkos::Experimental::SYCLDeviceUSMSpace,
                       void>::~SharedAllocationRecord() {
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    SharedAllocationHeader header;
    Kokkos::Impl::DeepCopy<Kokkos::Experimental::SYCLDeviceUSMSpace,
                           Kokkos::HostSpace>(&header, RecordBase::m_alloc_ptr,
                                              sizeof(SharedAllocationHeader));

    Kokkos::Profiling::deallocateData(
        Kokkos::Profiling::make_space_handle(
            Kokkos::Experimental::SYCLDeviceUSMSpace::name()),
        header.m_label, data(), size());
  }

  m_space.deallocate(SharedAllocationRecord<void, void>::m_alloc_ptr,
                     SharedAllocationRecord<void, void>::m_alloc_size);
}

//----------------------------------------------------------------------------

void* SharedAllocationRecord<Kokkos::Experimental::SYCLDeviceUSMSpace, void>::
    allocate_tracked(const Kokkos::Experimental::SYCLDeviceUSMSpace& arg_space,
                     const std::string& arg_alloc_label,
                     const size_t arg_alloc_size) {
  if (!arg_alloc_size) return nullptr;

  SharedAllocationRecord* const r =
      allocate(arg_space, arg_alloc_label, arg_alloc_size);

  RecordBase::increment(r);

  return r->data();
}

void SharedAllocationRecord<Kokkos::Experimental::SYCLDeviceUSMSpace,
                            void>::deallocate_tracked(void* const
                                                          arg_alloc_ptr) {
  if (arg_alloc_ptr != nullptr) {
    SharedAllocationRecord* const r = get_record(arg_alloc_ptr);

    RecordBase::decrement(r);
  }
}

void* SharedAllocationRecord<Kokkos::Experimental::SYCLDeviceUSMSpace, void>::
    reallocate_tracked(void* const arg_alloc_ptr, const size_t arg_alloc_size) {
  SharedAllocationRecord* const r_old = get_record(arg_alloc_ptr);
  SharedAllocationRecord* const r_new =
      allocate(r_old->m_space, r_old->get_label(), arg_alloc_size);

  Kokkos::Impl::DeepCopy<Kokkos::Experimental::SYCLDeviceUSMSpace,
                         Kokkos::Experimental::SYCLDeviceUSMSpace>(
      r_new->data(), r_old->data(), std::min(r_old->size(), r_new->size()));

  RecordBase::increment(r_new);
  RecordBase::decrement(r_old);

  return r_new->data();
}

//----------------------------------------------------------------------------

SharedAllocationRecord<Kokkos::Experimental::SYCLDeviceUSMSpace, void>*
SharedAllocationRecord<Kokkos::Experimental::SYCLDeviceUSMSpace,
                       void>::get_record(void* alloc_ptr) {
  using Header = SharedAllocationHeader;
  using RecordSYCL =
      SharedAllocationRecord<Kokkos::Experimental::SYCLDeviceUSMSpace, void>;

  // Copy the header from the allocation
  Header head;

  Header const* const head_sycl =
      alloc_ptr ? Header::get_header(alloc_ptr) : nullptr;

  if (alloc_ptr) {
    Kokkos::Impl::DeepCopy<Kokkos::HostSpace,
                           Kokkos::Experimental::SYCLDeviceUSMSpace>(
        &head, head_sycl, sizeof(SharedAllocationHeader));
  }

  RecordSYCL* const record =
      alloc_ptr ? static_cast<RecordSYCL*>(head.m_record) : nullptr;

  if (!alloc_ptr || record->m_alloc_ptr != head_sycl) {
    Kokkos::Impl::throw_runtime_exception(
        std::string("Kokkos::Impl::SharedAllocationRecord< "
                    "Kokkos::Experimental::SYCLDeviceUSMSpace "
                    ", void >::get_record ERROR"));
  }

  return record;
}

// Iterate records to print orphaned memory ...
void SharedAllocationRecord<Kokkos::Experimental::SYCLDeviceUSMSpace, void>::
    print_records(std::ostream& s,
                  const Kokkos::Experimental::SYCLDeviceUSMSpace&,
                  bool detail) {
#ifdef KOKKOS_ENABLE_DEBUG
  SharedAllocationRecord<void, void>* r = &s_root_record;

  char buffer[256];

  SharedAllocationHeader head;

  if (detail) {
    do {
      if (r->m_alloc_ptr) {
        Kokkos::Impl::DeepCopy<Kokkos::HostSpace,
                               Kokkos::Experimental::SYCLDeviceUSMSpace>(
            &head, r->m_alloc_ptr, sizeof(SharedAllocationHeader));
      } else {
        head.m_label[0] = 0;
      }

      // Formatting dependent on sizeof(uintptr_t)
      const char* format_string;

      if (sizeof(uintptr_t) == sizeof(unsigned long)) {
        format_string =
            "SYCL addr( 0x%.12lx ) list( 0x%.12lx 0x%.12lx ) extent[ 0x%.12lx "
            "+ %.8ld ] count(%d) dealloc(0x%.12lx) %s\n";
      } else if (sizeof(uintptr_t) == sizeof(unsigned long long)) {
        format_string =
            "SYCL addr( 0x%.12llx ) list( 0x%.12llx 0x%.12llx ) extent[ "
            "0x%.12llx + %.8ld ] count(%d) dealloc(0x%.12llx) %s\n";
      }

      snprintf(buffer, 256, format_string, reinterpret_cast<uintptr_t>(r),
               reinterpret_cast<uintptr_t>(r->m_prev),
               reinterpret_cast<uintptr_t>(r->m_next),
               reinterpret_cast<uintptr_t>(r->m_alloc_ptr), r->m_alloc_size,
               r->m_count, reinterpret_cast<uintptr_t>(r->m_dealloc),
               head.m_label);
      s << buffer;
      r = r->m_next;
    } while (r != &s_root_record);
  } else {
    do {
      if (r->m_alloc_ptr) {
        Kokkos::Impl::DeepCopy<Kokkos::HostSpace,
                               Kokkos::Experimental::SYCLDeviceUSMSpace>(
            &head, r->m_alloc_ptr, sizeof(SharedAllocationHeader));

        // Formatting dependent on sizeof(uintptr_t)
        const char* format_string;

        if (sizeof(uintptr_t) == sizeof(unsigned long)) {
          format_string = "SYCL [ 0x%.12lx + %ld ] %s\n";
        } else if (sizeof(uintptr_t) == sizeof(unsigned long long)) {
          format_string = "SYCL [ 0x%.12llx + %ld ] %s\n";
        }

        snprintf(buffer, 256, format_string,
                 reinterpret_cast<uintptr_t>(r->data()), r->size(),
                 head.m_label);
      } else {
        snprintf(buffer, 256, "SYCL [ 0 + 0 ]\n");
      }
      s << buffer;
      r = r->m_next;
    } while (r != &s_root_record);
  }
#else
  (void)s;
  (void)detail;
  throw_runtime_exception(
      "Kokkos::Impl::SharedAllocationRecord<SYCLDeviceUSMSpace>::print_records"
      " only works with KOKKOS_ENABLE_DEBUG enabled");
#endif
}

}  // namespace Impl
}  // namespace Kokkos
