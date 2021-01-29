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

#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <atomic>
#include <Kokkos_Macros.hpp>

#include <Kokkos_Core.hpp>
#include <Kokkos_HIP.hpp>
#include <Kokkos_HIP_Space.hpp>

#include <impl/Kokkos_Error.hpp>
#include <impl/Kokkos_MemorySpace.hpp>

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
namespace Kokkos {
namespace Impl {

namespace {
hipStream_t get_deep_copy_stream() {
  static hipStream_t s = nullptr;
  if (s == nullptr) {
    HIP_SAFE_CALL(hipStreamCreate(&s));
  }
  return s;
}
}  // namespace

DeepCopy<Kokkos::Experimental::HIPSpace, Kokkos::Experimental::HIPSpace,
         Kokkos::Experimental::HIP>::DeepCopy(void* dst, const void* src,
                                              size_t n) {
  HIP_SAFE_CALL(hipMemcpy(dst, src, n, hipMemcpyDefault));
}

DeepCopy<HostSpace, Kokkos::Experimental::HIPSpace,
         Kokkos::Experimental::HIP>::DeepCopy(void* dst, const void* src,
                                              size_t n) {
  HIP_SAFE_CALL(hipMemcpy(dst, src, n, hipMemcpyDefault));
}

DeepCopy<Kokkos::Experimental::HIPSpace, HostSpace,
         Kokkos::Experimental::HIP>::DeepCopy(void* dst, const void* src,
                                              size_t n) {
  HIP_SAFE_CALL(hipMemcpy(dst, src, n, hipMemcpyDefault));
}

DeepCopy<Kokkos::Experimental::HIPSpace, Kokkos::Experimental::HIPSpace,
         Kokkos::Experimental::HIP>::DeepCopy(const Kokkos::Experimental::HIP&
                                                  instance,
                                              void* dst, const void* src,
                                              size_t n) {
  HIP_SAFE_CALL(
      hipMemcpyAsync(dst, src, n, hipMemcpyDefault, instance.hip_stream()));
}

DeepCopy<HostSpace, Kokkos::Experimental::HIPSpace, Kokkos::Experimental::HIP>::
    DeepCopy(const Kokkos::Experimental::HIP& instance, void* dst,
             const void* src, size_t n) {
  HIP_SAFE_CALL(
      hipMemcpyAsync(dst, src, n, hipMemcpyDefault, instance.hip_stream()));
}

DeepCopy<Kokkos::Experimental::HIPSpace, HostSpace, Kokkos::Experimental::HIP>::
    DeepCopy(const Kokkos::Experimental::HIP& instance, void* dst,
             const void* src, size_t n) {
  HIP_SAFE_CALL(
      hipMemcpyAsync(dst, src, n, hipMemcpyDefault, instance.hip_stream()));
}

DeepCopy<Kokkos::Experimental::HIPHostPinnedSpace,
         Kokkos::Experimental::HIPHostPinnedSpace,
         Kokkos::Experimental::HIP>::DeepCopy(void* dst, const void* src,
                                              size_t n) {
  HIP_SAFE_CALL(hipMemcpy(dst, src, n, hipMemcpyDefault));
}

DeepCopy<HostSpace, Kokkos::Experimental::HIPHostPinnedSpace,
         Kokkos::Experimental::HIP>::DeepCopy(void* dst, const void* src,
                                              size_t n) {
  HIP_SAFE_CALL(hipMemcpy(dst, src, n, hipMemcpyDefault));
}

DeepCopy<Kokkos::Experimental::HIPHostPinnedSpace, HostSpace,
         Kokkos::Experimental::HIP>::DeepCopy(void* dst, const void* src,
                                              size_t n) {
  HIP_SAFE_CALL(hipMemcpy(dst, src, n, hipMemcpyDefault));
}

DeepCopy<Kokkos::Experimental::HIPHostPinnedSpace,
         Kokkos::Experimental::HIPHostPinnedSpace, Kokkos::Experimental::HIP>::
    DeepCopy(const Kokkos::Experimental::HIP& instance, void* dst,
             const void* src, size_t n) {
  HIP_SAFE_CALL(
      hipMemcpyAsync(dst, src, n, hipMemcpyDefault, instance.hip_stream()));
}

DeepCopy<HostSpace, Kokkos::Experimental::HIPHostPinnedSpace,
         Kokkos::Experimental::HIP>::DeepCopy(const Kokkos::Experimental::HIP&
                                                  instance,
                                              void* dst, const void* src,
                                              size_t n) {
  HIP_SAFE_CALL(
      hipMemcpyAsync(dst, src, n, hipMemcpyDefault, instance.hip_stream()));
}

DeepCopy<Kokkos::Experimental::HIPHostPinnedSpace, HostSpace,
         Kokkos::Experimental::HIP>::DeepCopy(const Kokkos::Experimental::HIP&
                                                  instance,
                                              void* dst, const void* src,
                                              size_t n) {
  HIP_SAFE_CALL(
      hipMemcpyAsync(dst, src, n, hipMemcpyDefault, instance.hip_stream()));
}

void DeepCopyAsyncHIP(void* dst, void const* src, size_t n) {
  hipStream_t s = get_deep_copy_stream();
  HIP_SAFE_CALL(hipMemcpyAsync(dst, src, n, hipMemcpyDefault, s));
  HIP_SAFE_CALL(hipStreamSynchronize(s));
}

}  // namespace Impl
}  // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {

void Experimental::HIPSpace::access_error() {
  const std::string msg(
      "Kokkos::Experimental::HIPSpace::access_error attempt to execute "
      "Experimental::HIP function from non-HIP space");
  Kokkos::Impl::throw_runtime_exception(msg);
}

void Experimental::HIPSpace::access_error(const void* const) {
  const std::string msg(
      "Kokkos::Experimental::HIPSpace::access_error attempt to execute "
      "Experimental::HIP function from non-HIP space");
  Kokkos::Impl::throw_runtime_exception(msg);
}

}  // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Experimental {

HIPSpace::HIPSpace() : m_device(HIP().hip_device()) {}

HIPHostPinnedSpace::HIPHostPinnedSpace() {}

void* HIPSpace::allocate(const size_t arg_alloc_size) const {
  return allocate("[unlabeled]", arg_alloc_size);
}
void* HIPSpace::allocate(

    const char* arg_label, const size_t arg_alloc_size,
    const size_t arg_logical_size) const {
  return impl_allocate(arg_label, arg_alloc_size, arg_logical_size);
}
void* HIPSpace::impl_allocate(

    const char* arg_label, const size_t arg_alloc_size,
    const size_t arg_logical_size,
    const Kokkos::Tools::SpaceHandle arg_handle) const {
  void* ptr = nullptr;

  auto const error_code = hipMalloc(&ptr, arg_alloc_size);
  if (error_code != hipSuccess) {
    // This is the only way to clear the last error, which we should do here
    // since we're turning it into an exception here
    (void)hipGetLastError();
    throw HIPRawMemoryAllocationFailure(
        arg_alloc_size, error_code,
        RawMemoryAllocationFailure::AllocationMechanism::HIPMalloc);
  }
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    const size_t reported_size =
        (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;
    Kokkos::Profiling::allocateData(arg_handle, arg_label, ptr, reported_size);
  }

  return ptr;
}

void* HIPHostPinnedSpace::allocate(const size_t arg_alloc_size) const {
  return allocate("[unlabeled]", arg_alloc_size);
}
void* HIPHostPinnedSpace::allocate(const char* arg_label,
                                   const size_t arg_alloc_size,
                                   const size_t arg_logical_size) const {
  return impl_allocate(arg_label, arg_alloc_size, arg_logical_size);
}
void* HIPHostPinnedSpace::impl_allocate(
    const char* arg_label, const size_t arg_alloc_size,
    const size_t arg_logical_size,
    const Kokkos::Tools::SpaceHandle arg_handle) const {
  void* ptr = nullptr;

  auto const error_code = hipHostMalloc(&ptr, arg_alloc_size);
  if (error_code != hipSuccess) {
    // This is the only way to clear the last error, which we should do here
    // since we're turning it into an exception here
    (void)hipGetLastError();
    throw HIPRawMemoryAllocationFailure(
        arg_alloc_size, error_code,
        RawMemoryAllocationFailure::AllocationMechanism::HIPHostMalloc);
  }
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    const size_t reported_size =
        (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;
    Kokkos::Profiling::allocateData(arg_handle, arg_label, ptr, reported_size);
  }

  return ptr;
}
void HIPSpace::deallocate(void* const arg_alloc_ptr,
                          const size_t arg_alloc_size) const {
  deallocate("[unlabeled]", arg_alloc_ptr, arg_alloc_size);
}
void HIPSpace::deallocate(const char* arg_label, void* const arg_alloc_ptr,
                          const size_t arg_alloc_size,
                          const size_t arg_logical_size) const {
  impl_deallocate(arg_label, arg_alloc_ptr, arg_alloc_size, arg_logical_size);
}
void HIPSpace::impl_deallocate(
    const char* arg_label, void* const arg_alloc_ptr,
    const size_t arg_alloc_size, const size_t arg_logical_size,
    const Kokkos::Tools::SpaceHandle arg_handle) const {
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    const size_t reported_size =
        (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;
    Kokkos::Profiling::deallocateData(arg_handle, arg_label, arg_alloc_ptr,
                                      reported_size);
  }
  HIP_SAFE_CALL(hipFree(arg_alloc_ptr));
}

void HIPHostPinnedSpace::deallocate(void* const arg_alloc_ptr,
                                    const size_t arg_alloc_size) const {
  deallocate("[unlabeled]", arg_alloc_ptr, arg_alloc_size);
}

void HIPHostPinnedSpace::deallocate(const char* arg_label,
                                    void* const arg_alloc_ptr,
                                    const size_t arg_alloc_size,
                                    const size_t arg_logical_size) const {
  impl_deallocate(arg_label, arg_alloc_ptr, arg_alloc_size, arg_logical_size);
}
void HIPHostPinnedSpace::impl_deallocate(
    const char* arg_label, void* const arg_alloc_ptr,
    const size_t arg_alloc_size, const size_t arg_logical_size,
    const Kokkos::Tools::SpaceHandle arg_handle) const {
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    const size_t reported_size =
        (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;
    Kokkos::Profiling::deallocateData(arg_handle, arg_label, arg_alloc_ptr,
                                      reported_size);
  }
  HIP_SAFE_CALL(hipHostFree(arg_alloc_ptr));
}

}  // namespace Experimental
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

#ifdef KOKKOS_ENABLE_DEBUG
SharedAllocationRecord<void, void>
    SharedAllocationRecord<Kokkos::Experimental::HIPSpace, void>::s_root_record;

SharedAllocationRecord<void, void> SharedAllocationRecord<
    Kokkos::Experimental::HIPHostPinnedSpace, void>::s_root_record;
#endif

std::string SharedAllocationRecord<Kokkos::Experimental::HIPSpace,
                                   void>::get_label() const {
  SharedAllocationHeader header;

  Kokkos::Impl::DeepCopy<Kokkos::HostSpace, Kokkos::Experimental::HIPSpace>(
      &header, RecordBase::head(), sizeof(SharedAllocationHeader));

  return std::string(header.m_label);
}

std::string SharedAllocationRecord<Kokkos::Experimental::HIPHostPinnedSpace,
                                   void>::get_label() const {
  return std::string(RecordBase::head()->m_label);
}

SharedAllocationRecord<Kokkos::Experimental::HIPSpace, void>*
SharedAllocationRecord<Kokkos::Experimental::HIPSpace, void>::allocate(
    const Kokkos::Experimental::HIPSpace& arg_space,
    const std::string& arg_label, const size_t arg_alloc_size) {
  return new SharedAllocationRecord(arg_space, arg_label, arg_alloc_size);
}

SharedAllocationRecord<Kokkos::Experimental::HIPHostPinnedSpace, void>*
SharedAllocationRecord<Kokkos::Experimental::HIPHostPinnedSpace, void>::
    allocate(const Kokkos::Experimental::HIPHostPinnedSpace& arg_space,
             const std::string& arg_label, const size_t arg_alloc_size) {
  return new SharedAllocationRecord(arg_space, arg_label, arg_alloc_size);
}

void SharedAllocationRecord<Kokkos::Experimental::HIPSpace, void>::deallocate(
    SharedAllocationRecord<void, void>* arg_rec) {
  delete static_cast<SharedAllocationRecord*>(arg_rec);
}

void SharedAllocationRecord<Kokkos::Experimental::HIPHostPinnedSpace, void>::
    deallocate(SharedAllocationRecord<void, void>* arg_rec) {
  delete static_cast<SharedAllocationRecord*>(arg_rec);
}

SharedAllocationRecord<Kokkos::Experimental::HIPSpace,
                       void>::~SharedAllocationRecord() {
  const char* label = nullptr;
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    SharedAllocationHeader header;
    Kokkos::Impl::DeepCopy<Kokkos::Experimental::HIPSpace, HostSpace>(
        &header, RecordBase::m_alloc_ptr, sizeof(SharedAllocationHeader));
    label = header.label();
  }
  auto alloc_size = SharedAllocationRecord<void, void>::m_alloc_size;
  m_space.deallocate(label, SharedAllocationRecord<void, void>::m_alloc_ptr,
                     alloc_size, (alloc_size - sizeof(SharedAllocationHeader)));
}

SharedAllocationRecord<Kokkos::Experimental::HIPHostPinnedSpace,
                       void>::~SharedAllocationRecord() {
  m_space.deallocate(RecordBase::m_alloc_ptr->m_label,
                     SharedAllocationRecord<void, void>::m_alloc_ptr,
                     SharedAllocationRecord<void, void>::m_alloc_size);
}

SharedAllocationRecord<Kokkos::Experimental::HIPSpace, void>::
    SharedAllocationRecord(
        const Kokkos::Experimental::HIPSpace& arg_space,
        const std::string& arg_label, const size_t arg_alloc_size,
        const SharedAllocationRecord<void, void>::function_type arg_dealloc)
    // Pass through allocated [ SharedAllocationHeader , user_memory ]
    // Pass through deallocation function
    : SharedAllocationRecord<void, void>(
#ifdef KOKKOS_ENABLE_DEBUG
          &SharedAllocationRecord<Kokkos::Experimental::HIPSpace,
                                  void>::s_root_record,
#endif
          Kokkos::Impl::checked_allocation_with_header(arg_space, arg_label,
                                                       arg_alloc_size),
          sizeof(SharedAllocationHeader) + arg_alloc_size, arg_dealloc),
      m_space(arg_space) {

  SharedAllocationHeader header;

  // Fill in the Header information
  header.m_record = static_cast<SharedAllocationRecord<void, void>*>(this);

  strncpy(header.m_label, arg_label.c_str(),
          SharedAllocationHeader::maximum_label_length);
  // Set last element zero, in case c_str is too long
  header.m_label[SharedAllocationHeader::maximum_label_length - 1] = (char)0;

  // Copy to device memory
  Kokkos::Impl::DeepCopy<Kokkos::Experimental::HIPSpace, HostSpace>(
      RecordBase::m_alloc_ptr, &header, sizeof(SharedAllocationHeader));
}

SharedAllocationRecord<Kokkos::Experimental::HIPHostPinnedSpace, void>::
    SharedAllocationRecord(
        const Kokkos::Experimental::HIPHostPinnedSpace& arg_space,
        const std::string& arg_label, const size_t arg_alloc_size,
        const SharedAllocationRecord<void, void>::function_type arg_dealloc)
    // Pass through allocated [ SharedAllocationHeader , user_memory ]
    // Pass through deallocation function
    : SharedAllocationRecord<void, void>(
#ifdef KOKKOS_ENABLE_DEBUG
          &SharedAllocationRecord<Kokkos::Experimental::HIPHostPinnedSpace,
                                  void>::s_root_record,
#endif
          Kokkos::Impl::checked_allocation_with_header(arg_space, arg_label,
                                                       arg_alloc_size),
          sizeof(SharedAllocationHeader) + arg_alloc_size, arg_dealloc),
      m_space(arg_space) {
  // Fill in the Header information, directly accessible via host pinned memory

  RecordBase::m_alloc_ptr->m_record = this;

  strncpy(RecordBase::m_alloc_ptr->m_label, arg_label.c_str(),
          SharedAllocationHeader::maximum_label_length);
  // Set last element zero, in case c_str is too long
  RecordBase::m_alloc_ptr
      ->m_label[SharedAllocationHeader::maximum_label_length - 1] = (char)0;
}

//----------------------------------------------------------------------------

void* SharedAllocationRecord<Kokkos::Experimental::HIPSpace, void>::
    allocate_tracked(const Kokkos::Experimental::HIPSpace& arg_space,
                     const std::string& arg_alloc_label,
                     const size_t arg_alloc_size) {
  if (!arg_alloc_size) return nullptr;

  SharedAllocationRecord* const r =
      allocate(arg_space, arg_alloc_label, arg_alloc_size);

  RecordBase::increment(r);

  return r->data();
}

void SharedAllocationRecord<Kokkos::Experimental::HIPSpace,
                            void>::deallocate_tracked(void* const
                                                          arg_alloc_ptr) {
  if (arg_alloc_ptr != nullptr) {
    SharedAllocationRecord* const r = get_record(arg_alloc_ptr);

    RecordBase::decrement(r);
  }
}

void* SharedAllocationRecord<Kokkos::Experimental::HIPSpace, void>::
    reallocate_tracked(void* const arg_alloc_ptr, const size_t arg_alloc_size) {
  SharedAllocationRecord* const r_old = get_record(arg_alloc_ptr);
  SharedAllocationRecord* const r_new =
      allocate(r_old->m_space, r_old->get_label(), arg_alloc_size);

  Kokkos::Impl::DeepCopy<Kokkos::Experimental::HIPSpace,
                         Kokkos::Experimental::HIPSpace>(
      r_new->data(), r_old->data(), std::min(r_old->size(), r_new->size()));

  RecordBase::increment(r_new);
  RecordBase::decrement(r_old);

  return r_new->data();
}

void* SharedAllocationRecord<Kokkos::Experimental::HIPHostPinnedSpace, void>::
    allocate_tracked(const Kokkos::Experimental::HIPHostPinnedSpace& arg_space,
                     const std::string& arg_alloc_label,
                     const size_t arg_alloc_size) {
  if (!arg_alloc_size) return nullptr;

  SharedAllocationRecord* const r =
      allocate(arg_space, arg_alloc_label, arg_alloc_size);

  RecordBase::increment(r);

  return r->data();
}

void SharedAllocationRecord<Kokkos::Experimental::HIPHostPinnedSpace,
                            void>::deallocate_tracked(void* const
                                                          arg_alloc_ptr) {
  if (arg_alloc_ptr) {
    SharedAllocationRecord* const r = get_record(arg_alloc_ptr);

    RecordBase::decrement(r);
  }
}

void* SharedAllocationRecord<Kokkos::Experimental::HIPHostPinnedSpace, void>::
    reallocate_tracked(void* const arg_alloc_ptr, const size_t arg_alloc_size) {
  SharedAllocationRecord* const r_old = get_record(arg_alloc_ptr);
  SharedAllocationRecord* const r_new =
      allocate(r_old->m_space, r_old->get_label(), arg_alloc_size);

  using HIPHostPinnedSpace = Kokkos::Experimental::HIPHostPinnedSpace;
  Kokkos::Impl::DeepCopy<HIPHostPinnedSpace, HIPHostPinnedSpace>(
      r_new->data(), r_old->data(), std::min(r_old->size(), r_new->size()));

  RecordBase::increment(r_new);
  RecordBase::decrement(r_old);

  return r_new->data();
}

//----------------------------------------------------------------------------

SharedAllocationRecord<Kokkos::Experimental::HIPSpace, void>*
SharedAllocationRecord<Kokkos::Experimental::HIPSpace, void>::get_record(
    void* alloc_ptr) {
  using Header = SharedAllocationHeader;
  using RecordHIP =
      SharedAllocationRecord<Kokkos::Experimental::HIPSpace, void>;

  // Copy the header from the allocation
  Header head;

  Header const* const head_hip =
      alloc_ptr ? Header::get_header(alloc_ptr) : nullptr;

  if (alloc_ptr) {
    Kokkos::Impl::DeepCopy<HostSpace, Kokkos::Experimental::HIPSpace>(
        &head, head_hip, sizeof(SharedAllocationHeader));
  }

  RecordHIP* const record =
      alloc_ptr ? static_cast<RecordHIP*>(head.m_record) : nullptr;

  if (!alloc_ptr || record->m_alloc_ptr != head_hip) {
    Kokkos::Impl::throw_runtime_exception(std::string(
        "Kokkos::Impl::SharedAllocationRecord< Kokkos::Experimental::HIPSpace "
        ", void >::get_record ERROR"));
  }

  return record;
}

SharedAllocationRecord<Kokkos::Experimental::HIPHostPinnedSpace, void>*
SharedAllocationRecord<Kokkos::Experimental::HIPHostPinnedSpace,
                       void>::get_record(void* alloc_ptr) {
  using Header = SharedAllocationHeader;
  using RecordHIP =
      SharedAllocationRecord<Kokkos::Experimental::HIPHostPinnedSpace, void>;

  Header* const h =
      alloc_ptr ? reinterpret_cast<Header*>(alloc_ptr) - 1 : nullptr;

  if (!alloc_ptr || h->m_record->m_alloc_ptr != h) {
    Kokkos::Impl::throw_runtime_exception(std::string(
        "Kokkos::Impl::SharedAllocationRecord< "
        "Kokkos::Experimental::HIPHostPinnedSpace , void >::get_record ERROR"));
  }

  return static_cast<RecordHIP*>(h->m_record);
}

// Iterate records to print orphaned memory ...
void SharedAllocationRecord<Kokkos::Experimental::HIPSpace, void>::
    print_records(std::ostream& s, const Kokkos::Experimental::HIPSpace&,
                  bool detail) {
#ifdef KOKKOS_ENABLE_DEBUG
  SharedAllocationRecord<void, void>* r = &s_root_record;

  char buffer[256];

  SharedAllocationHeader head;

  if (detail) {
    do {
      if (r->m_alloc_ptr) {
        Kokkos::Impl::DeepCopy<HostSpace, Kokkos::Experimental::HIPSpace>(
            &head, r->m_alloc_ptr, sizeof(SharedAllocationHeader));
      } else {
        head.m_label[0] = 0;
      }

      // Formatting dependent on sizeof(uintptr_t)
      const char* format_string;

      if (sizeof(uintptr_t) == sizeof(unsigned long)) {
        format_string =
            "HIP addr( 0x%.12lx ) list( 0x%.12lx 0x%.12lx ) extent[ 0x%.12lx + "
            "%.8ld ] count(%d) dealloc(0x%.12lx) %s\n";
      } else if (sizeof(uintptr_t) == sizeof(unsigned long long)) {
        format_string =
            "HIP addr( 0x%.12llx ) list( 0x%.12llx 0x%.12llx ) extent[ "
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
        Kokkos::Impl::DeepCopy<HostSpace, Kokkos::Experimental::HIPSpace>(
            &head, r->m_alloc_ptr, sizeof(SharedAllocationHeader));

        // Formatting dependent on sizeof(uintptr_t)
        const char* format_string;

        if (sizeof(uintptr_t) == sizeof(unsigned long)) {
          format_string = "HIP [ 0x%.12lx + %ld ] %s\n";
        } else if (sizeof(uintptr_t) == sizeof(unsigned long long)) {
          format_string = "HIP [ 0x%.12llx + %ld ] %s\n";
        }

        snprintf(buffer, 256, format_string,
                 reinterpret_cast<uintptr_t>(r->data()), r->size(),
                 head.m_label);
      } else {
        snprintf(buffer, 256, "HIP [ 0 + 0 ]\n");
      }
      s << buffer;
      r = r->m_next;
    } while (r != &s_root_record);
  }
#else
  (void)s;
  (void)detail;
  throw_runtime_exception(
      "Kokkos::Impl::SharedAllocationRecord<HIPSpace>::print_records"
      " only works with KOKKOS_ENABLE_DEBUG enabled");
#endif
}

}  // namespace Impl
}  // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
namespace Kokkos {
namespace Impl {
int get_gpu(const InitArguments& args);
}
namespace Experimental {

int HIP::concurrency() {
  auto const& prop = hip_device_prop();
  return prop.maxThreadsPerMultiProcessor * prop.multiProcessorCount;
}
int HIP::impl_is_initialized() {
  return Impl::HIPInternal::singleton().is_initialized();
}

void HIP::impl_initialize(const HIP::SelectDevice config) {
  Impl::HIPInternal::singleton().initialize(config.hip_device_id);
}

void HIP::impl_finalize() { Impl::HIPInternal::singleton().finalize(); }

HIP::HIP()
    : m_space_instance(&Impl::HIPInternal::singleton()), m_counter(nullptr) {
  Impl::HIPInternal::singleton().verify_is_initialized(
      "HIP instance constructor");
}

HIP::HIP(hipStream_t const stream)
    : m_space_instance(new Impl::HIPInternal), m_counter(new int(1)) {
  Impl::HIPInternal::singleton().verify_is_initialized(
      "HIP instance constructor");
  m_space_instance->initialize(Impl::HIPInternal::singleton().m_hipDev, stream);
}

KOKKOS_FUNCTION HIP::HIP(HIP&& other) noexcept {
  m_space_instance       = other.m_space_instance;
  other.m_space_instance = nullptr;
  m_counter              = other.m_counter;
  other.m_counter        = nullptr;
}

KOKKOS_FUNCTION HIP::HIP(HIP const& other)
    : m_space_instance(other.m_space_instance), m_counter(other.m_counter) {
#ifndef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HIP_GPU
  if (m_counter) Kokkos::atomic_add(m_counter, 1);
#endif
}

KOKKOS_FUNCTION HIP& HIP::operator=(HIP&& other) noexcept {
  m_space_instance       = other.m_space_instance;
  other.m_space_instance = nullptr;
  m_counter              = other.m_counter;
  other.m_counter        = nullptr;

  return *this;
}

KOKKOS_FUNCTION HIP& HIP::operator=(HIP const& other) {
  m_space_instance = other.m_space_instance;
  m_counter        = other.m_counter;
#ifndef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HIP_GPU
  if (m_counter) Kokkos::atomic_add(m_counter, 1);
#endif

  return *this;
}

KOKKOS_FUNCTION HIP::~HIP() noexcept {
#ifndef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HIP_GPU
  if (m_counter == nullptr) return;
  int const count = Kokkos::atomic_fetch_sub(m_counter, 1);
  if (count == 1) {
    delete m_counter;
    m_space_instance->finalize();
    delete m_space_instance;
  }
#endif
}

void HIP::print_configuration(std::ostream& s, const bool) {
  Impl::HIPInternal::singleton().print_configuration(s);
}

void HIP::impl_static_fence() { HIP_SAFE_CALL(hipDeviceSynchronize()); }

void HIP::fence() const { m_space_instance->fence(); }

hipStream_t HIP::hip_stream() const { return m_space_instance->m_stream; }

int HIP::hip_device() const { return impl_internal_space_instance()->m_hipDev; }

hipDeviceProp_t const& HIP::hip_device_prop() {
  return Impl::HIPInternal::singleton().m_deviceProp;
}

const char* HIP::name() { return "HIP"; }

}  // namespace Experimental

namespace Impl {

int g_hip_space_factory_initialized =
    initialize_space_factory<HIPSpaceInitializer>("150_HIP");

void HIPSpaceInitializer::initialize(const InitArguments& args) {
  int use_gpu = Impl::get_gpu(args);

  if (std::is_same<Kokkos::Experimental::HIP,
                   Kokkos::DefaultExecutionSpace>::value ||
      0 < use_gpu) {
    if (use_gpu > -1) {
      Kokkos::Experimental::HIP::impl_initialize(
          Kokkos::Experimental::HIP::SelectDevice(use_gpu));
    } else {
      Kokkos::Experimental::HIP::impl_initialize();
    }
  }
}

void HIPSpaceInitializer::finalize(const bool all_spaces) {
  if (std::is_same<Kokkos::Experimental::HIP,
                   Kokkos::DefaultExecutionSpace>::value ||
      all_spaces) {
    if (Kokkos::Experimental::HIP::impl_is_initialized())
      Kokkos::Experimental::HIP::impl_finalize();
  }
}

void HIPSpaceInitializer::fence() {
  Kokkos::Experimental::HIP::impl_static_fence();
}

void HIPSpaceInitializer::print_configuration(std::ostream& msg,
                                              const bool detail) {
  msg << "Devices:" << std::endl;
  msg << "  KOKKOS_ENABLE_HIP: ";
  msg << "yes" << std::endl;

  msg << "HIP Options:" << std::endl;
  msg << "  KOKKOS_ENABLE_HIP_RELOCATABLE_DEVICE_CODE: ";
#ifdef KOKKOS_ENABLE_HIP_RELOCATABLE_DEVICE_CODE
  msg << "yes" << std::endl;
#else
  msg << "no" << std::endl;
#endif

  msg << "\nRuntime Configuration:" << std::endl;
  Experimental::HIP::print_configuration(msg, detail);
}

}  // namespace Impl
}  // namespace Kokkos
