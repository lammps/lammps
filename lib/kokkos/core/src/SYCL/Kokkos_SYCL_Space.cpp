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

#include <Kokkos_Macros.hpp>

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
auto USM_memcpy(sycl::queue& q, void* dst, const void* src, size_t n) {
  return q.memcpy(dst, src, n);
}

void USM_memcpy(Kokkos::Experimental::Impl::SYCLInternal& space, void* dst,
                const void* src, size_t n) {
  (void)USM_memcpy(*space.m_queue, dst, src, n);
}

void USM_memcpy(void* dst, const void* src, size_t n) {
  Experimental::SYCL().fence();
  auto event = USM_memcpy(
      *Experimental::Impl::SYCLInternal::singleton().m_queue, dst, src, n);
  Experimental::Impl::SYCLInternal::fence(event);
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

SYCLDeviceUSMSpace::SYCLDeviceUSMSpace()
    : m_queue(*SYCL().impl_internal_space_instance()->m_queue) {}
SYCLDeviceUSMSpace::SYCLDeviceUSMSpace(sycl::queue queue)
    : m_queue(std::move(queue)) {}

SYCLSharedUSMSpace::SYCLSharedUSMSpace()
    : m_queue(*SYCL().impl_internal_space_instance()->m_queue) {}
SYCLSharedUSMSpace::SYCLSharedUSMSpace(sycl::queue queue)
    : m_queue(std::move(queue)) {}

void* allocate_sycl(
    const char* arg_label, const size_t arg_alloc_size,
    const size_t arg_logical_size, const Kokkos::Tools::SpaceHandle arg_handle,
    const RawMemoryAllocationFailure::AllocationMechanism failure_tag,
    const sycl::usm::alloc allocation_kind, const sycl::queue& queue) {
  void* const hostPtr = sycl::malloc(arg_alloc_size, queue, allocation_kind);

  if (hostPtr == nullptr)
    throw RawMemoryAllocationFailure(
        arg_alloc_size, 1, RawMemoryAllocationFailure::FailureMode::Unknown,
        failure_tag);

  if (Kokkos::Profiling::profileLibraryLoaded()) {
    const size_t reported_size =
        (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;
    Kokkos::Profiling::allocateData(arg_handle, arg_label, hostPtr,
                                    reported_size);
  }

  return hostPtr;
}

void* SYCLDeviceUSMSpace::allocate(const size_t arg_alloc_size) const {
  return allocate("[unlabeled]", arg_alloc_size);
}

void* SYCLDeviceUSMSpace::allocate(const char* arg_label,
                                   const size_t arg_alloc_size,
                                   const size_t arg_logical_size) const {
  return allocate_sycl(
      arg_label, arg_alloc_size, arg_logical_size,
      Kokkos::Tools::make_space_handle(name()),
      RawMemoryAllocationFailure::AllocationMechanism::SYCLMallocDevice,
      sycl::usm::alloc::device, m_queue);
}

void* SYCLSharedUSMSpace::allocate(const size_t arg_alloc_size) const {
  return allocate("[unlabeled]", arg_alloc_size);
}
void* SYCLSharedUSMSpace::allocate(const char* arg_label,
                                   const size_t arg_alloc_size,
                                   const size_t arg_logical_size) const {
  return allocate_sycl(
      arg_label, arg_alloc_size, arg_logical_size,
      Kokkos::Tools::make_space_handle(name()),
      RawMemoryAllocationFailure::AllocationMechanism::SYCLMallocShared,
      sycl::usm::alloc::shared, m_queue);
}

void sycl_deallocate(const char* arg_label, void* const arg_alloc_ptr,
                     const size_t arg_alloc_size, const size_t arg_logical_size,
                     const Kokkos::Tools::SpaceHandle arg_handle,
                     const sycl::queue& queue) {
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    const size_t reported_size =
        (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;
    Kokkos::Profiling::deallocateData(arg_handle, arg_label, arg_alloc_ptr,
                                      reported_size);
  }

  sycl::free(arg_alloc_ptr, queue);
}

void SYCLDeviceUSMSpace::deallocate(void* const arg_alloc_ptr,
                                    const size_t arg_alloc_size) const {
  deallocate("[unlabeled]", arg_alloc_ptr, arg_alloc_size);
}
void SYCLDeviceUSMSpace::deallocate(const char* arg_label,
                                    void* const arg_alloc_ptr,
                                    const size_t arg_alloc_size,
                                    const size_t arg_logical_size) const {
  sycl_deallocate(arg_label, arg_alloc_ptr, arg_alloc_size, arg_logical_size,
                  Kokkos::Tools::make_space_handle(name()), m_queue);
}

void SYCLSharedUSMSpace::deallocate(void* const arg_alloc_ptr,
                                    const size_t arg_alloc_size) const {
  deallocate("[unlabeled]", arg_alloc_ptr, arg_alloc_size);
}

void SYCLSharedUSMSpace::deallocate(const char* arg_label,
                                    void* const arg_alloc_ptr,
                                    const size_t arg_alloc_size,
                                    const size_t arg_logical_size) const {
  sycl_deallocate(arg_label, arg_alloc_ptr, arg_alloc_size, arg_logical_size,
                  Kokkos::Tools::make_space_handle(name()), m_queue);
}

}  // namespace Experimental
}  // namespace Kokkos

namespace Kokkos {
namespace Impl {

#ifdef KOKKOS_ENABLE_DEBUG
SharedAllocationRecord<void, void> SharedAllocationRecord<
    Kokkos::Experimental::SYCLDeviceUSMSpace, void>::s_root_record;

SharedAllocationRecord<void, void> SharedAllocationRecord<
    Kokkos::Experimental::SYCLSharedUSMSpace, void>::s_root_record;
#endif

SharedAllocationRecord<Kokkos::Experimental::SYCLDeviceUSMSpace, void>::
    SharedAllocationRecord(
        const Kokkos::Experimental::SYCLDeviceUSMSpace& space,
        const std::string& label, const size_t size,
        const SharedAllocationRecord<void, void>::function_type dealloc)
    // Pass through allocated [ SharedAllocationHeader , user_memory ]
    // Pass through deallocation function
    : base_t(
#ifdef KOKKOS_ENABLE_DEBUG
          &SharedAllocationRecord<Kokkos::Experimental::SYCLDeviceUSMSpace,
                                  void>::s_root_record,
#endif
          Kokkos::Impl::checked_allocation_with_header(space, label, size),
          sizeof(SharedAllocationHeader) + size, dealloc),
      m_space(space) {
  SharedAllocationHeader header;

  this->base_t::_fill_host_accessible_header_info(header, label);

  // Copy to device memory
  Kokkos::Impl::DeepCopy<Kokkos::Experimental::SYCLDeviceUSMSpace, HostSpace>(
      RecordBase::m_alloc_ptr, &header, sizeof(SharedAllocationHeader));
}

SharedAllocationRecord<Kokkos::Experimental::SYCLSharedUSMSpace, void>::
    SharedAllocationRecord(
        const Kokkos::Experimental::SYCLSharedUSMSpace& arg_space,
        const std::string& arg_label, const size_t arg_alloc_size,
        const SharedAllocationRecord<void, void>::function_type arg_dealloc)
    // Pass through allocated [ SharedAllocationHeader , user_memory ]
    // Pass through deallocation function
    : base_t(
#ifdef KOKKOS_ENABLE_DEBUG
          &SharedAllocationRecord<Kokkos::Experimental::SYCLSharedUSMSpace,
                                  void>::s_root_record,
#endif
          Impl::checked_allocation_with_header(arg_space, arg_label,
                                               arg_alloc_size),
          sizeof(SharedAllocationHeader) + arg_alloc_size, arg_dealloc),
      m_space(arg_space) {

  this->base_t::_fill_host_accessible_header_info(*base_t::m_alloc_ptr,
                                                  arg_label);
}

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

SharedAllocationRecord<Kokkos::Experimental::SYCLDeviceUSMSpace,
                       void>::~SharedAllocationRecord() {
  const char* label = nullptr;
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    SharedAllocationHeader header;
    Kokkos::Impl::DeepCopy<Kokkos::Experimental::SYCLDeviceUSMSpace,
                           Kokkos::HostSpace>(&header, RecordBase::m_alloc_ptr,
                                              sizeof(SharedAllocationHeader));
    label = header.label();
  }
  const auto alloc_size = SharedAllocationRecord<void, void>::m_alloc_size;
  m_space.deallocate(label, SharedAllocationRecord<void, void>::m_alloc_ptr,
                     alloc_size, alloc_size - sizeof(SharedAllocationHeader));
}

SharedAllocationRecord<Kokkos::Experimental::SYCLSharedUSMSpace,
                       void>::~SharedAllocationRecord() {
  const char* label = nullptr;
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    label = RecordBase::m_alloc_ptr->m_label;
  }
  const auto alloc_size = SharedAllocationRecord<void, void>::m_alloc_size;
  m_space.deallocate(label, SharedAllocationRecord<void, void>::m_alloc_ptr,
                     alloc_size, alloc_size - sizeof(SharedAllocationHeader));
}

//----------------------------------------------------------------------------

}  // namespace Impl
}  // namespace Kokkos

//==============================================================================
// <editor-fold desc="Explicit instantiations of CRTP Base classes"> {{{1

#include <impl/Kokkos_SharedAlloc_timpl.hpp>

namespace Kokkos {
namespace Impl {

// To avoid additional compilation cost for something that's (mostly?) not
// performance sensitive, we explicity instantiate these CRTP base classes here,
// where we have access to the associated *_timpl.hpp header files.
template class HostInaccessibleSharedAllocationRecordCommon<
    Kokkos::Experimental::SYCLDeviceUSMSpace>;
template class SharedAllocationRecordCommon<
    Kokkos::Experimental::SYCLDeviceUSMSpace>;
template class SharedAllocationRecordCommon<
    Kokkos::Experimental::SYCLSharedUSMSpace>;

}  // namespace Impl
}  // namespace Kokkos

// </editor-fold> end Explicit instantiations of CRTP Base classes }}}1
//==============================================================================
