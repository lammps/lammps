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

/* only compile this file if ROCM is enabled for Kokkos */
#ifdef KOKKOS_ENABLE_ROCM

#include <Kokkos_Core.hpp>
#include <Kokkos_ROCm.hpp>
#include <Kokkos_ROCmSpace.hpp>

#include <impl/Kokkos_Error.hpp>

#if defined(KOKKOS_ENABLE_PROFILING)
#include <impl/Kokkos_Profiling_Interface.hpp>
#endif

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
#define ROCM_SAFE_CALL
namespace Kokkos {
namespace Impl {
using namespace hc;

DeepCopy<Kokkos::Experimental::ROCmSpace, Kokkos::Experimental::ROCmSpace,
         Kokkos::Experimental::ROCm>::DeepCopy(void* dst, const void* src,
                                               size_t n) {
  hc::accelerator acc;
  hc::accelerator_view av = acc.get_default_view();
  av.copy(src, dst, n);
}

DeepCopy<HostSpace, Kokkos::Experimental::ROCmSpace,
         Kokkos::Experimental::ROCm>::DeepCopy(void* dst, const void* src,
                                               size_t n) {
  hc::accelerator acc;
  hc::accelerator_view av = acc.get_default_view();
  av.copy(src, dst, n);
}

DeepCopy<Kokkos::Experimental::ROCmSpace, HostSpace,
         Kokkos::Experimental::ROCm>::DeepCopy(void* dst, const void* src,
                                               size_t n) {
  hc::accelerator acc;
  hc::accelerator_view av = acc.get_default_view();
  av.copy(src, dst, n);
}

DeepCopy<Kokkos::Experimental::ROCmSpace, Kokkos::Experimental::ROCmSpace,
         Kokkos::Experimental::ROCm>::DeepCopy(const Kokkos::Experimental::ROCm&
                                                   instance,
                                               void* dst, const void* src,
                                               size_t n) {
  hc::accelerator acc;
  hc::accelerator_view av = acc.get_default_view();
  av.copy(src, dst, n);
}

DeepCopy<HostSpace, Kokkos::Experimental::ROCmSpace,
         Kokkos::Experimental::ROCm>::DeepCopy(const Kokkos::Experimental::ROCm&
                                                   instance,
                                               void* dst, const void* src,
                                               size_t n) {
  hc::accelerator acc;
  hc::accelerator_view av = acc.get_default_view();
  av.copy(src, dst, n);
}

DeepCopy<Kokkos::Experimental::ROCmSpace, HostSpace,
         Kokkos::Experimental::ROCm>::DeepCopy(const Kokkos::Experimental::ROCm&
                                                   instance,
                                               void* dst, const void* src,
                                               size_t n) {
  hc::accelerator acc;
  hc::accelerator_view av = acc.get_default_view();
  av.copy(src, dst, n);
}

DeepCopy<Kokkos::Experimental::ROCmHostPinnedSpace,
         Kokkos::Experimental::ROCmHostPinnedSpace,
         Kokkos::Experimental::ROCm>::DeepCopy(void* dst, const void* src,
                                               size_t n) {
  hc::accelerator acc;
  hc::accelerator_view av = acc.get_default_view();
  av.copy(src, dst, n);
}

DeepCopy<HostSpace, Kokkos::Experimental::ROCmHostPinnedSpace,
         Kokkos::Experimental::ROCm>::DeepCopy(void* dst, const void* src,
                                               size_t n) {
  hc::accelerator acc;
  hc::accelerator_view av = acc.get_default_view();
  av.copy(src, dst, n);
}

DeepCopy<Kokkos::Experimental::ROCmHostPinnedSpace, HostSpace,
         Kokkos::Experimental::ROCm>::DeepCopy(void* dst, const void* src,
                                               size_t n) {
  hc::accelerator acc;
  hc::accelerator_view av = acc.get_default_view();
  av.copy(src, dst, n);
}

DeepCopy<Kokkos::Experimental::ROCmHostPinnedSpace,
         Kokkos::Experimental::ROCmHostPinnedSpace,
         Kokkos::Experimental::ROCm>::DeepCopy(const Kokkos::Experimental::ROCm&
                                                   instance,
                                               void* dst, const void* src,
                                               size_t n) {
  hc::accelerator acc;
  hc::accelerator_view av = acc.get_default_view();
  av.copy(src, dst, n);
}

DeepCopy<HostSpace, Kokkos::Experimental::ROCmHostPinnedSpace,
         Kokkos::Experimental::ROCm>::DeepCopy(const Kokkos::Experimental::ROCm&
                                                   instance,
                                               void* dst, const void* src,
                                               size_t n) {
  hc::accelerator acc;
  hc::accelerator_view av = acc.get_default_view();
  av.copy(src, dst, n);
}

DeepCopy<Kokkos::Experimental::ROCmHostPinnedSpace, HostSpace,
         Kokkos::Experimental::ROCm>::DeepCopy(const Kokkos::Experimental::ROCm&
                                                   instance,
                                               void* dst, const void* src,
                                               size_t n) {
  hc::accelerator acc;
  hc::accelerator_view av = acc.get_default_view();
  av.copy(src, dst, n);
}

hc::completion_future DeepCopyAsyncROCm(void* dst, const void* src, size_t n) {
  hc::accelerator acc;
  hc::accelerator_view av = acc.get_default_view();
  return (av.copy_async(src, dst, n));
}

}  // namespace Impl
}  // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {

void Experimental::ROCmSpace::access_error() {
  const std::string msg(
      "Kokkos::Experimental::ROCmSpace::access_error attempt to execute "
      "Experimental::ROCm function from non-ROCm space");
  Kokkos::Impl::throw_runtime_exception(msg);
}

void Experimental::ROCmSpace::access_error(const void* const) {
  const std::string msg(
      "Kokkos::Experimental::ROCmSpace::access_error attempt to execute "
      "Experimental::ROCm function from non-ROCm space");
  Kokkos::Impl::throw_runtime_exception(msg);
}

}  // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Experimental {

ROCmSpace::ROCmSpace() : m_device(ROCm().rocm_device()) {}

ROCmHostPinnedSpace::ROCmHostPinnedSpace() {}

void* ROCmSpace::allocate(const size_t arg_alloc_size) const {
  void* ptr = Kokkos::Impl::rocm_device_allocate(arg_alloc_size);
  return ptr;
}

void* Experimental::ROCmHostPinnedSpace::allocate(
    const size_t arg_alloc_size) const {
  void* ptr = Kokkos::Impl::rocm_hostpinned_allocate(arg_alloc_size);
  return ptr;
}

void ROCmSpace::deallocate(void* const arg_alloc_ptr,
                           const size_t /* arg_alloc_size */) const {
  Kokkos::Impl::rocm_device_free(arg_alloc_ptr);
}

void Experimental::ROCmHostPinnedSpace::deallocate(
    void* const arg_alloc_ptr, const size_t /* arg_alloc_size */) const {
  Kokkos::Impl::rocm_device_free(arg_alloc_ptr);
}

}  // namespace Experimental
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

#ifdef KOKKOS_DEBUG
SharedAllocationRecord<void, void> SharedAllocationRecord<
    Kokkos::Experimental::ROCmSpace, void>::s_root_record;

SharedAllocationRecord<void, void> SharedAllocationRecord<
    Kokkos::Experimental::ROCmHostPinnedSpace, void>::s_root_record;
#endif

std::string SharedAllocationRecord<Kokkos::Experimental::ROCmSpace,
                                   void>::get_label() const {
  SharedAllocationHeader header;

  Kokkos::Impl::DeepCopy<Kokkos::HostSpace, Kokkos::Experimental::ROCmSpace>(
      &header, RecordBase::head(), sizeof(SharedAllocationHeader));

  return std::string(header.m_label);
}

std::string SharedAllocationRecord<Kokkos::Experimental::ROCmHostPinnedSpace,
                                   void>::get_label() const {
  return std::string(RecordBase::head()->m_label);
}

SharedAllocationRecord<Kokkos::Experimental::ROCmSpace, void>*
SharedAllocationRecord<Kokkos::Experimental::ROCmSpace, void>::allocate(
    const Kokkos::Experimental::ROCmSpace& arg_space,
    const std::string& arg_label, const size_t arg_alloc_size) {
  return new SharedAllocationRecord(arg_space, arg_label, arg_alloc_size);
}

SharedAllocationRecord<Kokkos::Experimental::ROCmHostPinnedSpace, void>*
SharedAllocationRecord<Kokkos::Experimental::ROCmHostPinnedSpace, void>::
    allocate(const Kokkos::Experimental::ROCmHostPinnedSpace& arg_space,
             const std::string& arg_label, const size_t arg_alloc_size) {
  return new SharedAllocationRecord(arg_space, arg_label, arg_alloc_size);
}

void SharedAllocationRecord<Kokkos::Experimental::ROCmSpace, void>::deallocate(
    SharedAllocationRecord<void, void>* arg_rec) {
  delete static_cast<SharedAllocationRecord*>(arg_rec);
}

void SharedAllocationRecord<Kokkos::Experimental::ROCmHostPinnedSpace, void>::
    deallocate(SharedAllocationRecord<void, void>* arg_rec) {
  delete static_cast<SharedAllocationRecord*>(arg_rec);
}

SharedAllocationRecord<Kokkos::Experimental::ROCmSpace,
                       void>::~SharedAllocationRecord() {
#if defined(KOKKOS_ENABLE_PROFILING)
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    SharedAllocationHeader header;
    Kokkos::Impl::DeepCopy<Kokkos::Experimental::ROCmSpace, HostSpace>(
        &header, RecordBase::m_alloc_ptr, sizeof(SharedAllocationHeader));

    Kokkos::Profiling::deallocateData(
        Kokkos::Profiling::SpaceHandle(Kokkos::Experimental::ROCmSpace::name()),
        header.m_label, data(), size());
  }
#endif

  m_space.deallocate(SharedAllocationRecord<void, void>::m_alloc_ptr,
                     SharedAllocationRecord<void, void>::m_alloc_size);
}

SharedAllocationRecord<Kokkos::Experimental::ROCmHostPinnedSpace,
                       void>::~SharedAllocationRecord() {
#if defined(KOKKOS_ENABLE_PROFILING)
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    Kokkos::Profiling::deallocateData(
        Kokkos::Profiling::SpaceHandle(
            Kokkos::Experimental::ROCmHostPinnedSpace::name()),
        RecordBase::m_alloc_ptr->m_label, data(), size());
  }
#endif

  m_space.deallocate(SharedAllocationRecord<void, void>::m_alloc_ptr,
                     SharedAllocationRecord<void, void>::m_alloc_size);
}

SharedAllocationRecord<Kokkos::Experimental::ROCmSpace, void>::
    SharedAllocationRecord(
        const Kokkos::Experimental::ROCmSpace& arg_space,
        const std::string& arg_label, const size_t arg_alloc_size,
        const SharedAllocationRecord<void, void>::function_type arg_dealloc)
    // Pass through allocated [ SharedAllocationHeader , user_memory ]
    // Pass through deallocation function
    : SharedAllocationRecord<void, void>(
#ifdef KOKKOS_DEBUG
          &SharedAllocationRecord<Kokkos::Experimental::ROCmSpace,
                                  void>::s_root_record,
#endif
          reinterpret_cast<SharedAllocationHeader*>(arg_space.allocate(
              sizeof(SharedAllocationHeader) + arg_alloc_size)),
          sizeof(SharedAllocationHeader) + arg_alloc_size, arg_dealloc),
      m_space(arg_space) {
#if defined(KOKKOS_ENABLE_PROFILING)
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    Kokkos::Profiling::allocateData(
        Kokkos::Profiling::SpaceHandle(arg_space.name()), arg_label, data(),
        arg_alloc_size);
  }
#endif

  SharedAllocationHeader header;

  // Fill in the Header information
  header.m_record = static_cast<SharedAllocationRecord<void, void>*>(this);

  strncpy(header.m_label, arg_label.c_str(),
          SharedAllocationHeader::maximum_label_length);
  // Set last element zero, in case c_str is too long
  header.m_label[SharedAllocationHeader::maximum_label_length - 1] = (char)0;

  // Copy to device memory
  Kokkos::Impl::DeepCopy<Kokkos::Experimental::ROCmSpace, HostSpace>(
      RecordBase::m_alloc_ptr, &header, sizeof(SharedAllocationHeader));
}

SharedAllocationRecord<Kokkos::Experimental::ROCmHostPinnedSpace, void>::
    SharedAllocationRecord(
        const Kokkos::Experimental::ROCmHostPinnedSpace& arg_space,
        const std::string& arg_label, const size_t arg_alloc_size,
        const SharedAllocationRecord<void, void>::function_type arg_dealloc)
    // Pass through allocated [ SharedAllocationHeader , user_memory ]
    // Pass through deallocation function
    : SharedAllocationRecord<void, void>(
#ifdef KOKKOS_DEBUG
          &SharedAllocationRecord<Kokkos::Experimental::ROCmHostPinnedSpace,
                                  void>::s_root_record,
#endif
          reinterpret_cast<SharedAllocationHeader*>(arg_space.allocate(
              sizeof(SharedAllocationHeader) + arg_alloc_size)),
          sizeof(SharedAllocationHeader) + arg_alloc_size, arg_dealloc),
      m_space(arg_space) {
#if defined(KOKKOS_ENABLE_PROFILING)
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    Kokkos::Profiling::allocateData(
        Kokkos::Profiling::SpaceHandle(arg_space.name()), arg_label, data(),
        arg_alloc_size);
  }
#endif
  // Fill in the Header information, directly accessible via host pinned memory

  RecordBase::m_alloc_ptr->m_record = this;

  strncpy(RecordBase::m_alloc_ptr->m_label, arg_label.c_str(),
          SharedAllocationHeader::maximum_label_length);
  // Set last element zero, in case c_str is too long
  RecordBase::m_alloc_ptr
      ->m_label[SharedAllocationHeader::maximum_label_length - 1] = (char)0;
}

//----------------------------------------------------------------------------

void* SharedAllocationRecord<Kokkos::Experimental::ROCmSpace, void>::
    allocate_tracked(const Kokkos::Experimental::ROCmSpace& arg_space,
                     const std::string& arg_alloc_label,
                     const size_t arg_alloc_size) {
  if (!arg_alloc_size) return nullptr;

  SharedAllocationRecord* const r =
      allocate(arg_space, arg_alloc_label, arg_alloc_size);

  RecordBase::increment(r);

  return r->data();
}

void SharedAllocationRecord<Kokkos::Experimental::ROCmSpace,
                            void>::deallocate_tracked(void* const
                                                          arg_alloc_ptr) {
  if (arg_alloc_ptr != 0) {
    SharedAllocationRecord* const r = get_record(arg_alloc_ptr);

    RecordBase::decrement(r);
  }
}

void* SharedAllocationRecord<Kokkos::Experimental::ROCmSpace, void>::
    reallocate_tracked(void* const arg_alloc_ptr, const size_t arg_alloc_size) {
  SharedAllocationRecord* const r_old = get_record(arg_alloc_ptr);
  SharedAllocationRecord* const r_new =
      allocate(r_old->m_space, r_old->get_label(), arg_alloc_size);

  Kokkos::Impl::DeepCopy<Kokkos::Experimental::ROCmSpace,
                         Kokkos::Experimental::ROCmSpace>(
      r_new->data(), r_old->data(), std::min(r_old->size(), r_new->size()));

  RecordBase::increment(r_new);
  RecordBase::decrement(r_old);

  return r_new->data();
}

#if 0
void * SharedAllocationRecord< Kokkos::Experimental::ROCmHostPinnedSpace , void >::
allocate_tracked( const Kokkos::Experimental::ROCmHostPinnedSpace & arg_space
                , const std::string & arg_alloc_label
                , const size_t arg_alloc_size )
{
  if ( ! arg_alloc_size ) return (void *) 0 ;

  SharedAllocationRecord * const r =
    allocate( arg_space , arg_alloc_label , arg_alloc_size );

  RecordBase::increment( r );

  return r->data();
}

void SharedAllocationRecord< Kokkos::Experimental::ROCmHostPinnedSpace , void >::
deallocate_tracked( void * const arg_alloc_ptr )
{
  if ( arg_alloc_ptr != 0 ) {
    SharedAllocationRecord * const r = get_record( arg_alloc_ptr );

    RecordBase::decrement( r );
  }
}

void * SharedAllocationRecord< Kokkos::Experimental::ROCmHostPinnedSpace , void >::
reallocate_tracked( void * const arg_alloc_ptr
                  , const size_t arg_alloc_size )
{
  SharedAllocationRecord * const r_old = get_record( arg_alloc_ptr );
  SharedAllocationRecord * const r_new = allocate( r_old->m_space , r_old->get_label() , arg_alloc_size );

  Kokkos::Impl::DeepCopy<Experimental::ROCmHostPinnedSpace,Experimental::ROCmHostPinnedSpace>( r_new->data() , r_old->data()
                                             , std::min( r_old->size() , r_new->size() ) );

  RecordBase::increment( r_new );
  RecordBase::decrement( r_old );

  return r_new->data();
}
#endif

//----------------------------------------------------------------------------

SharedAllocationRecord<Kokkos::Experimental::ROCmSpace, void>*
SharedAllocationRecord<Kokkos::Experimental::ROCmSpace, void>::get_record(
    void* alloc_ptr) {
  using Header     = SharedAllocationHeader;
  using RecordBase = SharedAllocationRecord<void, void>;
  using RecordROCm =
      SharedAllocationRecord<Kokkos::Experimental::ROCmSpace, void>;

  // Copy the header from the allocation
  Header head;

  Header const* const head_rocm =
      alloc_ptr ? Header::get_header(alloc_ptr) : (Header*)0;

  if (alloc_ptr) {
    Kokkos::Impl::DeepCopy<HostSpace, Kokkos::Experimental::ROCmSpace>(
        &head, head_rocm, sizeof(SharedAllocationHeader));
  }

  RecordROCm* const record =
      alloc_ptr ? static_cast<RecordROCm*>(head.m_record) : (RecordROCm*)0;

  if (!alloc_ptr || record->m_alloc_ptr != head_rocm) {
    Kokkos::Impl::throw_runtime_exception(std::string(
        "Kokkos::Impl::SharedAllocationRecord< Kokkos::Experimental::ROCmSpace "
        ", void >::get_record ERROR"));
  }

  return record;
}

#if 0
SharedAllocationRecord< Kokkos::Experimental::ROCmHostPinnedSpace , void > *
SharedAllocationRecord< Kokkos::Experimental::ROCmHostPinnedSpace , void >::get_record( void * alloc_ptr )
{
  using Header     = SharedAllocationHeader ;
  using RecordROCm = SharedAllocationRecord< Kokkos::Experimental::ROCmHostPinnedSpace , void > ;

  Header * const h = alloc_ptr ? reinterpret_cast< Header * >( alloc_ptr ) - 1 : (Header *) 0 ;

  if ( ! alloc_ptr || h->m_record->m_alloc_ptr != h ) {
    Kokkos::Impl::throw_runtime_exception( std::string("Kokkos::Impl::SharedAllocationRecord< Kokkos::Experimental::ROCmHostPinnedSpace , void >::get_record ERROR" ) );
  }

  return static_cast< RecordROCm * >( h->m_record );
}
#endif

// Iterate records to print orphaned memory ...
void SharedAllocationRecord<Kokkos::Experimental::ROCmSpace, void>::
    print_records(std::ostream& s, const Kokkos::Experimental::ROCmSpace& space,
                  bool detail) {
#ifdef KOKKOS_DEBUG
  SharedAllocationRecord<void, void>* r = &s_root_record;

  char buffer[256];

  SharedAllocationHeader head;

  if (detail) {
    do {
      if (r->m_alloc_ptr) {
        Kokkos::Impl::DeepCopy<HostSpace, Kokkos::Experimental::ROCmSpace>(
            &head, r->m_alloc_ptr, sizeof(SharedAllocationHeader));
      } else {
        head.m_label[0] = 0;
      }

      // Formatting dependent on sizeof(uintptr_t)
      const char* format_string;

      if (sizeof(uintptr_t) == sizeof(unsigned long)) {
        format_string =
            "ROCm addr( 0x%.12lx ) list( 0x%.12lx 0x%.12lx ) extent[ 0x%.12lx "
            "+ %.8ld ] count(%d) dealloc(0x%.12lx) %s\n";
      } else if (sizeof(uintptr_t) == sizeof(unsigned long long)) {
        format_string =
            "ROCm addr( 0x%.12llx ) list( 0x%.12llx 0x%.12llx ) extent[ "
            "0x%.12llx + %.8ld ] count(%d) dealloc(0x%.12llx) %s\n";
      }

      snprintf(buffer, 256, format_string, reinterpret_cast<uintptr_t>(r),
               reinterpret_cast<uintptr_t>(r->m_prev),
               reinterpret_cast<uintptr_t>(r->m_next),
               reinterpret_cast<uintptr_t>(r->m_alloc_ptr), r->m_alloc_size,
               r->m_count, reinterpret_cast<uintptr_t>(r->m_dealloc),
               head.m_label);
      std::cout << buffer;
      r = r->m_next;
    } while (r != &s_root_record);
  } else {
    do {
      if (r->m_alloc_ptr) {
        Kokkos::Impl::DeepCopy<HostSpace, Kokkos::Experimental::ROCmSpace>(
            &head, r->m_alloc_ptr, sizeof(SharedAllocationHeader));

        // Formatting dependent on sizeof(uintptr_t)
        const char* format_string;

        if (sizeof(uintptr_t) == sizeof(unsigned long)) {
          format_string = "ROCm [ 0x%.12lx + %ld ] %s\n";
        } else if (sizeof(uintptr_t) == sizeof(unsigned long long)) {
          format_string = "ROCm [ 0x%.12llx + %ld ] %s\n";
        }

        snprintf(buffer, 256, format_string,
                 reinterpret_cast<uintptr_t>(r->data()), r->size(),
                 head.m_label);
      } else {
        snprintf(buffer, 256, "ROCm [ 0 + 0 ]\n");
      }
      std::cout << buffer;
      r = r->m_next;
    } while (r != &s_root_record);
  }
#else
  throw_runtime_exception(
      "Kokkos::Impl::SharedAllocationRecord<ROCmSpace>::print_records"
      " only works with KOKKOS_DEBUG enabled");
#endif
}

}  // namespace Impl
}  // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
namespace Kokkos {
namespace {

void* rocm_resize_scratch_space(size_t bytes, bool force_shrink) {
  static void* ptr           = nullptr;
  static size_t current_size = 0;
  if (current_size == 0) {
    current_size = bytes;
    ptr          = Kokkos::kokkos_malloc<Kokkos::Experimental::ROCmSpace>(
        "ROCmSpace::ScratchMemory", current_size);
  }
  if (bytes > current_size) {
    current_size = bytes;
    ptr          = Kokkos::kokkos_realloc<Kokkos::Experimental::ROCmSpace>(ptr,
                                                                  current_size);
  }
  if ((bytes < current_size) && (force_shrink)) {
    current_size = bytes;
    Kokkos::kokkos_free<Kokkos::Experimental::ROCmSpace>(ptr);
    ptr = Kokkos::kokkos_malloc<Kokkos::Experimental::ROCmSpace>(
        "ROCmSpace::ScratchMemory", current_size);
  }
  return ptr;
}

}  // namespace
}  // namespace Kokkos

#endif  // KOKKOS_ENABLE_ROCM
