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

#include <algorithm>
#include <omp.h>
#include <Kokkos_Macros.hpp>

/*--------------------------------------------------------------------------*/

#include <stdlib.h>
#include <stdint.h>
#include <memory.h>

#include <iostream>
#include <sstream>
#include <cstring>

#include <Kokkos_OpenMPTargetSpace.hpp>
#include <impl/Kokkos_Error.hpp>
#include <Kokkos_Atomic.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {
/* Default allocation mechanism */
OpenMPTargetSpace::OpenMPTargetSpace() {}

void *OpenMPTargetSpace::allocate(const size_t arg_alloc_size) const {
  static_assert(sizeof(void *) == sizeof(uintptr_t),
                "Error sizeof(void*) != sizeof(uintptr_t)");

  void *ptr;

  ptr = omp_target_alloc(arg_alloc_size, omp_get_default_device());

  return ptr;
}

void OpenMPTargetSpace::deallocate(void *const arg_alloc_ptr,
                                   const size_t /*arg_alloc_size*/) const {
  if (arg_alloc_ptr) {
    omp_target_free(arg_alloc_ptr, omp_get_default_device());
  }
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
  m_space.deallocate(SharedAllocationRecord<void, void>::m_alloc_ptr,
                     SharedAllocationRecord<void, void>::m_alloc_size);
}

// TODO: Implement deep copy back see CudaSpace
std::string SharedAllocationRecord<Kokkos::Experimental::OpenMPTargetSpace,
                                   void>::get_label() const {
  return std::string("OpenMPTargetAllocation");
}

void SharedAllocationRecord<Kokkos::Experimental::OpenMPTargetSpace,
                            void>::deallocate(SharedAllocationRecord<void, void>
                                                  *arg_rec) {
  delete static_cast<SharedAllocationRecord *>(arg_rec);
}

SharedAllocationRecord<Kokkos::Experimental::OpenMPTargetSpace, void>::
    SharedAllocationRecord(
        const Kokkos::Experimental::OpenMPTargetSpace &arg_space,
        const std::string &arg_label, const size_t arg_alloc_size,
        const SharedAllocationRecord<void, void>::function_type arg_dealloc)
    // Pass through allocated [ SharedAllocationHeader , user_memory ]
    // Pass through deallocation function
    : SharedAllocationRecord<void, void>(
#ifdef KOKKOS_ENABLE_DEBUG
          &SharedAllocationRecord<Kokkos::Experimental::OpenMPTargetSpace,
                                  void>::s_root_record,
#endif
          reinterpret_cast<SharedAllocationHeader *>(arg_space.allocate(
              sizeof(SharedAllocationHeader) + arg_alloc_size)),
          sizeof(SharedAllocationHeader) + arg_alloc_size, arg_dealloc),
      m_space(arg_space) {
  SharedAllocationHeader header;

  header.m_record = static_cast<SharedAllocationRecord<void, void> *>(this);

  strncpy(header.m_label, arg_label.c_str(),
          SharedAllocationHeader::maximum_label_length);
  // Set last element zero, in case c_str is too long
  header.m_label[SharedAllocationHeader::maximum_label_length - 1] = (char)0;
  // TODO DeepCopy
  // DeepCopy
  Kokkos::Impl::DeepCopy<Experimental::OpenMPTargetSpace, HostSpace>(
      RecordBase::m_alloc_ptr, &header, sizeof(SharedAllocationHeader));
}

//----------------------------------------------------------------------------

void *SharedAllocationRecord<Kokkos::Experimental::OpenMPTargetSpace, void>::
    allocate_tracked(const Kokkos::Experimental::OpenMPTargetSpace &arg_space,
                     const std::string &arg_alloc_label,
                     const size_t arg_alloc_size) {
  if (!arg_alloc_size) return nullptr;

  SharedAllocationRecord *const r =
      allocate(arg_space, arg_alloc_label, arg_alloc_size);

  RecordBase::increment(r);

  return r->data();
}

void SharedAllocationRecord<Kokkos::Experimental::OpenMPTargetSpace,
                            void>::deallocate_tracked(void *const
                                                          arg_alloc_ptr) {
  if (arg_alloc_ptr != nullptr) {
    SharedAllocationRecord *const r = get_record(arg_alloc_ptr);

    RecordBase::decrement(r);
  }
}

void *SharedAllocationRecord<Kokkos::Experimental::OpenMPTargetSpace, void>::
    reallocate_tracked(void *const arg_alloc_ptr, const size_t arg_alloc_size) {
  SharedAllocationRecord *const r_old = get_record(arg_alloc_ptr);
  SharedAllocationRecord *const r_new =
      allocate(r_old->m_space, r_old->get_label(), arg_alloc_size);

  // Kokkos::Impl::DeepCopy<OpenMPTargetSpace,OpenMPTargetSpace>( r_new->data()
  // , r_old->data()
  //                                           , std::min( r_old->size() ,
  //                                           r_new->size() ) );

  RecordBase::increment(r_new);
  RecordBase::decrement(r_old);

  return r_new->data();
}

SharedAllocationRecord<Kokkos::Experimental::OpenMPTargetSpace, void>
    *SharedAllocationRecord<Kokkos::Experimental::OpenMPTargetSpace,
                            void>::get_record(void *alloc_ptr) {
  using Header = SharedAllocationHeader;
  using RecordHost =
      SharedAllocationRecord<Kokkos::Experimental::OpenMPTargetSpace, void>;

  if (alloc_ptr) {
    Header head;
    const Header *const head_ompt = Header::get_header(alloc_ptr);

    Kokkos::Impl::DeepCopy<HostSpace, Experimental::OpenMPTargetSpace>(
        &head, head_ompt, sizeof(SharedAllocationHeader));

    RecordHost *record = static_cast<RecordHost *>(head.m_record);
    if (record->m_alloc_ptr == head_ompt) {
      return record;
    }
  }
  Kokkos::Impl::throw_runtime_exception(std::string(
      "Kokkos::Experimental::Impl::SharedAllocationRecord< "
      "Kokkos::Experimental::OpenMPTargetSpace , void >::get_record ERROR"));
  return nullptr;
}

// Iterate records to print orphaned memory ...
void SharedAllocationRecord<Kokkos::Experimental::OpenMPTargetSpace, void>::
    print_records(std::ostream &s,
                  const Kokkos::Experimental::OpenMPTargetSpace &,
                  bool detail) {
#ifdef KOKKOS_ENABLE_DEBUG
  SharedAllocationRecord<void, void>::print_host_accessible_records(
      s, "OpenMPTargetSpace", &s_root_record, detail);
#else
  (void)s;
  (void)detail;
  throw_runtime_exception(
      "SharedAllocationRecord<OpenMPTargetSpace>::print_records"
      " only works with KOKKOS_ENABLE_DEBUG enabled");
#endif
}

}  // namespace Impl
}  // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Impl {

template <class>
struct ViewOperatorBoundsErrorAbort;

template <>
struct ViewOperatorBoundsErrorAbort<Kokkos::Experimental::OpenMPTargetSpace> {
  static void apply(const size_t rank, const size_t n0, const size_t n1,
                    const size_t n2, const size_t n3, const size_t n4,
                    const size_t n5, const size_t n6, const size_t n7,
                    const size_t i0, const size_t i1, const size_t i2,
                    const size_t i3, const size_t i4, const size_t i5,
                    const size_t i6, const size_t i7);
};

void ViewOperatorBoundsErrorAbort<Kokkos::Experimental::OpenMPTargetSpace>::
    apply(const size_t rank, const size_t n0, const size_t n1, const size_t n2,
          const size_t n3, const size_t n4, const size_t n5, const size_t n6,
          const size_t n7, const size_t i0, const size_t i1, const size_t i2,
          const size_t i3, const size_t i4, const size_t i5, const size_t i6,
          const size_t i7) {
  printf(
      "View operator bounds error : rank(%lu) "
      "dim(%lu,%lu,%lu,%lu,%lu,%lu,%lu,%lu) "
      "index(%lu,%lu,%lu,%lu,%lu,%lu,%lu,%lu)",
      rank, n0, n1, n2, n3, n4, n5, n6, n7, i0, i1, i2, i3, i4, i5, i6, i7);
  // Kokkos::Impl::throw_runtime_exception( buffer );
}

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
