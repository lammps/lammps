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

#include <cstddef>
#include <cstdlib>
#include <cstdint>
#include <cstring>

#include <iostream>
#include <sstream>
#include <cstring>
#include <algorithm>

#include <Kokkos_HBWSpace.hpp>
#include <impl/Kokkos_Error.hpp>
#include <impl/Kokkos_MemorySpace.hpp>
#include <Kokkos_Atomic.hpp>
#ifdef KOKKOS_ENABLE_HBWSPACE
#include <memkind.h>
#endif

#include <impl/Kokkos_Tools.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
#ifdef KOKKOS_ENABLE_HBWSPACE
#define MEMKIND_TYPE MEMKIND_HBW  // hbw_get_kind(HBW_PAGESIZE_4KB)

/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Experimental {

/* Default allocation mechanism */
HBWSpace::HBWSpace() : m_alloc_mech(HBWSpace::STD_MALLOC) {
  printf("Init\n");
  setenv("MEMKIND_HBW_NODES", "1", 0);
}

/* Default allocation mechanism */
HBWSpace::HBWSpace(const HBWSpace::AllocationMechanism &arg_alloc_mech)
    : m_alloc_mech(HBWSpace::STD_MALLOC) {
  printf("Init2\n");
  setenv("MEMKIND_HBW_NODES", "1", 0);
  if (arg_alloc_mech == STD_MALLOC) {
    m_alloc_mech = HBWSpace::STD_MALLOC;
  }
}

void *HBWSpace::allocate(const size_t arg_alloc_size) const {
  return allocate("[unlabeled]", arg_alloc_size);
}
void *HBWSpace::allocate(const char *arg_label, const size_t arg_alloc_size,
                         const size_t arg_logical_size) const {
  return impl_allocate(arg_label, arg_alloc_size, arg_logical_size);
}
void *HBWSpace::impl_allocate(
    const char *arg_label, const size_t arg_alloc_size,
    const size_t arg_logical_size,
    const Kokkos::Tools::SpaceHandle arg_handle) const {
  static_assert(sizeof(void *) == sizeof(uintptr_t),
                "Error sizeof(void*) != sizeof(uintptr_t)");

  static_assert(
      Kokkos::Impl::power_of_two<Kokkos::Impl::MEMORY_ALIGNMENT>::value,
      "Memory alignment must be power of two");

  constexpr uintptr_t alignment      = Kokkos::Impl::MEMORY_ALIGNMENT;
  constexpr uintptr_t alignment_mask = alignment - 1;

  void *ptr = nullptr;

  if (arg_alloc_size) {
    if (m_alloc_mech == STD_MALLOC) {
      // Over-allocate to and round up to guarantee proper alignment.
      size_t size_padded = arg_alloc_size + sizeof(void *) + alignment;

      void *alloc_ptr = memkind_malloc(MEMKIND_TYPE, size_padded);

      if (alloc_ptr) {
        uintptr_t address = reinterpret_cast<uintptr_t>(alloc_ptr);

        // offset enough to record the alloc_ptr
        address += sizeof(void *);
        uintptr_t rem    = address % alignment;
        uintptr_t offset = rem ? (alignment - rem) : 0u;
        address += offset;
        ptr = reinterpret_cast<void *>(address);
        // record the alloc'd pointer
        address -= sizeof(void *);
        *reinterpret_cast<void **>(address) = alloc_ptr;
      }
    }
  }

  if ((ptr == nullptr) || (reinterpret_cast<uintptr_t>(ptr) == ~uintptr_t(0)) ||
      (reinterpret_cast<uintptr_t>(ptr) & alignment_mask)) {
    std::ostringstream msg;
    msg << "Kokkos::Experimental::HBWSpace::allocate[ ";
    switch (m_alloc_mech) {
      case STD_MALLOC: msg << "STD_MALLOC"; break;
      case POSIX_MEMALIGN: msg << "POSIX_MEMALIGN"; break;
      case POSIX_MMAP: msg << "POSIX_MMAP"; break;
      case INTEL_MM_ALLOC: msg << "INTEL_MM_ALLOC"; break;
    }
    msg << " ]( " << arg_alloc_size << " ) FAILED";
    if (ptr == nullptr) {
      msg << " nullptr";
    } else {
      msg << " NOT ALIGNED " << ptr;
    }

    std::cerr << msg.str() << std::endl;
    std::cerr.flush();

    Kokkos::Impl::throw_runtime_exception(msg.str());
  }
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    const size_t reported_size =
        (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;
    Kokkos::Profiling::allocateData(arg_handle, arg_label, ptr, reported_size);
  }

  return ptr;
}

void HBWSpace::deallocate(void *const arg_alloc_ptr,
                          const size_t arg_alloc_size) const {
  deallocate("[unlabeled]", arg_alloc_ptr, arg_alloc_size);
}
void HBWSpace::deallocate(const char *arg_label, void *const arg_alloc_ptr,
                          const size_t arg_alloc_size,
                          const size_t arg_logical_size) const {
  impl_deallocate(arg_label, arg_alloc_ptr, arg_alloc_size, arg_logical_size);
}
void HBWSpace::impl_deallocate(
    const char *arg_label, void *const arg_alloc_ptr,
    const size_t arg_alloc_size, const size_t arg_logical_size,
    const Kokkos::Tools::SpaceHandle arg_handle) const {
  if (arg_alloc_ptr) {
    if (Kokkos::Profiling::profileLibraryLoaded()) {
      const size_t reported_size =
          (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;
      Kokkos::Profiling::deallocateData(arg_handle, arg_label, arg_alloc_ptr,
                                        reported_size);
    }

    if (m_alloc_mech == STD_MALLOC) {
      void *alloc_ptr = *(reinterpret_cast<void **>(arg_alloc_ptr) - 1);
      memkind_free(MEMKIND_TYPE, alloc_ptr);
    }
  }
}

}  // namespace Experimental
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

#ifdef KOKKOS_ENABLE_DEBUG
SharedAllocationRecord<void, void>
    SharedAllocationRecord<Kokkos::Experimental::HBWSpace, void>::s_root_record;
#endif

void SharedAllocationRecord<Kokkos::Experimental::HBWSpace, void>::deallocate(
    SharedAllocationRecord<void, void> *arg_rec) {
  delete static_cast<SharedAllocationRecord *>(arg_rec);
}

SharedAllocationRecord<Kokkos::Experimental::HBWSpace,
                       void>::~SharedAllocationRecord()
#if defined( \
    KOKKOS_IMPL_INTEL_WORKAROUND_NOEXCEPT_SPECIFICATION_VIRTUAL_FUNCTION)
    noexcept
#endif
{

  m_space.deallocate(RecordBase::m_alloc_ptr->m_label,
                     SharedAllocationRecord<void, void>::m_alloc_ptr,
                     SharedAllocationRecord<void, void>::m_alloc_size,
                     (SharedAllocationRecord<void, void>::m_alloc_size -
                      sizeof(SharedAllocationHeader)));
}

SharedAllocationRecord<Kokkos::Experimental::HBWSpace, void>::
    SharedAllocationRecord(
        const Kokkos::Experimental::HBWSpace &arg_space,
        const std::string &arg_label, const size_t arg_alloc_size,
        const SharedAllocationRecord<void, void>::function_type arg_dealloc)
    // Pass through allocated [ SharedAllocationHeader , user_memory ]
    // Pass through deallocation function
    : SharedAllocationRecord<void, void>(
#ifdef KOKKOS_ENABLE_DEBUG
          &SharedAllocationRecord<Kokkos::Experimental::HBWSpace,
                                  void>::s_root_record,
#endif
          Impl::checked_allocation_with_header(arg_space, arg_label,
                                               arg_alloc_size),
          sizeof(SharedAllocationHeader) + arg_alloc_size, arg_dealloc),
      m_space(arg_space) {
  // Fill in the Header information
  RecordBase::m_alloc_ptr->m_record =
      static_cast<SharedAllocationRecord<void, void> *>(this);

  strncpy(RecordBase::m_alloc_ptr->m_label, arg_label.c_str(),
          SharedAllocationHeader::maximum_label_length);
  // Set last element zero, in case c_str is too long
  RecordBase::m_alloc_ptr
      ->m_label[SharedAllocationHeader::maximum_label_length - 1] = (char)0;
}

//----------------------------------------------------------------------------

void *
SharedAllocationRecord<Kokkos::Experimental::HBWSpace, void>::allocate_tracked(
    const Kokkos::Experimental::HBWSpace &arg_space,
    const std::string &arg_alloc_label, const size_t arg_alloc_size) {
  if (!arg_alloc_size) return nullptr;

  SharedAllocationRecord *const r =
      allocate(arg_space, arg_alloc_label, arg_alloc_size);

  RecordBase::increment(r);

  return r->data();
}

void SharedAllocationRecord<Kokkos::Experimental::HBWSpace,
                            void>::deallocate_tracked(void *const
                                                          arg_alloc_ptr) {
  if (arg_alloc_ptr != nullptr) {
    SharedAllocationRecord *const r = get_record(arg_alloc_ptr);

    RecordBase::decrement(r);
  }
}

void *SharedAllocationRecord<Kokkos::Experimental::HBWSpace, void>::
    reallocate_tracked(void *const arg_alloc_ptr, const size_t arg_alloc_size) {
  SharedAllocationRecord *const r_old = get_record(arg_alloc_ptr);
  SharedAllocationRecord *const r_new =
      allocate(r_old->m_space, r_old->get_label(), arg_alloc_size);

  Kokkos::Impl::DeepCopy<Kokkos::Experimental::HBWSpace,
                         Kokkos::Experimental::HBWSpace>(
      r_new->data(), r_old->data(), std::min(r_old->size(), r_new->size()));

  RecordBase::increment(r_new);
  RecordBase::decrement(r_old);

  return r_new->data();
}

SharedAllocationRecord<Kokkos::Experimental::HBWSpace, void>
    *SharedAllocationRecord<Kokkos::Experimental::HBWSpace, void>::get_record(
        void *alloc_ptr) {
  using Header = SharedAllocationHeader;
  using RecordHost =
      SharedAllocationRecord<Kokkos::Experimental::HBWSpace, void>;

  SharedAllocationHeader const *const head =
      alloc_ptr ? Header::get_header(alloc_ptr) : nullptr;
  RecordHost *const record =
      head ? static_cast<RecordHost *>(head->m_record) : nullptr;

  if (!alloc_ptr || record->m_alloc_ptr != head) {
    Kokkos::Impl::throw_runtime_exception(std::string(
        "Kokkos::Impl::SharedAllocationRecord< Kokkos::Experimental::HBWSpace "
        ", void >::get_record ERROR"));
  }

  return record;
}

// Iterate records to print orphaned memory ...
void SharedAllocationRecord<Kokkos::Experimental::HBWSpace, void>::
    print_records(std::ostream &s, const Kokkos::Experimental::HBWSpace &space,
                  bool detail) {
#ifdef KOKKOS_ENABLE_DEBUG
  SharedAllocationRecord<void, void>::print_host_accessible_records(
      s, "HBWSpace", &s_root_record, detail);
#else
  throw_runtime_exception(
      "SharedAllocationRecord<HBWSpace>::print_records"
      " only works with KOKKOS_ENABLE_DEBUG enabled");
#endif
}

}  // namespace Impl
}  // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace Experimental {
namespace {
const unsigned HBW_SPACE_ATOMIC_MASK     = 0xFFFF;
const unsigned HBW_SPACE_ATOMIC_XOR_MASK = 0x5A39;
static int HBW_SPACE_ATOMIC_LOCKS[HBW_SPACE_ATOMIC_MASK + 1];
}  // namespace

namespace Impl {
void init_lock_array_hbw_space() {
  static int is_initialized = 0;
  if (!is_initialized)
    for (int i = 0; i < static_cast<int>(HBW_SPACE_ATOMIC_MASK + 1); i++)
      HBW_SPACE_ATOMIC_LOCKS[i] = 0;
}

bool lock_address_hbw_space(void *ptr) {
  return 0 == atomic_compare_exchange(
                  &HBW_SPACE_ATOMIC_LOCKS[((size_t(ptr) >> 2) &
                                           HBW_SPACE_ATOMIC_MASK) ^
                                          HBW_SPACE_ATOMIC_XOR_MASK],
                  0, 1);
}

void unlock_address_hbw_space(void *ptr) {
  atomic_exchange(
      &HBW_SPACE_ATOMIC_LOCKS[((size_t(ptr) >> 2) & HBW_SPACE_ATOMIC_MASK) ^
                              HBW_SPACE_ATOMIC_XOR_MASK],
      0);
}

}  // namespace Impl
}  // namespace Experimental
}  // namespace Kokkos
#endif
