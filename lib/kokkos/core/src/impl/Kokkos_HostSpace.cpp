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

#include <cstdio>
#include <algorithm>
#include <Kokkos_Macros.hpp>
#include <impl/Kokkos_Error.hpp>
#include <impl/Kokkos_MemorySpace.hpp>
#if defined(KOKKOS_ENABLE_PROFILING)
#include <impl/Kokkos_Profiling_Interface.hpp>
#endif

/*--------------------------------------------------------------------------*/

#if defined(__INTEL_COMPILER) && !defined(KOKKOS_ENABLE_CUDA)

// Intel specialized allocator does not interoperate with CUDA memory allocation

#define KOKKOS_ENABLE_INTEL_MM_ALLOC

#endif

/*--------------------------------------------------------------------------*/

#if defined(KOKKOS_ENABLE_POSIX_MEMALIGN)

#include <unistd.h>
#include <sys/mman.h>

/* mmap flags for private anonymous memory allocation */

#if defined(MAP_ANONYMOUS) && defined(MAP_PRIVATE)
#define KOKKOS_IMPL_POSIX_MMAP_FLAGS (MAP_PRIVATE | MAP_ANONYMOUS)
#elif defined(MAP_ANON) && defined(MAP_PRIVATE)
#define KOKKOS_IMPL_POSIX_MMAP_FLAGS (MAP_PRIVATE | MAP_ANON)
#endif

// mmap flags for huge page tables
// the Cuda driver does not interoperate with MAP_HUGETLB
#if defined(KOKKOS_IMPL_POSIX_MMAP_FLAGS)
#if defined(MAP_HUGETLB) && !defined(KOKKOS_ENABLE_CUDA)
#define KOKKOS_IMPL_POSIX_MMAP_FLAGS_HUGE \
  (KOKKOS_IMPL_POSIX_MMAP_FLAGS | MAP_HUGETLB)
#else
#define KOKKOS_IMPL_POSIX_MMAP_FLAGS_HUGE KOKKOS_IMPL_POSIX_MMAP_FLAGS
#endif
#endif

#endif

/*--------------------------------------------------------------------------*/

#include <cstddef>
#include <cstdlib>
#include <cstdint>
#include <cstring>

#include <iostream>
#include <sstream>
#include <cstring>

#include <Kokkos_HostSpace.hpp>
#include <impl/Kokkos_Error.hpp>
#include <Kokkos_Atomic.hpp>

#if (defined(KOKKOS_ENABLE_ASM) || defined(KOKKOS_ENABLE_TM)) && \
    defined(KOKKOS_ENABLE_ISA_X86_64) && !defined(KOKKOS_COMPILER_PGI)
#include <immintrin.h>
#endif

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

/* Default allocation mechanism */
HostSpace::HostSpace()
    : m_alloc_mech(
#if defined(KOKKOS_ENABLE_INTEL_MM_ALLOC)
          HostSpace::INTEL_MM_ALLOC
#elif defined(KOKKOS_IMPL_POSIX_MMAP_FLAGS)
          HostSpace::POSIX_MMAP
#elif defined(KOKKOS_ENABLE_POSIX_MEMALIGN)
          HostSpace::POSIX_MEMALIGN
#else
          HostSpace::STD_MALLOC
#endif
      ) {
}

/* Default allocation mechanism */
HostSpace::HostSpace(const HostSpace::AllocationMechanism &arg_alloc_mech)
    : m_alloc_mech(HostSpace::STD_MALLOC) {
  if (arg_alloc_mech == STD_MALLOC) {
    m_alloc_mech = HostSpace::STD_MALLOC;
  }
#if defined(KOKKOS_ENABLE_INTEL_MM_ALLOC)
  else if (arg_alloc_mech == HostSpace::INTEL_MM_ALLOC) {
    m_alloc_mech = HostSpace::INTEL_MM_ALLOC;
  }
#elif defined(KOKKOS_ENABLE_POSIX_MEMALIGN)
  else if (arg_alloc_mech == HostSpace::POSIX_MEMALIGN) {
    m_alloc_mech = HostSpace::POSIX_MEMALIGN;
  }
#elif defined(KOKKOS_IMPL_POSIX_MMAP_FLAGS)
  else if (arg_alloc_mech == HostSpace::POSIX_MMAP) {
    m_alloc_mech = HostSpace::POSIX_MMAP;
  }
#endif
  else {
    const char *const mech =
        (arg_alloc_mech == HostSpace::INTEL_MM_ALLOC)
            ? "INTEL_MM_ALLOC"
            : ((arg_alloc_mech == HostSpace::POSIX_MEMALIGN)
                   ? "POSIX_MEMALIGN"
                   : ((arg_alloc_mech == HostSpace::POSIX_MMAP) ? "POSIX_MMAP"
                                                                : ""));

    std::string msg;
    msg.append("Kokkos::HostSpace ");
    msg.append(mech);
    msg.append(" is not available");
    Kokkos::Impl::throw_runtime_exception(msg);
  }
}

void *HostSpace::allocate(const size_t arg_alloc_size) const {
  static_assert(sizeof(void *) == sizeof(uintptr_t),
                "Error sizeof(void*) != sizeof(uintptr_t)");

  static_assert(
      Kokkos::Impl::is_integral_power_of_two(Kokkos::Impl::MEMORY_ALIGNMENT),
      "Memory alignment must be power of two");

  constexpr uintptr_t alignment      = Kokkos::Impl::MEMORY_ALIGNMENT;
  constexpr uintptr_t alignment_mask = alignment - 1;

  void *ptr = nullptr;

  if (arg_alloc_size) {
    if (m_alloc_mech == STD_MALLOC) {
      // Over-allocate to and round up to guarantee proper alignment.
      size_t size_padded = arg_alloc_size + sizeof(void *) + alignment;

      void *alloc_ptr = malloc(size_padded);

      if (alloc_ptr) {
        auto address = reinterpret_cast<uintptr_t>(alloc_ptr);

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
#if defined(KOKKOS_ENABLE_INTEL_MM_ALLOC)
    else if (m_alloc_mech == INTEL_MM_ALLOC) {
      ptr = _mm_malloc(arg_alloc_size, alignment);
    }
#endif

#if defined(KOKKOS_ENABLE_POSIX_MEMALIGN)
    else if (m_alloc_mech == POSIX_MEMALIGN) {
      posix_memalign(&ptr, alignment, arg_alloc_size);
    }
#endif

#if defined(KOKKOS_IMPL_POSIX_MMAP_FLAGS)
    else if (m_alloc_mech == POSIX_MMAP) {
      constexpr size_t use_huge_pages = (1u << 27);
      constexpr int prot              = PROT_READ | PROT_WRITE;
      const int flags                 = arg_alloc_size < use_huge_pages
                            ? KOKKOS_IMPL_POSIX_MMAP_FLAGS
                            : KOKKOS_IMPL_POSIX_MMAP_FLAGS_HUGE;

      // read write access to private memory

      ptr =
          mmap(nullptr /* address hint, if nullptr OS kernel chooses address */
               ,
               arg_alloc_size /* size in bytes */
               ,
               prot /* memory protection */
               ,
               flags /* visibility of updates */
               ,
               -1 /* file descriptor */
               ,
               0 /* offset */
          );

      /* Associated reallocation:
             ptr = mremap( old_ptr , old_size , new_size , MREMAP_MAYMOVE );
      */
    }
#endif
  }

  if ((ptr == nullptr) || (reinterpret_cast<uintptr_t>(ptr) == ~uintptr_t(0)) ||
      (reinterpret_cast<uintptr_t>(ptr) & alignment_mask)) {
    Experimental::RawMemoryAllocationFailure::FailureMode failure_mode =
        Experimental::RawMemoryAllocationFailure::FailureMode::
            AllocationNotAligned;
    if (ptr == nullptr) {
      failure_mode = Experimental::RawMemoryAllocationFailure::FailureMode::
          OutOfMemoryError;
    }

    Experimental::RawMemoryAllocationFailure::AllocationMechanism alloc_mec =
        Experimental::RawMemoryAllocationFailure::AllocationMechanism::
            StdMalloc;
    switch (m_alloc_mech) {
      case STD_MALLOC: break;  // default
      case POSIX_MEMALIGN:
        alloc_mec = Experimental::RawMemoryAllocationFailure::
            AllocationMechanism::PosixMemAlign;
        break;
      case POSIX_MMAP:
        alloc_mec = Experimental::RawMemoryAllocationFailure::
            AllocationMechanism::PosixMMap;
        break;
      case INTEL_MM_ALLOC:
        alloc_mec = Experimental::RawMemoryAllocationFailure::
            AllocationMechanism::IntelMMAlloc;
        break;
    }

    throw Kokkos::Experimental::RawMemoryAllocationFailure(
        arg_alloc_size, alignment, failure_mode, alloc_mec);
  }

  return ptr;
}

void HostSpace::deallocate(void *const arg_alloc_ptr, const size_t
#if defined(KOKKOS_IMPL_POSIX_MMAP_FLAGS)
                                                          arg_alloc_size
#endif
                           ) const {
  if (arg_alloc_ptr) {
    if (m_alloc_mech == STD_MALLOC) {
      void *alloc_ptr = *(reinterpret_cast<void **>(arg_alloc_ptr) - 1);
      free(alloc_ptr);
    }
#if defined(KOKKOS_ENABLE_INTEL_MM_ALLOC)
    else if (m_alloc_mech == INTEL_MM_ALLOC) {
      _mm_free(arg_alloc_ptr);
    }
#endif

#if defined(KOKKOS_ENABLE_POSIX_MEMALIGN)
    else if (m_alloc_mech == POSIX_MEMALIGN) {
      free(arg_alloc_ptr);
    }
#endif

#if defined(KOKKOS_IMPL_POSIX_MMAP_FLAGS)
    else if (m_alloc_mech == POSIX_MMAP) {
      munmap(arg_alloc_ptr, arg_alloc_size);
    }
#endif
  }
}

}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

#ifdef KOKKOS_DEBUG
SharedAllocationRecord<void, void>
    SharedAllocationRecord<Kokkos::HostSpace, void>::s_root_record;
#endif

void SharedAllocationRecord<Kokkos::HostSpace, void>::deallocate(
    SharedAllocationRecord<void, void> *arg_rec) {
  delete static_cast<SharedAllocationRecord *>(arg_rec);
}

SharedAllocationRecord<Kokkos::HostSpace, void>::~SharedAllocationRecord()
#if defined( \
    KOKKOS_IMPL_INTEL_WORKAROUND_NOEXCEPT_SPECIFICATION_VIRTUAL_FUNCTION)
    noexcept
#endif
{
#if defined(KOKKOS_ENABLE_PROFILING)
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    Kokkos::Profiling::deallocateData(
        Kokkos::Profiling::SpaceHandle(Kokkos::HostSpace::name()),
        RecordBase::m_alloc_ptr->m_label, data(), size());
  }
#endif

  m_space.deallocate(SharedAllocationRecord<void, void>::m_alloc_ptr,
                     SharedAllocationRecord<void, void>::m_alloc_size);
}

SharedAllocationHeader *_do_allocation(Kokkos::HostSpace const &space,
                                       std::string const &label,
                                       size_t alloc_size) {
  try {
    return reinterpret_cast<SharedAllocationHeader *>(
        space.allocate(alloc_size));
  } catch (Experimental::RawMemoryAllocationFailure const &failure) {
    if (failure.failure_mode() == Experimental::RawMemoryAllocationFailure::
                                      FailureMode::AllocationNotAligned) {
      // TODO: delete the misaligned memory
    }

    std::cerr << "Kokkos failed to allocate memory for label \"" << label
              << "\".  Allocation using MemorySpace named \"" << space.name()
              << " failed with the following error:  ";
    failure.print_error_message(std::cerr);
    std::cerr.flush();
    Kokkos::Impl::throw_runtime_exception("Memory allocation failure");
  }
  return nullptr;  // unreachable
}

SharedAllocationRecord<Kokkos::HostSpace, void>::SharedAllocationRecord(
    const Kokkos::HostSpace &arg_space, const std::string &arg_label,
    const size_t arg_alloc_size,
    const SharedAllocationRecord<void, void>::function_type arg_dealloc)
    // Pass through allocated [ SharedAllocationHeader , user_memory ]
    // Pass through deallocation function
    : SharedAllocationRecord<void, void>(
#ifdef KOKKOS_DEBUG
          &SharedAllocationRecord<Kokkos::HostSpace, void>::s_root_record,
#endif
          Impl::checked_allocation_with_header(arg_space, arg_label,
                                               arg_alloc_size),
          sizeof(SharedAllocationHeader) + arg_alloc_size, arg_dealloc),
      m_space(arg_space) {
#if defined(KOKKOS_ENABLE_PROFILING)
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    Kokkos::Profiling::allocateData(
        Kokkos::Profiling::SpaceHandle(arg_space.name()), arg_label, data(),
        arg_alloc_size);
  }
#endif
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

void *SharedAllocationRecord<Kokkos::HostSpace, void>::allocate_tracked(
    const Kokkos::HostSpace &arg_space, const std::string &arg_alloc_label,
    const size_t arg_alloc_size) {
  if (!arg_alloc_size) return nullptr;

  SharedAllocationRecord *const r =
      allocate(arg_space, arg_alloc_label, arg_alloc_size);

  RecordBase::increment(r);

  return r->data();
}

void SharedAllocationRecord<Kokkos::HostSpace, void>::deallocate_tracked(
    void *const arg_alloc_ptr) {
  if (arg_alloc_ptr != nullptr) {
    SharedAllocationRecord *const r = get_record(arg_alloc_ptr);

    RecordBase::decrement(r);
  }
}

void *SharedAllocationRecord<Kokkos::HostSpace, void>::reallocate_tracked(
    void *const arg_alloc_ptr, const size_t arg_alloc_size) {
  SharedAllocationRecord *const r_old = get_record(arg_alloc_ptr);
  SharedAllocationRecord *const r_new =
      allocate(r_old->m_space, r_old->get_label(), arg_alloc_size);

  Kokkos::Impl::DeepCopy<HostSpace, HostSpace>(
      r_new->data(), r_old->data(), std::min(r_old->size(), r_new->size()));

  RecordBase::increment(r_new);
  RecordBase::decrement(r_old);

  return r_new->data();
}

SharedAllocationRecord<Kokkos::HostSpace, void> *
SharedAllocationRecord<Kokkos::HostSpace, void>::get_record(void *alloc_ptr) {
  typedef SharedAllocationHeader Header;
  typedef SharedAllocationRecord<Kokkos::HostSpace, void> RecordHost;

  SharedAllocationHeader const *const head =
      alloc_ptr ? Header::get_header(alloc_ptr) : nullptr;
  RecordHost *const record =
      head ? static_cast<RecordHost *>(head->m_record) : nullptr;

  if (!alloc_ptr || record->m_alloc_ptr != head) {
    Kokkos::Impl::throw_runtime_exception(
        std::string("Kokkos::Impl::SharedAllocationRecord< Kokkos::HostSpace , "
                    "void >::get_record ERROR"));
  }

  return record;
}

// Iterate records to print orphaned memory ...
#ifdef KOKKOS_DEBUG
void SharedAllocationRecord<Kokkos::HostSpace, void>::print_records(
    std::ostream &s, const Kokkos::HostSpace &, bool detail) {
  SharedAllocationRecord<void, void>::print_host_accessible_records(
      s, "HostSpace", &s_root_record, detail);
}
#else
void SharedAllocationRecord<Kokkos::HostSpace, void>::print_records(
    std::ostream &, const Kokkos::HostSpace &, bool) {
  throw_runtime_exception(
      "SharedAllocationRecord<HostSpace>::print_records only works with "
      "KOKKOS_DEBUG enabled");
}
#endif

}  // namespace Impl
}  // namespace Kokkos

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace Kokkos {
namespace {
const unsigned HOST_SPACE_ATOMIC_MASK     = 0xFFFF;
const unsigned HOST_SPACE_ATOMIC_XOR_MASK = 0x5A39;
static int HOST_SPACE_ATOMIC_LOCKS[HOST_SPACE_ATOMIC_MASK + 1];
}  // namespace

namespace Impl {
void init_lock_array_host_space() {
  static int is_initialized = 0;
  if (!is_initialized)
    for (int i = 0; i < static_cast<int>(HOST_SPACE_ATOMIC_MASK + 1); i++)
      HOST_SPACE_ATOMIC_LOCKS[i] = 0;
}

bool lock_address_host_space(void *ptr) {
#if defined(KOKKOS_ENABLE_ISA_X86_64) && defined(KOKKOS_ENABLE_TM) && \
    !defined(KOKKOS_COMPILER_PGI)
  const unsigned status = _xbegin();

  if (_XBEGIN_STARTED == status) {
    const int val =
        HOST_SPACE_ATOMIC_LOCKS[((size_t(ptr) >> 2) & HOST_SPACE_ATOMIC_MASK) ^
                                HOST_SPACE_ATOMIC_XOR_MASK];

    if (0 == val) {
      HOST_SPACE_ATOMIC_LOCKS[((size_t(ptr) >> 2) & HOST_SPACE_ATOMIC_MASK) ^
                              HOST_SPACE_ATOMIC_XOR_MASK] = 1;
    } else {
      _xabort(1);
    }

    _xend();

    return 1;
  } else {
#endif
    return 0 == atomic_compare_exchange(
                    &HOST_SPACE_ATOMIC_LOCKS[((size_t(ptr) >> 2) &
                                              HOST_SPACE_ATOMIC_MASK) ^
                                             HOST_SPACE_ATOMIC_XOR_MASK],
                    0, 1);
#if defined(KOKKOS_ENABLE_ISA_X86_64) && defined(KOKKOS_ENABLE_TM) && \
    !defined(KOKKOS_COMPILER_PGI)
  }
#endif
}

void unlock_address_host_space(void *ptr) {
#if defined(KOKKOS_ENABLE_ISA_X86_64) && defined(KOKKOS_ENABLE_TM) && \
    !defined(KOKKOS_COMPILER_PGI)
  const unsigned status = _xbegin();

  if (_XBEGIN_STARTED == status) {
    HOST_SPACE_ATOMIC_LOCKS[((size_t(ptr) >> 2) & HOST_SPACE_ATOMIC_MASK) ^
                            HOST_SPACE_ATOMIC_XOR_MASK] = 0;
  } else {
#endif
    atomic_exchange(
        &HOST_SPACE_ATOMIC_LOCKS[((size_t(ptr) >> 2) & HOST_SPACE_ATOMIC_MASK) ^
                                 HOST_SPACE_ATOMIC_XOR_MASK],
        0);
#if defined(KOKKOS_ENABLE_ISA_X86_64) && defined(KOKKOS_ENABLE_TM) && \
    !defined(KOKKOS_COMPILER_PGI)
  }
#endif
}

}  // namespace Impl
}  // namespace Kokkos
