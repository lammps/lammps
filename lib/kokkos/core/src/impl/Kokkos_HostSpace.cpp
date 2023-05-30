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

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#endif

#include <Kokkos_Macros.hpp>

#include <impl/Kokkos_Error.hpp>
#include <impl/Kokkos_MemorySpace.hpp>
#include <impl/Kokkos_Tools.hpp>

/*--------------------------------------------------------------------------*/

#if defined(KOKKOS_COMPILER_INTEL) && !defined(KOKKOS_ENABLE_CUDA)

// Intel specialized allocator does not interoperate with CUDA memory allocation

#define KOKKOS_ENABLE_INTEL_MM_ALLOC

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

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

/* Default allocation mechanism */
HostSpace::HostSpace()
    : m_alloc_mech(
#if defined(KOKKOS_ENABLE_INTEL_MM_ALLOC)
          HostSpace::INTEL_MM_ALLOC
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
#endif
  else {
    const char *const mech =
        (arg_alloc_mech == HostSpace::INTEL_MM_ALLOC)
            ? "INTEL_MM_ALLOC"
            : ((arg_alloc_mech == HostSpace::POSIX_MMAP) ? "POSIX_MMAP" : "");

    std::string msg;
    msg.append("Kokkos::HostSpace ");
    msg.append(mech);
    msg.append(" is not available");
    Kokkos::Impl::throw_runtime_exception(msg);
  }
}

void *HostSpace::allocate(const size_t arg_alloc_size) const {
  return allocate("[unlabeled]", arg_alloc_size);
}
void *HostSpace::allocate(const char *arg_label, const size_t arg_alloc_size,
                          const size_t

                              arg_logical_size) const {
  return impl_allocate(arg_label, arg_alloc_size, arg_logical_size);
}
void *HostSpace::impl_allocate(
    const char *arg_label, const size_t arg_alloc_size,
    const size_t arg_logical_size,
    const Kokkos::Tools::SpaceHandle arg_handle) const {
  const size_t reported_size =
      (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;
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
  if (Kokkos::Profiling::profileLibraryLoaded()) {
    Kokkos::Profiling::allocateData(arg_handle, arg_label, ptr, reported_size);
  }
  return ptr;
}

void HostSpace::deallocate(void *const arg_alloc_ptr,
                           const size_t arg_alloc_size) const {
  deallocate("[unlabeled]", arg_alloc_ptr, arg_alloc_size);
}

void HostSpace::deallocate(const char *arg_label, void *const arg_alloc_ptr,
                           const size_t arg_alloc_size,
                           const size_t

                               arg_logical_size) const {
  impl_deallocate(arg_label, arg_alloc_ptr, arg_alloc_size, arg_logical_size);
}
void HostSpace::impl_deallocate(
    const char *arg_label, void *const arg_alloc_ptr,
    const size_t arg_alloc_size, const size_t arg_logical_size,
    const Kokkos::Tools::SpaceHandle arg_handle) const {
  if (arg_alloc_ptr) {
    Kokkos::fence("HostSpace::impl_deallocate before free");
    size_t reported_size =
        (arg_logical_size > 0) ? arg_logical_size : arg_alloc_size;
    if (Kokkos::Profiling::profileLibraryLoaded()) {
      Kokkos::Profiling::deallocateData(arg_handle, arg_label, arg_alloc_ptr,
                                        reported_size);
    }
    if (m_alloc_mech == STD_MALLOC) {
      void *alloc_ptr = *(reinterpret_cast<void **>(arg_alloc_ptr) - 1);
      free(alloc_ptr);
    }
#if defined(KOKKOS_ENABLE_INTEL_MM_ALLOC)
    else if (m_alloc_mech == INTEL_MM_ALLOC) {
      _mm_free(arg_alloc_ptr);
    }
#endif
  }
}

}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

#ifdef KOKKOS_ENABLE_DEBUG
SharedAllocationRecord<void, void>
    SharedAllocationRecord<Kokkos::HostSpace, void>::s_root_record;
#endif

SharedAllocationRecord<Kokkos::HostSpace, void>::~SharedAllocationRecord()
#if defined( \
    KOKKOS_IMPL_INTEL_WORKAROUND_NOEXCEPT_SPECIFICATION_VIRTUAL_FUNCTION)
    noexcept
#endif
{
  m_space.deallocate(m_label.c_str(),
                     SharedAllocationRecord<void, void>::m_alloc_ptr,
                     SharedAllocationRecord<void, void>::m_alloc_size,
                     (SharedAllocationRecord<void, void>::m_alloc_size -
                      sizeof(SharedAllocationHeader)));
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
    : base_t(
#ifdef KOKKOS_ENABLE_DEBUG
          &SharedAllocationRecord<Kokkos::HostSpace, void>::s_root_record,
#endif
          Impl::checked_allocation_with_header(arg_space, arg_label,
                                               arg_alloc_size),
          sizeof(SharedAllocationHeader) + arg_alloc_size, arg_dealloc,
          arg_label),
      m_space(arg_space) {
  this->base_t::_fill_host_accessible_header_info(*RecordBase::m_alloc_ptr,
                                                  arg_label);
}

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
  return 0 == atomic_compare_exchange(
                  &HOST_SPACE_ATOMIC_LOCKS[((size_t(ptr) >> 2) &
                                            HOST_SPACE_ATOMIC_MASK) ^
                                           HOST_SPACE_ATOMIC_XOR_MASK],
                  0, 1);
}

void unlock_address_host_space(void *ptr) {
  atomic_exchange(
      &HOST_SPACE_ATOMIC_LOCKS[((size_t(ptr) >> 2) & HOST_SPACE_ATOMIC_MASK) ^
                               HOST_SPACE_ATOMIC_XOR_MASK],
      0);
}

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
template class SharedAllocationRecordCommon<Kokkos::HostSpace>;

}  // end namespace Impl
}  // end namespace Kokkos

// </editor-fold> end Explicit instantiations of CRTP Base classes }}}1
//==============================================================================
