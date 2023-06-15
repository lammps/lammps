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
#include <Kokkos_Macros.hpp>
#ifndef KOKKOS_ENABLE_DEPRECATED_CODE_3
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#else
KOKKOS_IMPL_WARNING("Including non-public Kokkos header files is not allowed.")
#endif
#endif
#ifndef KOKKOS_SCRATCHSPACE_HPP
#define KOKKOS_SCRATCHSPACE_HPP

#include <cstdio>
#include <cstddef>
#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_Concepts.hpp>

/*--------------------------------------------------------------------------*/

namespace Kokkos {

/** \brief  Scratch memory space associated with an execution space.
 *
 */
template <class ExecSpace>
class ScratchMemorySpace {
  static_assert(
      is_execution_space<ExecSpace>::value,
      "Instantiating ScratchMemorySpace on non-execution-space type.");

 public:
  // Minimal overalignment used by view scratch allocations
  constexpr static int ALIGN = 8;

 private:
  mutable char* m_iter_L0 = nullptr;
  mutable char* m_iter_L1 = nullptr;
  char* m_end_L0          = nullptr;
  char* m_end_L1          = nullptr;

  mutable int m_multiplier    = 0;
  mutable int m_offset        = 0;
  mutable int m_default_level = 0;

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  constexpr static int DEFAULT_ALIGNMENT_MASK = ALIGN - 1;
#endif

 public:
  //! Tag this class as a memory space
  using memory_space    = ScratchMemorySpace<ExecSpace>;
  using execution_space = ExecSpace;
  //! This execution space preferred device_type
  using device_type = Kokkos::Device<execution_space, memory_space>;

  using array_layout = typename ExecSpace::array_layout;
  using size_type    = typename ExecSpace::size_type;

  static constexpr const char* name() { return "ScratchMemorySpace"; }

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  // This function is unused
  template <typename IntType>
  KOKKOS_DEPRECATED KOKKOS_INLINE_FUNCTION static constexpr IntType align(
      const IntType& size) {
    return (size + DEFAULT_ALIGNMENT_MASK) & ~DEFAULT_ALIGNMENT_MASK;
  }
#endif

  template <typename IntType>
  KOKKOS_INLINE_FUNCTION void* get_shmem(const IntType& size,
                                         int level = -1) const {
    return get_shmem_common</*alignment_requested*/ false>(size, 1, level);
  }

  template <typename IntType>
  KOKKOS_INLINE_FUNCTION void* get_shmem_aligned(const IntType& size,
                                                 const ptrdiff_t alignment,
                                                 int level = -1) const {
    return get_shmem_common</*alignment_requested*/ true>(size, alignment,
                                                          level);
  }

 private:
  template <bool alignment_requested, typename IntType>
  KOKKOS_INLINE_FUNCTION void* get_shmem_common(const IntType& size,
                                                const ptrdiff_t alignment,
                                                int level = -1) const {
    if (level == -1) level = m_default_level;
    auto& m_iter = (level == 0) ? m_iter_L0 : m_iter_L1;
    auto& m_end  = (level == 0) ? m_end_L0 : m_end_L1;

    if (alignment_requested) {
      const ptrdiff_t missalign = size_t(m_iter) % alignment;
      if (missalign) m_iter += alignment - missalign;
    }

    // This is each thread's start pointer for its allocation
    // Note: for team scratch m_offset is 0, since every
    // thread will get back the same shared pointer
    void* tmp           = m_iter + m_offset * size;
    ptrdiff_t increment = size * m_multiplier;

    // increment m_iter first and decrement it again if not
    // enough memory was available. In the non-failing path
    // this will save instructions.
    m_iter += increment;

    if (m_end < m_iter) {
      // Request did overflow: reset the base team ptr, and
      // return nullptr
      m_iter -= increment;
      tmp = nullptr;
#ifdef KOKKOS_ENABLE_DEBUG
      // mfh 23 Jun 2015: printf call consumes 25 registers
      // in a CUDA build, so only print in debug mode.  The
      // function still returns nullptr if not enough memory.
      KOKKOS_IMPL_DO_NOT_USE_PRINTF(
          "ScratchMemorySpace<...>::get_shmem: Failed to allocate "
          "%ld byte(s); remaining capacity is %ld byte(s)\n",
          long(size), long(m_end - m_iter));
#endif  // KOKKOS_ENABLE_DEBUG
    }
    return tmp;
  }

 public:
  KOKKOS_DEFAULTED_FUNCTION
  ScratchMemorySpace() = default;

  template <typename IntType>
  KOKKOS_INLINE_FUNCTION ScratchMemorySpace(void* ptr_L0,
                                            const IntType& size_L0,
                                            void* ptr_L1           = nullptr,
                                            const IntType& size_L1 = 0)
      : m_iter_L0(static_cast<char*>(ptr_L0)),
        m_iter_L1(static_cast<char*>(ptr_L1)),
        m_end_L0(static_cast<char*>(ptr_L0) + size_L0),
        m_end_L1(static_cast<char*>(ptr_L1) + size_L1),
        m_multiplier(1),
        m_offset(0),
        m_default_level(0) {}

  KOKKOS_INLINE_FUNCTION
  const ScratchMemorySpace& set_team_thread_mode(const int& level,
                                                 const int& multiplier,
                                                 const int& offset) const {
    m_default_level = level;
    m_multiplier    = multiplier;
    m_offset        = offset;
    return *this;
  }
};

}  // namespace Kokkos

#endif /* #ifndef KOKKOS_SCRATCHSPACE_HPP */
