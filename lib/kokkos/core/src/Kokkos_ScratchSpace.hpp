//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#include <Kokkos_Macros.hpp>
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
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
  KOKKOS_INLINE_FUNCTION void* get_shmem_common(
      const IntType& size, [[maybe_unused]] const ptrdiff_t alignment,
      int level = -1) const {
    if (level == -1) level = m_default_level;
    auto& m_iter    = (level == 0) ? m_iter_L0 : m_iter_L1;
    auto m_iter_old = m_iter;
    if constexpr (alignment_requested) {
      const ptrdiff_t missalign = size_t(m_iter) % alignment;
      if (missalign) m_iter += alignment - missalign;
    }

    // This is each thread's start pointer for its allocation
    // Note: for team scratch m_offset is 0, since every
    // thread will get back the same shared pointer
    void* tmp           = m_iter + m_offset * size;
    uintptr_t increment = size * m_multiplier;

    // Cast to uintptr_t to avoid problems with pointer arithmetic using SYCL
    const auto end_iter =
        reinterpret_cast<uintptr_t>((level == 0) ? m_end_L0 : m_end_L1);
    auto current_iter = reinterpret_cast<uintptr_t>(m_iter);
    auto capacity     = end_iter - current_iter;

    if (increment > capacity) {
      // Request did overflow: return nullptr and reset m_iter
      m_iter = m_iter_old;
      tmp    = nullptr;
#ifdef KOKKOS_ENABLE_DEBUG
      // mfh 23 Jun 2015: printf call consumes 25 registers
      // in a CUDA build, so only print in debug mode.  The
      // function still returns nullptr if not enough memory.
      Kokkos::printf(
          "ScratchMemorySpace<...>::get_shmem: Failed to allocate "
          "%ld byte(s); remaining capacity is %ld byte(s)\n",
          long(size), long(capacity));
#endif  // KOKKOS_ENABLE_DEBUG
    } else {
      m_iter += increment;
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
