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
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#endif

#include <impl/Kokkos_Error.hpp>

#include <ostream>
#include <sstream>
#include <cstdint>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/* Verify size constraints:
 *   min_block_alloc_size <= max_block_alloc_size
 *   max_block_alloc_size <= min_superblock_size
 *   min_superblock_size  <= max_superblock_size
 *   min_superblock_size  <= min_total_alloc_size
 *   min_superblock_size  <= min_block_alloc_size *
 *                           max_block_per_superblock
 */
void memory_pool_bounds_verification(size_t min_block_alloc_size,
                                     size_t max_block_alloc_size,
                                     size_t min_superblock_size,
                                     size_t max_superblock_size,
                                     size_t max_block_per_superblock,
                                     size_t min_total_alloc_size) {
  const size_t max_superblock = min_block_alloc_size * max_block_per_superblock;

  if ((size_t(max_superblock_size) < min_superblock_size) ||
      (min_total_alloc_size < min_superblock_size) ||
      (max_superblock < min_superblock_size) ||
      (min_superblock_size < max_block_alloc_size) ||
      (max_block_alloc_size < min_block_alloc_size)) {
    std::ostringstream msg;

    msg << "Kokkos::MemoryPool size constraint violation";

    if (size_t(max_superblock_size) < min_superblock_size) {
      msg << " : max_superblock_size(" << max_superblock_size
          << ") < min_superblock_size(" << min_superblock_size << ")";
    }

    if (min_total_alloc_size < min_superblock_size) {
      msg << " : min_total_alloc_size(" << min_total_alloc_size
          << ") < min_superblock_size(" << min_superblock_size << ")";
    }

    if (max_superblock < min_superblock_size) {
      msg << " : max_superblock(" << max_superblock
          << ") < min_superblock_size(" << min_superblock_size << ")";
    }

    if (min_superblock_size < max_block_alloc_size) {
      msg << " : min_superblock_size(" << min_superblock_size
          << ") < max_block_alloc_size(" << max_block_alloc_size << ")";
    }

    if (max_block_alloc_size < min_block_alloc_size) {
      msg << " : max_block_alloc_size(" << max_block_alloc_size
          << ") < min_block_alloc_size(" << min_block_alloc_size << ")";
    }

    Kokkos::Impl::throw_runtime_exception(msg.str());
  }
}

// This has way too many parameters, but it is entirely for moving the iostream
// inclusion out of the header file with as few changes as possible
void _print_memory_pool_state(std::ostream& s, uint32_t const* sb_state_ptr,
                              int32_t sb_count, uint32_t sb_size_lg2,
                              uint32_t sb_state_size, uint32_t state_shift,
                              uint32_t state_used_mask) {
  s << "pool_size(" << (size_t(sb_count) << sb_size_lg2) << ")"
    << " superblock_size(" << (1LU << sb_size_lg2) << ")" << std::endl;

  for (int32_t i = 0; i < sb_count; ++i, sb_state_ptr += sb_state_size) {
    if (*sb_state_ptr) {
      const uint32_t block_count_lg2 = (*sb_state_ptr) >> state_shift;
      const uint32_t block_size_lg2  = sb_size_lg2 - block_count_lg2;
      const uint32_t block_count     = 1u << block_count_lg2;
      const uint32_t block_used      = (*sb_state_ptr) & state_used_mask;

      s << "Superblock[ " << i << " / " << sb_count << " ] {"
        << " block_size(" << (1 << block_size_lg2) << ")"
        << " block_count( " << block_used << " / " << block_count << " )"
        << std::endl;
    }
  }
}

}  // namespace Impl
}  // namespace Kokkos
