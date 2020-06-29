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

#include <ostream>
#include <sstream>
#include <impl/Kokkos_Error.hpp>

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

}  // namespace Impl
}  // namespace Kokkos
