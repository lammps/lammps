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

#ifndef KOKKOS_CONCURRENTBITSET_HPP
#define KOKKOS_CONCURRENTBITSET_HPP

#include <stdint.h>
#include <Kokkos_Atomic.hpp>
#include <impl/Kokkos_BitOps.hpp>
#include <impl/Kokkos_ClockTic.hpp>

namespace Kokkos {
namespace Impl {

struct concurrent_bitset {
 public:
  // 32 bits per integer value

  enum : uint32_t { bits_per_int_lg2 = 5 };
  enum : uint32_t { bits_per_int_mask = (1 << bits_per_int_lg2) - 1 };

  // Buffer is uint32_t[ buffer_bound ]
  //   [ uint32_t { state_header | used_count } , uint32_t bits[*] ]
  //
  //  Maximum bit count is 33 million (1u<<25):
  //
  //  - Maximum bit set size occupies 1 Mbyte
  //
  //  - State header can occupy bits [30-26]
  //    which can be the bit_count_lg2
  //
  //  - Accept at least 33 million concurrent calls to 'acquire'
  //    before risking an overflow race condition on a full bitset.

  enum : uint32_t { max_bit_count_lg2 = 25 };
  enum : uint32_t { max_bit_count = 1u << max_bit_count_lg2 };
  enum : uint32_t { state_shift = 26 };
  enum : uint32_t { state_used_mask = (1 << state_shift) - 1 };
  enum : uint32_t { state_header_mask = uint32_t(0x001f) << state_shift };

  KOKKOS_INLINE_FUNCTION static constexpr uint32_t buffer_bound_lg2(
      uint32_t const bit_bound_lg2) noexcept {
    return bit_bound_lg2 <= max_bit_count_lg2
               ? 1 + (1u << (bit_bound_lg2 > bits_per_int_lg2
                                 ? bit_bound_lg2 - bits_per_int_lg2
                                 : 0))
               : 0;
  }

  /**\brief  Initialize bitset buffer */
  KOKKOS_INLINE_FUNCTION static constexpr uint32_t buffer_bound(
      uint32_t const bit_bound) noexcept {
    return bit_bound <= max_bit_count
               ? 1 + (bit_bound >> bits_per_int_lg2) +
                     (bit_bound & bits_per_int_mask ? 1 : 0)
               : 0;
  }

  /**\brief  Claim any bit within the bitset bound.
   *
   *  Return : ( which_bit , bit_count )
   *
   *  if success then
   *    bit_count is the atomic-count of claimed > 0
   *    which_bit is the claimed bit >= 0
   *  else if attempt failed due to filled buffer
   *    bit_count == which_bit == -1
   *  else if attempt failed due to non-matching state_header
   *    bit_count == which_bit == -2
   *  else if attempt failed due to max_bit_count_lg2 < bit_bound_lg2
   *                             or invalid state_header
   *                             or (1u << bit_bound_lg2) <= bit
   *    bit_count == which_bit == -3
   *  endif
   *
   *  Recommended to have hint
   *    bit = Kokkos::Impl::clock_tic() & ((1u<<bit_bound_lg2) - 1)
   */
  KOKKOS_INLINE_FUNCTION static Kokkos::pair<int, int> acquire_bounded_lg2(
      uint32_t volatile *const buffer, uint32_t const bit_bound_lg2,
      uint32_t bit = 0 /* optional hint */
      ,
      uint32_t const state_header = 0 /* optional header */
      ) noexcept {
    using type = Kokkos::pair<int, int>;

    const uint32_t bit_bound  = 1 << bit_bound_lg2;
    const uint32_t word_count = bit_bound >> bits_per_int_lg2;

    if ((max_bit_count_lg2 < bit_bound_lg2) ||
        (state_header & ~state_header_mask) || (bit_bound < bit)) {
      return type(-3, -3);
    }

    // Use potentially two fetch_add to avoid CAS loop.
    // Could generate "racing" failure-to-acquire
    // when is full at the atomic_fetch_add(+1)
    // then a release occurs before the atomic_fetch_add(-1).

    const uint32_t state = (uint32_t)Kokkos::atomic_fetch_add(
        reinterpret_cast<volatile int *>(buffer), 1);

    const uint32_t state_error = state_header != (state & state_header_mask);

    const uint32_t state_bit_used = state & state_used_mask;

    if (state_error || (bit_bound <= state_bit_used)) {
      Kokkos::atomic_fetch_add(reinterpret_cast<volatile int *>(buffer), -1);
      return state_error ? type(-2, -2) : type(-1, -1);
    }

    // Do not update bit until count is visible:

    Kokkos::memory_fence();

    // There is a zero bit available somewhere,
    // now find the (first) available bit and set it.

    while (1) {
      const uint32_t word = bit >> bits_per_int_lg2;
      const uint32_t mask = 1u << (bit & bits_per_int_mask);
      const uint32_t prev = Kokkos::atomic_fetch_or(buffer + word + 1, mask);

      if (!(prev & mask)) {
        // Successfully claimed 'result.first' by
        // atomically setting that bit.
        return type(bit, state_bit_used + 1);
      }

      // Failed race to set the selected bit
      // Find a new bit to try.

      const int j = Kokkos::Impl::bit_first_zero(prev);

      if (0 <= j) {
        bit = (word << bits_per_int_lg2) | uint32_t(j);
      } else {
        bit = ((word + 1) < word_count ? ((word + 1) << bits_per_int_lg2) : 0) |
              (bit & bits_per_int_mask);
      }
    }
  }

  /**\brief  Claim any bit within the bitset bound.
   *
   *  Return : ( which_bit , bit_count )
   *
   *  if success then
   *    bit_count is the atomic-count of claimed > 0
   *    which_bit is the claimed bit >= 0
   *  else if attempt failed due to filled buffer
   *    bit_count == which_bit == -1
   *  else if attempt failed due to non-matching state_header
   *    bit_count == which_bit == -2
   *  else if attempt failed due to max_bit_count_lg2 < bit_bound_lg2
   *                             or invalid state_header
   *                             or bit_bound <= bit
   *    bit_count == which_bit == -3
   *  endif
   *
   *  Recommended to have hint
   *    bit = Kokkos::Impl::clock_tic() % bit_bound
   */
  KOKKOS_INLINE_FUNCTION static Kokkos::pair<int, int> acquire_bounded(
      uint32_t volatile *const buffer, uint32_t const bit_bound,
      uint32_t bit = 0 /* optional hint */
      ,
      uint32_t const state_header = 0 /* optional header */
      ) noexcept {
    using type = Kokkos::pair<int, int>;

    if ((max_bit_count < bit_bound) || (state_header & ~state_header_mask) ||
        (bit_bound <= bit)) {
      return type(-3, -3);
    }

    const uint32_t word_count = bit_bound >> bits_per_int_lg2;

    // Use potentially two fetch_add to avoid CAS loop.
    // Could generate "racing" failure-to-acquire
    // when is full at the atomic_fetch_add(+1)
    // then a release occurs before the atomic_fetch_add(-1).

    const uint32_t state = (uint32_t)Kokkos::atomic_fetch_add(
        reinterpret_cast<volatile int *>(buffer), 1);

    const uint32_t state_error = state_header != (state & state_header_mask);

    const uint32_t state_bit_used = state & state_used_mask;

    if (state_error || (bit_bound <= state_bit_used)) {
      Kokkos::atomic_fetch_add(reinterpret_cast<volatile int *>(buffer), -1);
      return state_error ? type(-2, -2) : type(-1, -1);
    }

    // Do not update bit until count is visible:

    Kokkos::memory_fence();

    // There is a zero bit available somewhere,
    // now find the (first) available bit and set it.

    while (1) {
      const uint32_t word = bit >> bits_per_int_lg2;
      const uint32_t mask = 1u << (bit & bits_per_int_mask);
      const uint32_t prev = Kokkos::atomic_fetch_or(buffer + word + 1, mask);

      if (!(prev & mask)) {
        // Successfully claimed 'result.first' by
        // atomically setting that bit.
        // Flush the set operation. Technically this only needs to be acquire/
        // release semantics and not sequentially consistent, but for now
        // we'll just do this.
        Kokkos::memory_fence();
        return type(bit, state_bit_used + 1);
      }

      // Failed race to set the selected bit
      // Find a new bit to try.

      const int j = Kokkos::Impl::bit_first_zero(prev);

      if (0 <= j) {
        bit = (word << bits_per_int_lg2) | uint32_t(j);
      }

      if ((j < 0) || (bit_bound <= bit)) {
        bit = ((word + 1) < word_count ? ((word + 1) << bits_per_int_lg2) : 0) |
              (bit & bits_per_int_mask);
      }
    }
  }

  /**\brief
   *
   *  Requires: 'bit' previously acquired and has not yet been released.
   *
   *  Returns:
   *    0 <= used count after successful release
   *    -1 bit was already released
   *    -2 state_header error
   */
  KOKKOS_INLINE_FUNCTION static int release(
      uint32_t volatile *const buffer, uint32_t const bit,
      uint32_t const state_header = 0 /* optional header */
      ) noexcept {
    if (state_header != (state_header_mask & *buffer)) {
      return -2;
    }

    const uint32_t mask = 1u << (bit & bits_per_int_mask);
    const uint32_t prev =
        Kokkos::atomic_fetch_and(buffer + (bit >> bits_per_int_lg2) + 1, ~mask);

    if (!(prev & mask)) {
      return -1;
    }

    // Do not update count until bit clear is visible
    Kokkos::memory_fence();

    const int count =
        Kokkos::atomic_fetch_add(reinterpret_cast<volatile int *>(buffer), -1);

    // Flush the store-release
    Kokkos::memory_fence();

    return (count & state_used_mask) - 1;
  }

  /**\brief
   *
   *  Requires: Bit within bounds and not already set.
   *
   *  Returns:
   *    0 <= used count after successful release
   *    -1 bit was already released
   *    -2 bit or state_header error
   */
  KOKKOS_INLINE_FUNCTION static int set(
      uint32_t volatile *const buffer, uint32_t const bit,
      uint32_t const state_header = 0 /* optional header */
      ) noexcept {
    if (state_header != (state_header_mask & *buffer)) {
      return -2;
    }

    const uint32_t mask = 1u << (bit & bits_per_int_mask);
    const uint32_t prev =
        Kokkos::atomic_fetch_or(buffer + (bit >> bits_per_int_lg2) + 1, mask);

    if (!(prev & mask)) {
      return -1;
    }

    // Do not update count until bit clear is visible
    Kokkos::memory_fence();

    const int count =
        Kokkos::atomic_fetch_add(reinterpret_cast<volatile int *>(buffer), -1);

    return (count & state_used_mask) - 1;
  }
};

}  // namespace Impl
}  // namespace Kokkos

#endif /* #ifndef KOKKOS_CONCURRENTBITSET_HPP */
