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

#ifndef KOKKOS_BITSET_IMPL_HPP
#define KOKKOS_BITSET_IMPL_HPP

#include <Kokkos_Macros.hpp>
#include <impl/Kokkos_BitOps.hpp>
#include <cstdint>

#include <cstdio>
#include <climits>
#include <iostream>
#include <iomanip>

namespace Kokkos {
namespace Impl {

KOKKOS_FORCEINLINE_FUNCTION
unsigned rotate_right(unsigned i, int r) {
  constexpr int size = static_cast<int>(sizeof(unsigned) * CHAR_BIT);
  return r ? ((i >> r) | (i << (size - r))) : i;
}

template <typename Bitset>
struct BitsetCount {
  using bitset_type = Bitset;
  using execution_space =
      typename bitset_type::execution_space::execution_space;
  using size_type  = typename bitset_type::size_type;
  using value_type = size_type;

  bitset_type m_bitset;

  BitsetCount(bitset_type const& bitset) : m_bitset(bitset) {}

  size_type apply() const {
    size_type count = 0u;
    parallel_reduce("Kokkos::Impl::BitsetCount::apply",
                    m_bitset.m_blocks.extent(0), *this, count);
    return count;
  }

  KOKKOS_INLINE_FUNCTION
  void init(value_type& count) const { count = 0u; }

  KOKKOS_INLINE_FUNCTION
  void join(value_type& count, const size_type& incr) const { count += incr; }

  KOKKOS_INLINE_FUNCTION
  void operator()(size_type i, value_type& count) const {
    count += bit_count(m_bitset.m_blocks[i]);
  }
};

}  // namespace Impl
}  // namespace Kokkos

#endif  // KOKKOS_BITSET_IMPL_HPP
