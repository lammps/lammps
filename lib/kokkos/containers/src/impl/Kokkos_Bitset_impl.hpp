// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_BITSET_IMPL_HPP
#define KOKKOS_BITSET_IMPL_HPP

#include <Kokkos_Macros.hpp>
#include <Kokkos_BitManipulation.hpp>
#include <cstdint>

#include <cstdio>
#include <climits>
#include <iomanip>

namespace Kokkos {
namespace Impl {

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
    count += Kokkos::Experimental::popcount_builtin(m_bitset.m_blocks[i]);
  }
};

}  // namespace Impl
}  // namespace Kokkos

#endif  // KOKKOS_BITSET_IMPL_HPP
