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

#ifndef KOKKOS_UNORDERED_MAP_IMPL_HPP
#define KOKKOS_UNORDERED_MAP_IMPL_HPP

#include <Kokkos_Core.hpp>
#include <cstdint>

#include <cstdio>
#include <climits>
#include <iostream>
#include <iomanip>

namespace Kokkos {
namespace Impl {

uint32_t find_hash_size(uint32_t size);

template <typename Map>
struct UnorderedMapRehash {
  using map_type        = Map;
  using const_map_type  = typename map_type::const_map_type;
  using execution_space = typename map_type::execution_space;
  using size_type       = typename map_type::size_type;

  map_type m_dst;
  const_map_type m_src;

  UnorderedMapRehash(map_type const& dst, const_map_type const& src)
      : m_dst(dst), m_src(src) {}

  void apply() const {
    parallel_for("Kokkos::Impl::UnorderedMapRehash::apply", m_src.capacity(),
                 *this);
  }

  template <typename Dummy = typename map_type::value_type>
  KOKKOS_INLINE_FUNCTION std::enable_if_t<std::is_void<Dummy>::value>
  operator()(size_type i) const {
    if (m_src.valid_at(i)) m_dst.insert(m_src.key_at(i));
  }

  template <typename Dummy = typename map_type::value_type>
  KOKKOS_INLINE_FUNCTION std::enable_if_t<!std::is_void<Dummy>::value>
  operator()(size_type i) const {
    if (m_src.valid_at(i)) m_dst.insert(m_src.key_at(i), m_src.value_at(i));
  }
};

template <typename UMap>
struct UnorderedMapErase {
  using map_type        = UMap;
  using execution_space = typename map_type::execution_space;
  using size_type       = typename map_type::size_type;
  using key_type        = typename map_type::key_type;
  using value_type      = typename map_type::impl_value_type;

  map_type m_map;

  UnorderedMapErase(map_type const& map) : m_map(map) {}

  void apply() const {
    parallel_for("Kokkos::Impl::UnorderedMapErase::apply",
                 m_map.m_hash_lists.extent(0), *this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(size_type i) const {
    const size_type invalid_index = map_type::invalid_index;

    size_type curr = m_map.m_hash_lists(i);
    size_type next = invalid_index;

    // remove erased head of the linked-list
    while (curr != invalid_index && !m_map.valid_at(curr)) {
      next                     = m_map.m_next_index[curr];
      m_map.m_next_index[curr] = invalid_index;
      m_map.m_keys[curr]       = key_type();
      if (m_map.is_set) m_map.m_values[curr] = value_type();
      curr                  = next;
      m_map.m_hash_lists(i) = next;
    }

    // if the list is non-empty and the head is valid
    if (curr != invalid_index && m_map.valid_at(curr)) {
      size_type prev = curr;
      curr           = m_map.m_next_index[prev];

      while (curr != invalid_index) {
        next = m_map.m_next_index[curr];
        if (m_map.valid_at(curr)) {
          prev = curr;
        } else {
          // remove curr from list
          m_map.m_next_index[prev] = next;
          m_map.m_next_index[curr] = invalid_index;
          m_map.m_keys[curr]       = key_type();
          if (map_type::is_set) m_map.m_values[curr] = value_type();
        }
        curr = next;
      }
    }
  }
};

template <typename UMap>
struct UnorderedMapHistogram {
  using map_type        = UMap;
  using execution_space = typename map_type::execution_space;
  using size_type       = typename map_type::size_type;

  using histogram_view      = View<int[100], typename map_type::device_type>;
  using host_histogram_view = typename histogram_view::HostMirror;

  map_type m_map;
  histogram_view m_length;
  histogram_view m_distance;
  histogram_view m_block_distance;

  UnorderedMapHistogram(map_type const& map)
      : m_map(map),
        m_length("UnorderedMap Histogram"),
        m_distance("UnorderedMap Histogram"),
        m_block_distance("UnorderedMap Histogram") {}

  void calculate() {
    parallel_for("Kokkos::Impl::UnorderedMapHistogram::calculate",
                 m_map.m_hash_lists.extent(0), *this);
  }

  void clear() {
    Kokkos::deep_copy(m_length, 0);
    Kokkos::deep_copy(m_distance, 0);
    Kokkos::deep_copy(m_block_distance, 0);
  }

  void print_length(std::ostream& out) {
    host_histogram_view host_copy =
        create_mirror_view_and_copy(Kokkos::HostSpace{}, m_length);

    for (int i = 0, size = host_copy.extent(0); i < size; ++i) {
      out << host_copy[i] << " , ";
    }
    out << "\b\b\b   " << std::endl;
  }

  void print_distance(std::ostream& out) {
    host_histogram_view host_copy =
        create_mirror_view_and_copy(Kokkos::HostSpace{}, m_distance);

    for (int i = 0, size = host_copy.extent(0); i < size; ++i) {
      out << host_copy[i] << " , ";
    }
    out << "\b\b\b   " << std::endl;
  }

  void print_block_distance(std::ostream& out) {
    host_histogram_view host_copy =
        create_mirror_view_and_copy(Kokkos::HostSpace{}, m_block_distance);

    for (int i = 0, size = host_copy.extent(0); i < size; ++i) {
      out << host_copy[i] << " , ";
    }
    out << "\b\b\b   " << std::endl;
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(size_type i) const {
    const size_type invalid_index = map_type::invalid_index;

    uint32_t length     = 0;
    size_type min_index = ~0u, max_index = 0;
    for (size_type curr = m_map.m_hash_lists(i); curr != invalid_index;
         curr           = m_map.m_next_index[curr]) {
      ++length;
      min_index = (curr < min_index) ? curr : min_index;
      max_index = (max_index < curr) ? curr : max_index;
    }

    size_type distance = (0u < length) ? max_index - min_index : 0u;
    size_type blocks   = (0u < length) ? max_index / 32u - min_index / 32u : 0u;

    // normalize data
    length   = length < 100u ? length : 99u;
    distance = distance < 100u ? distance : 99u;
    blocks   = blocks < 100u ? blocks : 99u;

    if (0u < length) {
      atomic_fetch_add(&m_length(length), 1);
      atomic_fetch_add(&m_distance(distance), 1);
      atomic_fetch_add(&m_block_distance(blocks), 1);
    }
  }
};

template <typename UMap>
struct UnorderedMapPrint {
  using map_type        = UMap;
  using execution_space = typename map_type::execution_space;
  using size_type       = typename map_type::size_type;

  map_type m_map;

  UnorderedMapPrint(map_type const& map) : m_map(map) {}

  void apply() {
    parallel_for("Kokkos::Impl::UnorderedMapPrint::apply",
                 m_map.m_hash_lists.extent(0), *this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(size_type i) const {
    const size_type invalid_index = map_type::invalid_index;

    uint32_t list = m_map.m_hash_lists(i);
    for (size_type curr = list, ii = 0; curr != invalid_index;
         curr = m_map.m_next_index[curr], ++ii) {
      KOKKOS_IMPL_DO_NOT_USE_PRINTF("%d[%d]: %d->%d\n", list, ii,
                                    m_map.key_at(curr), m_map.value_at(curr));
    }
  }
};

template <typename DKey, typename DValue, typename SKey, typename SValue>
struct UnorderedMapCanAssign : public std::false_type {};

template <typename Key, typename Value>
struct UnorderedMapCanAssign<Key, Value, Key, Value> : public std::true_type {};

template <typename Key, typename Value>
struct UnorderedMapCanAssign<const Key, Value, Key, Value>
    : public std::true_type {};

template <typename Key, typename Value>
struct UnorderedMapCanAssign<const Key, const Value, Key, Value>
    : public std::true_type {};

template <typename Key, typename Value>
struct UnorderedMapCanAssign<const Key, const Value, const Key, Value>
    : public std::true_type {};

}  // namespace Impl
}  // namespace Kokkos

#endif  // KOKKOS_UNORDERED_MAP_IMPL_HPP
