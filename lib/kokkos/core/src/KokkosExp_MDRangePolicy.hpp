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
#ifndef KOKKOS_CORE_EXP_MD_RANGE_POLICY_HPP
#define KOKKOS_CORE_EXP_MD_RANGE_POLICY_HPP

#include <initializer_list>

#include <Kokkos_Layout.hpp>
#include <Kokkos_Rank.hpp>
#include <Kokkos_Array.hpp>
#include <impl/KokkosExp_Host_IterateTile.hpp>
#include <Kokkos_ExecPolicy.hpp>
#include <type_traits>

namespace Kokkos {

// ------------------------------------------------------------------ //
// Moved to Kokkos_Layout.hpp for more general accessibility
/*
enum class Iterate
{
  Default, // Default for the device
  Left,    // Left indices stride fastest
  Right,   // Right indices stride fastest
};
*/

template <typename ExecSpace>
struct default_outer_direction {
  using type                     = Iterate;
  static constexpr Iterate value = Iterate::Right;
};

template <typename ExecSpace>
struct default_inner_direction {
  using type                     = Iterate;
  static constexpr Iterate value = Iterate::Right;
};

namespace Impl {
// NOTE the comparison below is encapsulated to silent warnings about pointless
// comparison of unsigned integer with zero
template <class T>
constexpr std::enable_if_t<!std::is_signed<T>::value, bool>
is_less_than_value_initialized_variable(T) {
  return false;
}

template <class T>
constexpr std::enable_if_t<std::is_signed<T>::value, bool>
is_less_than_value_initialized_variable(T arg) {
  return arg < T{};
}

// Checked narrowing conversion that calls abort if the cast changes the value
template <class To, class From>
constexpr To checked_narrow_cast(From arg) {
  constexpr const bool is_different_signedness =
      (std::is_signed<To>::value != std::is_signed<From>::value);
  auto const ret = static_cast<To>(arg);
  if (static_cast<From>(ret) != arg ||
      (is_different_signedness &&
       is_less_than_value_initialized_variable(arg) !=
           is_less_than_value_initialized_variable(ret))) {
    Kokkos::abort("unsafe narrowing conversion");
  }
  return ret;
}
// NOTE prefer C array U[M] to std::initalizer_list<U> so that the number of
// elements can be deduced (https://stackoverflow.com/q/40241370)
// NOTE for some unfortunate reason the policy bounds are stored as signed
// integer arrays (point_type which is Kokkos::Array<std::int64_t>) so we
// specify the index type (actual policy index_type from the traits) and check
// ahead of time that narrowing conversions will be safe.
template <class IndexType, class Array, class U, std::size_t M>
constexpr Array to_array_potentially_narrowing(const U (&init)[M]) {
  using T = typename Array::value_type;
  Array a{};
  constexpr std::size_t N = a.size();
  static_assert(M <= N, "");
  auto* ptr = a.data();
  // NOTE equivalent to
  // std::transform(std::begin(init), std::end(init), a.data(),
  //                [](U x) { return static_cast<T>(x); });
  // except that std::transform is not constexpr.
  for (auto x : init) {
    *ptr++ = checked_narrow_cast<T>(x);
    (void)checked_narrow_cast<IndexType>(x);  // see note above
  }
  return a;
}

// NOTE Making a copy even when std::is_same<Array, Kokkos::Array<U, M>>::value
// is true to reduce code complexity.  You may change this if you have a good
// reason to.  Intentionally not enabling std::array at this time but this may
// change too.
template <class IndexType, class NVCC_WONT_LET_ME_CALL_YOU_Array, class U,
          std::size_t M>
constexpr NVCC_WONT_LET_ME_CALL_YOU_Array to_array_potentially_narrowing(
    Kokkos::Array<U, M> const& other) {
  using T = typename NVCC_WONT_LET_ME_CALL_YOU_Array::value_type;
  NVCC_WONT_LET_ME_CALL_YOU_Array a{};
  constexpr std::size_t N = a.size();
  static_assert(M <= N, "");
  for (std::size_t i = 0; i < M; ++i) {
    a[i] = checked_narrow_cast<T>(other[i]);
    (void)checked_narrow_cast<IndexType>(other[i]);  // see note above
  }
  return a;
}

struct TileSizeProperties {
  int max_threads;
  int default_largest_tile_size;
  int default_tile_size;
  int max_total_tile_size;
};

template <typename ExecutionSpace>
TileSizeProperties get_tile_size_properties(const ExecutionSpace&) {
  // Host settings
  TileSizeProperties properties;
  properties.max_threads               = std::numeric_limits<int>::max();
  properties.default_largest_tile_size = 0;
  properties.default_tile_size         = 2;
  properties.max_total_tile_size       = std::numeric_limits<int>::max();
  return properties;
}

}  // namespace Impl

// multi-dimensional iteration pattern
template <typename... Properties>
struct MDRangePolicy : public Kokkos::Impl::PolicyTraits<Properties...> {
  using traits       = Kokkos::Impl::PolicyTraits<Properties...>;
  using range_policy = RangePolicy<Properties...>;

  typename traits::execution_space m_space;

  using impl_range_policy =
      RangePolicy<typename traits::execution_space,
                  typename traits::schedule_type, typename traits::index_type>;

  using execution_policy =
      MDRangePolicy<Properties...>;  // needed for is_execution_space
                                     // interrogation

  template <class... OtherProperties>
  friend struct MDRangePolicy;

  static_assert(!std::is_void<typename traits::iteration_pattern>::value,
                "Kokkos Error: MD iteration pattern not defined");

  using iteration_pattern = typename traits::iteration_pattern;
  using work_tag          = typename traits::work_tag;
  using launch_bounds     = typename traits::launch_bounds;
  using member_type       = typename range_policy::member_type;

  static constexpr int rank = iteration_pattern::rank;

  using index_type       = typename traits::index_type;
  using array_index_type = std::int64_t;
  using point_type = Kokkos::Array<array_index_type, rank>;  // was index_type
  using tile_type  = Kokkos::Array<array_index_type, rank>;
  // If point_type or tile_type is not templated on a signed integral type (if
  // it is unsigned), then if user passes in intializer_list of
  // runtime-determined values of signed integral type that are not const will
  // receive a compiler error due to an invalid case for implicit conversion -
  // "conversion from integer or unscoped enumeration type to integer type that
  // cannot represent all values of the original, except where source is a
  // constant expression whose value can be stored exactly in the target type"
  // This would require the user to either pass a matching index_type parameter
  // as template parameter to the MDRangePolicy or static_cast the individual
  // values

  point_type m_lower          = {};
  point_type m_upper          = {};
  tile_type m_tile            = {};
  point_type m_tile_end       = {};
  index_type m_num_tiles      = 1;
  index_type m_prod_tile_dims = 1;
  bool m_tune_tile_size       = false;

  static constexpr auto outer_direction =
      (iteration_pattern::outer_direction != Iterate::Default)
          ? iteration_pattern::outer_direction
          : default_outer_direction<typename traits::execution_space>::value;

  static constexpr auto inner_direction =
      iteration_pattern::inner_direction != Iterate::Default
          ? iteration_pattern::inner_direction
          : default_inner_direction<typename traits::execution_space>::value;

  static constexpr auto Right = Iterate::Right;
  static constexpr auto Left  = Iterate::Left;

  KOKKOS_INLINE_FUNCTION const typename traits::execution_space& space() const {
    return m_space;
  }

  MDRangePolicy() = default;

  template <typename LT, std::size_t LN, typename UT, std::size_t UN,
            typename TT = array_index_type, std::size_t TN = rank,
            typename = std::enable_if_t<std::is_integral<LT>::value &&
                                        std::is_integral<UT>::value &&
                                        std::is_integral<TT>::value>>
  MDRangePolicy(const LT (&lower)[LN], const UT (&upper)[UN],
                const TT (&tile)[TN] = {})
      : MDRangePolicy(
            Impl::to_array_potentially_narrowing<index_type, decltype(m_lower)>(
                lower),
            Impl::to_array_potentially_narrowing<index_type, decltype(m_upper)>(
                upper),
            Impl::to_array_potentially_narrowing<index_type, decltype(m_tile)>(
                tile)) {
    static_assert(
        LN == rank && UN == rank && TN <= rank,
        "MDRangePolicy: Constructor initializer lists have wrong size");
  }

  template <typename LT, std::size_t LN, typename UT, std::size_t UN,
            typename TT = array_index_type, std::size_t TN = rank,
            typename = std::enable_if_t<std::is_integral<LT>::value &&
                                        std::is_integral<UT>::value &&
                                        std::is_integral<TT>::value>>
  MDRangePolicy(const typename traits::execution_space& work_space,
                const LT (&lower)[LN], const UT (&upper)[UN],
                const TT (&tile)[TN] = {})
      : MDRangePolicy(
            work_space,
            Impl::to_array_potentially_narrowing<index_type, decltype(m_lower)>(
                lower),
            Impl::to_array_potentially_narrowing<index_type, decltype(m_upper)>(
                upper),
            Impl::to_array_potentially_narrowing<index_type, decltype(m_tile)>(
                tile)) {
    static_assert(
        LN == rank && UN == rank && TN <= rank,
        "MDRangePolicy: Constructor initializer lists have wrong size");
  }

  // NOTE: Keeping these two constructor despite the templated constructors
  // from Kokkos arrays for backwards compability to allow construction from
  // double-braced initializer lists.
  MDRangePolicy(point_type const& lower, point_type const& upper,
                tile_type const& tile = tile_type{})
      : MDRangePolicy(typename traits::execution_space(), lower, upper, tile) {}

  MDRangePolicy(const typename traits::execution_space& work_space,
                point_type const& lower, point_type const& upper,
                tile_type const& tile = tile_type{})
      : m_space(work_space), m_lower(lower), m_upper(upper), m_tile(tile) {
    init_helper(Impl::get_tile_size_properties(work_space));
  }

  template <typename T, std::size_t NT = rank,
            typename = std::enable_if_t<std::is_integral<T>::value>>
  MDRangePolicy(Kokkos::Array<T, rank> const& lower,
                Kokkos::Array<T, rank> const& upper,
                Kokkos::Array<T, NT> const& tile = Kokkos::Array<T, NT>{})
      : MDRangePolicy(typename traits::execution_space(), lower, upper, tile) {}

  template <typename T, std::size_t NT = rank,
            typename = std::enable_if_t<std::is_integral<T>::value>>
  MDRangePolicy(const typename traits::execution_space& work_space,
                Kokkos::Array<T, rank> const& lower,
                Kokkos::Array<T, rank> const& upper,
                Kokkos::Array<T, NT> const& tile = Kokkos::Array<T, NT>{})
      : MDRangePolicy(
            work_space,
            Impl::to_array_potentially_narrowing<index_type, decltype(m_lower)>(
                lower),
            Impl::to_array_potentially_narrowing<index_type, decltype(m_upper)>(
                upper),
            Impl::to_array_potentially_narrowing<index_type, decltype(m_tile)>(
                tile)) {}

  template <class... OtherProperties>
  MDRangePolicy(const MDRangePolicy<OtherProperties...> p)
      : traits(p),  // base class may contain data such as desired occupancy
        m_space(p.m_space),
        m_lower(p.m_lower),
        m_upper(p.m_upper),
        m_tile(p.m_tile),
        m_tile_end(p.m_tile_end),
        m_num_tiles(p.m_num_tiles),
        m_prod_tile_dims(p.m_prod_tile_dims),
        m_tune_tile_size(p.m_tune_tile_size) {}

  void impl_change_tile_size(const point_type& tile) {
    m_tile = tile;
    init_helper(Impl::get_tile_size_properties(m_space));
  }
  bool impl_tune_tile_size() const { return m_tune_tile_size; }

 private:
  void init_helper(Impl::TileSizeProperties properties) {
    m_prod_tile_dims = 1;
    int increment    = 1;
    int rank_start   = 0;
    int rank_end     = rank;
    if (inner_direction == Iterate::Right) {
      increment  = -1;
      rank_start = rank - 1;
      rank_end   = -1;
    }
    for (int i = rank_start; i != rank_end; i += increment) {
      const index_type length = m_upper[i] - m_lower[i];
      if (m_tile[i] <= 0) {
        m_tune_tile_size = true;
        if ((inner_direction == Iterate::Right && (i < rank - 1)) ||
            (inner_direction == Iterate::Left && (i > 0))) {
          if (m_prod_tile_dims * properties.default_tile_size <
              static_cast<index_type>(properties.max_total_tile_size)) {
            m_tile[i] = properties.default_tile_size;
          } else {
            m_tile[i] = 1;
          }
        } else {
          m_tile[i] = properties.default_largest_tile_size == 0
                          ? std::max<int>(length, 1)
                          : properties.default_largest_tile_size;
        }
      }
      m_tile_end[i] =
          static_cast<index_type>((length + m_tile[i] - 1) / m_tile[i]);
      m_num_tiles *= m_tile_end[i];
      m_prod_tile_dims *= m_tile[i];
    }
    if (m_prod_tile_dims > static_cast<index_type>(properties.max_threads)) {
      printf(" Product of tile dimensions exceed maximum limit: %d\n",
             static_cast<int>(properties.max_threads));
      Kokkos::abort(
          "ExecSpace Error: MDRange tile dims exceed maximum number "
          "of threads per block - choose smaller tile dims");
    }
  }
};

}  // namespace Kokkos

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_3
// For backward compatibility
namespace Kokkos {
namespace Experimental {
using Iterate KOKKOS_DEPRECATED = Kokkos::Iterate;
template <typename... Properties>
using MDRangePolicy KOKKOS_DEPRECATED = Kokkos::MDRangePolicy<Properties...>;
template <unsigned N, Kokkos::Iterate OuterDir = Kokkos::Iterate::Default,
          Kokkos::Iterate InnerDir = Kokkos::Iterate::Default>
using Rank KOKKOS_DEPRECATED = Kokkos::Rank<N, OuterDir, InnerDir>;
}  // namespace Experimental
}  // namespace Kokkos
#endif

#endif  // KOKKOS_CORE_EXP_MD_RANGE_POLICY_HPP
