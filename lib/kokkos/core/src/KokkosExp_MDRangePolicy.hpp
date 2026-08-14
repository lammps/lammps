// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#include <Kokkos_Macros.hpp>
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
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
#include <array>
#include <cmath>

namespace Kokkos {

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
constexpr std::enable_if_t<!std::is_signed_v<T>, bool>
is_less_than_value_initialized_variable(T) {
  return false;
}

template <class T>
constexpr std::enable_if_t<std::is_signed_v<T>, bool>
is_less_than_value_initialized_variable(T arg) {
  return arg < T{};
}

// Checked narrowing conversion that calls abort if the cast changes the value
template <class To, class From>
constexpr To checked_narrow_cast(From arg, std::size_t idx) {
  constexpr const bool is_different_signedness =
      (std::is_signed_v<To> != std::is_signed_v<From>);
  auto const ret = static_cast<To>(arg);  // NOLINT(bugprone-signed-char-misuse)
  if (static_cast<From>(ret) != arg ||
      (is_different_signedness &&
       is_less_than_value_initialized_variable(arg) !=
           is_less_than_value_initialized_variable(ret))) {
    auto msg =
        "Kokkos::MDRangePolicy bound type error: an unsafe implicit conversion "
        "is performed on a bound (" +
        std::to_string(arg) + ") in dimension (" + std::to_string(idx) +
        "), which may not preserve its original value.\n";
    Kokkos::abort(msg.c_str());
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
  static_assert(M <= N);
  auto* ptr = a.data();
  // NOTE equivalent to
  // std::transform(std::begin(init), std::end(init), a.data(),
  //                [](U x) { return static_cast<T>(x); });
  // except that std::transform is not constexpr.
  for (std::size_t i = 0; i < M; ++i) {
    *ptr++ = checked_narrow_cast<T>(init[i], i);
    (void)checked_narrow_cast<IndexType>(init[i], i);  // see note above
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
  static_assert(M <= N);
  for (std::size_t i = 0; i < M; ++i) {
    a[i] = checked_narrow_cast<T>(other[i], i);
    (void)checked_narrow_cast<IndexType>(other[i], i);  // see note above
  }
  return a;
}

struct TileSizeProperties {
  int max_threads;  // (per SM, CU)
  int max_total_tile_size;
  // For GPU backends: hardware limits for block dimensions
  std::array<int, 3> max_threads_dimensions;
};

template <typename ExecutionSpace>
TileSizeProperties get_tile_size_properties(const ExecutionSpace&) {
  // Host settings
  TileSizeProperties properties;
  properties.max_threads               = std::numeric_limits<int>::max();
  properties.max_total_tile_size       = std::numeric_limits<int>::max();
  properties.max_threads_dimensions[0] = std::numeric_limits<int>::max();
  properties.max_threads_dimensions[1] = std::numeric_limits<int>::max();
  properties.max_threads_dimensions[2] = std::numeric_limits<int>::max();
  return properties;
}

// Default tile size recommended (for MDRangePolicy)
template <typename ExecutionSpace>
struct TileSizeRecommended {
  template <typename Policy>
  static auto get(Policy const& policy);
};

// Recommend tile sizes for each rank of MDRangePolicy.
// Each rank is tiled with a default size of 2, except the innermost rank which
// is set to its full work range length.
template <typename ExecutionSpace>
template <typename Policy>
auto TileSizeRecommended<ExecutionSpace>::get(Policy const& policy) {
  constexpr auto InnerDirection = Policy::inner_direction;
  constexpr int Rank            = Policy::rank;

  using tile_type = Kokkos::Array<std::int64_t, Rank>;

  tile_type recommended_tile_sizes{};
  int default_tile_size   = 2;
  int max_total_tile_size = policy.max_total_tile_size();

  int inner_rank  = (InnerDirection == Iterate::Right) ? Rank - 1 : 0;
  int outer_bound = (InnerDirection == Iterate::Right) ? -1 : Rank;
  int iter_step   = (InnerDirection == Iterate::Right) ? -1 : 1;
  auto inner_work_range =
      policy.m_upper[inner_rank] - policy.m_lower[inner_rank];

  int prod_tile_size = 1;
  for (int i = inner_rank; i != outer_bound; i += iter_step) {
    int rank_tile_size = 1;
    if (prod_tile_size * default_tile_size <= max_total_tile_size) {
      rank_tile_size = default_tile_size;
    } else {
      rank_tile_size = 1;
    }
    if (i == inner_rank) {
      rank_tile_size = std::max<int>(inner_work_range, 1);
    }
    prod_tile_size *= rank_tile_size;
    recommended_tile_sizes[i] = rank_tile_size;
  }
  return recommended_tile_sizes;
}

}  // namespace Impl

// multi-dimensional iteration pattern
template <typename... Properties>
struct MDRangePolicy;

// Note: If MDRangePolicy has a primary template, implicit CTAD (deduction
// guides) are generated -> MDRangePolicy<> by some compilers, which is
// incorrect.  By making it a template specialization instead, no implicit CTAD
// is generated.  This works because there has to be at least one property
// specified (which is Rank<...>); otherwise, we'd get the static_assert
// "Kokkos::Error: MD iteration pattern not defined".  This template
// specialization uses <P, Properties...> in all places for correctness.
template <typename P, typename... Properties>
struct MDRangePolicy<P, Properties...>
    : public Kokkos::Impl::PolicyTraits<P, Properties...> {
  using traits          = Kokkos::Impl::PolicyTraits<P, Properties...>;
  using execution_space = typename traits::execution_space;
  using range_policy    = RangePolicy<P, Properties...>;

  using impl_range_policy =
      RangePolicy<execution_space, typename traits::schedule_type,
                  typename traits::index_type>;

  using execution_policy =
      MDRangePolicy<P, Properties...>;  // needed for is_execution_policy
                                        // interrogation

  template <class... OtherProperties>
  friend struct MDRangePolicy;

  static_assert(!std::is_void_v<typename traits::iteration_pattern>,
                "Kokkos Error: MD iteration pattern not defined");

  using iteration_pattern = typename traits::iteration_pattern;
  using work_tag          = typename traits::work_tag;
  using launch_bounds     = typename traits::launch_bounds;
  using member_type       = typename range_policy::member_type;

  static constexpr int rank = iteration_pattern::rank;
  static_assert(rank < 7, "Kokkos MDRangePolicy Error: Unsupported rank...");

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

  execution_space m_space;

  point_type m_lower                          = {};
  point_type m_upper                          = {};
  tile_type m_tile                            = {};
  point_type m_tile_end                       = {};
  index_type m_num_tiles                      = 1;
  index_type m_prod_tile_dims                 = 1;
  bool m_tune_tile_size                       = false;
  index_type m_max_total_tile_size            = 1;
  std::array<int, 3> m_max_threads_dimensions = {1, 1, 1};

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
            typename = std::enable_if_t<std::is_integral_v<LT> &&
                                        std::is_integral_v<UT> &&
                                        std::is_integral_v<TT>>>
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
            typename = std::enable_if_t<std::is_integral_v<LT> &&
                                        std::is_integral_v<UT> &&
                                        std::is_integral_v<TT>>>
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
    update_tiling_properties();
  }

  template <typename T, std::size_t NT = rank,
            typename = std::enable_if_t<std::is_integral_v<T>>>
  MDRangePolicy(Kokkos::Array<T, rank> const& lower,
                Kokkos::Array<T, rank> const& upper,
                Kokkos::Array<T, NT> const& tile = Kokkos::Array<T, NT>{})
      : MDRangePolicy(typename traits::execution_space(), lower, upper, tile) {}

  template <typename T, std::size_t NT = rank,
            typename = std::enable_if_t<std::is_integral_v<T>>>
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

  MDRangePolicy(const Impl::PolicyUpdate, const MDRangePolicy& other,
                typename traits::execution_space space)
      : MDRangePolicy(other) {
    this->m_space = std::move(space);
    // Reset auto-tuned tiles if the execution space changes since the computed
    // tile size may be different
    if (this->m_tune_tile_size) {
      this->m_tile = {};
    }
    update_tiling_properties();
  }

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
        m_tune_tile_size(p.m_tune_tile_size),
        m_max_total_tile_size(p.m_max_total_tile_size),
        m_max_threads_dimensions(p.m_max_threads_dimensions) {}

  void impl_change_tile_size(const point_type& tile) {
    this->m_tile = tile;
    this->update_tiling_properties();
  }

  bool impl_tune_tile_size() const { return m_tune_tile_size; }

  tile_type tile_size_recommended() const {
    return Kokkos::Impl::TileSizeRecommended<execution_space>::get(*this);
  }

  index_type max_total_tile_size() const { return m_max_total_tile_size; }

 private:
  void update_tiling_properties() {
    auto properties        = Impl::get_tile_size_properties(m_space);
    this->m_num_tiles      = 1;
    this->m_prod_tile_dims = 1;
    this->m_max_total_tile_size =
        static_cast<index_type>(properties.max_total_tile_size);
    this->m_max_threads_dimensions = properties.max_threads_dimensions;

    index_type effective_max_tile_size = this->m_max_total_tile_size;

    constexpr bool enforce_launch_bounds =
#if defined(KOKKOS_ENABLE_CUDA)
        std::is_same_v<execution_space, Kokkos::Cuda>;
#elif defined(KOKKOS_ENABLE_HIP)
        std::is_same_v<execution_space, Kokkos::HIP>;
#else
        false;
#endif

    if constexpr (enforce_launch_bounds && launch_bounds::maxTperB != 0) {
      effective_max_tile_size =
          std::min(effective_max_tile_size,
                   static_cast<index_type>(launch_bounds::maxTperB));
    }

    tile_type default_tile = this->tile_size_recommended();

    int inner_rank  = (inner_direction == Iterate::Right) ? rank - 1 : 0;
    int outer_bound = (inner_direction == Iterate::Right) ? -1 : rank;
    int iter_step   = (inner_direction == Iterate::Right) ? -1 : 1;

    for (int i = inner_rank; i != outer_bound; i += iter_step) {
      const index_type length = this->m_upper[i] - this->m_lower[i];

      if (this->m_upper[i] < this->m_lower[i]) {
        std::string msg =
            "Kokkos::MDRangePolicy bounds error: The lower bound (" +
            std::to_string(this->m_lower[i]) +
            ") is greater than its upper bound (" +
            std::to_string(this->m_upper[i]) + ") in dimension " +
            std::to_string(i) + ".\n";
        Kokkos::abort(msg.c_str());
      }

      // If tile size is not specified or <= 0 set to recommended tile size
      if (this->m_tile[i] <= 0) {
        this->m_tune_tile_size = true;
        // Set to recommended tile size if it fits within effective limit
        if (this->m_prod_tile_dims * default_tile[i] <=
            effective_max_tile_size) {
          this->m_tile[i] = default_tile[i];
        } else {
          // Try to fit within effective limit by reducing tile size
          while (default_tile[i] > 1 &&
                 this->m_prod_tile_dims * default_tile[i] >
                     effective_max_tile_size) {
            default_tile[i] >>= 1;
          }
          this->m_tile[i] = (default_tile[i] > 1) ? default_tile[i] : 1;
        }
      }

      this->m_tile_end[i] = static_cast<index_type>(
          (length + this->m_tile[i] - 1) / this->m_tile[i]);
      this->m_num_tiles *= this->m_tile_end[i];
      this->m_prod_tile_dims *= this->m_tile[i];
    }

    if constexpr (enforce_launch_bounds && launch_bounds::maxTperB != 0) {
      if (static_cast<index_type>(launch_bounds::maxTperB) <
          this->m_prod_tile_dims) {
        std::string msg =
            "Kokkos::MDRangePolicy tile dimensions error: Product of tile "
            "dimensions (" +
            std::to_string(static_cast<int>(this->m_prod_tile_dims)) +
            ") is greater than the maximum specified via LaunchBounds (" +
            std::to_string(launch_bounds::maxTperB) +
            ") - choose smaller tile dims\n";
        Kokkos::abort(msg.c_str());
      }
    }

    if (this->m_prod_tile_dims >
        static_cast<index_type>(this->m_max_total_tile_size)) {
      std::string msg =
          "Kokkos::MDRangePolicy tile dimensions error: Product of tile "
          "dimensions (" +
          std::to_string(static_cast<int>(this->m_prod_tile_dims)) +
          ") is greater than the maximum total tile size (" +
          std::to_string(static_cast<int>(this->m_max_total_tile_size)) +
          ") - choose smaller tile dims\n";
      Kokkos::abort(msg.c_str());
    }
  }
};

template <typename LT, size_t N, typename UT>
MDRangePolicy(const LT (&)[N], const UT (&)[N]) -> MDRangePolicy<Rank<N>>;

template <typename LT, size_t N, typename UT, typename TT, size_t TN>
MDRangePolicy(const LT (&)[N], const UT (&)[N], const TT (&)[TN])
    -> MDRangePolicy<Rank<N>>;

template <typename LT, size_t N, typename UT>
MDRangePolicy(DefaultExecutionSpace const&, const LT (&)[N], const UT (&)[N])
    -> MDRangePolicy<Rank<N>>;

template <typename LT, size_t N, typename UT, typename TT, size_t TN>
MDRangePolicy(DefaultExecutionSpace const&, const LT (&)[N], const UT (&)[N],
              const TT (&)[TN]) -> MDRangePolicy<Rank<N>>;

template <typename ES, typename LT, size_t N, typename UT,
          typename = std::enable_if_t<is_execution_space_v<ES>>>
MDRangePolicy(ES const&, const LT (&)[N], const UT (&)[N])
    -> MDRangePolicy<ES, Rank<N>>;

template <typename ES, typename LT, size_t N, typename UT, typename TT,
          size_t TN, typename = std::enable_if_t<is_execution_space_v<ES>>>
MDRangePolicy(ES const&, const LT (&)[N], const UT (&)[N], const TT (&)[TN])
    -> MDRangePolicy<ES, Rank<N>>;

template <typename T, size_t N>
MDRangePolicy(Array<T, N> const&, Array<T, N> const&) -> MDRangePolicy<Rank<N>>;

template <typename T, size_t N, size_t NT>
MDRangePolicy(Array<T, N> const&, Array<T, N> const&, Array<T, NT> const&)
    -> MDRangePolicy<Rank<N>>;

template <typename T, size_t N>
MDRangePolicy(DefaultExecutionSpace const&, Array<T, N> const&,
              Array<T, N> const&) -> MDRangePolicy<Rank<N>>;

template <typename T, size_t N, size_t NT>
MDRangePolicy(DefaultExecutionSpace const&, Array<T, N> const&,
              Array<T, N> const&, Array<T, NT> const&)
    -> MDRangePolicy<Rank<N>>;

template <typename ES, typename T, size_t N,
          typename = std::enable_if_t<is_execution_space_v<ES>>>
MDRangePolicy(ES const&, Array<T, N> const&, Array<T, N> const&)
    -> MDRangePolicy<ES, Rank<N>>;

template <typename ES, typename T, size_t N, size_t NT,
          typename = std::enable_if_t<is_execution_space_v<ES>>>
MDRangePolicy(ES const&, Array<T, N> const&, Array<T, N> const&,
              Array<T, NT> const&) -> MDRangePolicy<ES, Rank<N>>;

}  // namespace Kokkos

#endif  // KOKKOS_CORE_EXP_MD_RANGE_POLICY_HPP
