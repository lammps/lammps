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
#pragma once

#include "macros.hpp"
#include "extents.hpp"
#include "trait_backports.hpp"
#include "compressed_pair.hpp"
#include "utility.hpp"

#if !defined(MDSPAN_IMPL_USE_ATTRIBUTE_NO_UNIQUE_ADDRESS)
#  include "no_unique_address.hpp"
#endif

#include <array>
#include <type_traits>
#include <utility>

#ifdef __cpp_lib_span
#include <span>
#endif
#if defined(MDSPAN_IMPL_USE_CONCEPTS) && MDSPAN_HAS_CXX_20 && defined(__cpp_lib_concepts)
#  include <concepts>
#endif

namespace MDSPAN_IMPL_STANDARD_NAMESPACE {

struct layout_left {
  template<class Extents>
  class mapping;
};
struct layout_right {
  template<class Extents>
  class mapping;
};

namespace detail {
  template<class Layout, class Mapping>
  constexpr bool is_mapping_of =
    std::is_same<typename Layout::template mapping<typename Mapping::extents_type>, Mapping>::value;

#if defined(MDSPAN_IMPL_USE_CONCEPTS) && MDSPAN_HAS_CXX_20
#  if !defined(__cpp_lib_concepts)
  namespace internal {
  namespace detail {
  template <typename Tp, typename _Up>
  concept same_as = std::is_same_v<Tp, _Up>;
  } // namespace detail
  template <class T, class U>
  concept same_as = detail::same_as<T, U> && detail::same_as<U, T>;
  } // namespace internal
#  endif

  template<class M>
  concept layout_mapping_alike = requires {
    requires impl_is_extents<typename M::extents_type>::value;
#if defined(__cpp_lib_concepts)
    { M::is_always_strided() } -> std::same_as<bool>;
    { M::is_always_exhaustive() } -> std::same_as<bool>;
    { M::is_always_unique() } -> std::same_as<bool>;
#else
    { M::is_always_strided() } -> internal::same_as<bool>;
    { M::is_always_exhaustive() } -> internal::_ame_as<bool>;
    { M::is_always_unique() } -> internal::same_as<bool>;
#endif
    std::bool_constant<M::is_always_strided()>::value;
    std::bool_constant<M::is_always_exhaustive()>::value;
    std::bool_constant<M::is_always_unique()>::value;
  };
#endif

} // namespace detail

struct layout_stride {
  template <class Extents>
  class mapping
#if !defined(MDSPAN_IMPL_USE_ATTRIBUTE_NO_UNIQUE_ADDRESS)
    : private detail::no_unique_address_emulation<
        detail::impl_compressed_pair<
          Extents,
          detail::possibly_empty_array<typename Extents::index_type, Extents::rank()>
        >
      >
#endif
  {
  public:
    using extents_type = Extents;
    using index_type = typename extents_type::index_type;
    using size_type = typename extents_type::size_type;
    using rank_type = typename extents_type::rank_type;
    using layout_type = layout_stride;

    // This could be a `requires`, but I think it's better and clearer as a `static_assert`.
    static_assert(detail::impl_is_extents_v<Extents>,
                  MDSPAN_IMPL_STANDARD_NAMESPACE_STRING "::layout_stride::mapping must be instantiated with a specialization of " MDSPAN_IMPL_STANDARD_NAMESPACE_STRING "::extents.");


  private:

    //----------------------------------------------------------------------------

    using strides_storage_t = detail::possibly_empty_array<index_type, extents_type::rank()>;
    using member_pair_t = detail::impl_compressed_pair<extents_type, strides_storage_t>;

#if defined(MDSPAN_IMPL_USE_ATTRIBUTE_NO_UNIQUE_ADDRESS)
    MDSPAN_IMPL_NO_UNIQUE_ADDRESS member_pair_t m_members;
#else
    using base_t = detail::no_unique_address_emulation<member_pair_t>;
#endif

    MDSPAN_FORCE_INLINE_FUNCTION constexpr strides_storage_t const&
    strides_storage() const noexcept {
#if defined(MDSPAN_IMPL_USE_ATTRIBUTE_NO_UNIQUE_ADDRESS)
      return m_members.second();
#else
      return this->base_t::ref().second();
#endif
    }
    MDSPAN_FORCE_INLINE_FUNCTION MDSPAN_IMPL_CONSTEXPR_14 strides_storage_t&
    strides_storage() noexcept {
#if defined(MDSPAN_IMPL_USE_ATTRIBUTE_NO_UNIQUE_ADDRESS)
      return m_members.second();
#else
      return this->base_t::ref().second();
#endif
    }

    template<class SizeType, size_t ... Ep, size_t ... Idx>
    MDSPAN_IMPL_HOST_DEVICE
    constexpr index_type get_size(::MDSPAN_IMPL_STANDARD_NAMESPACE::extents<SizeType, Ep...>,std::integer_sequence<size_t, Idx...>) const {
      return MDSPAN_IMPL_FOLD_TIMES_RIGHT( static_cast<index_type>(extents().extent(Idx)), 1 );
    }

    //----------------------------------------------------------------------------

    template <class>
    friend class mapping;

    //----------------------------------------------------------------------------

    // Workaround for non-deducibility of the index sequence template parameter if it's given at the top level
    template <class>
    struct deduction_workaround;

    template <size_t... Idxs>
    struct deduction_workaround<std::index_sequence<Idxs...>>
    {
      template <class OtherExtents>
      MDSPAN_INLINE_FUNCTION
      static constexpr bool _eq_impl(mapping const& self, mapping<OtherExtents> const& other) noexcept {
        using common_t = std::common_type_t<index_type, typename OtherExtents::index_type>;
        return    MDSPAN_IMPL_FOLD_AND((static_cast<common_t>(self.stride(Idxs)) == static_cast<common_t>(other.stride(Idxs))) /* && ... */)
               && MDSPAN_IMPL_FOLD_AND((static_cast<common_t>(self.extents().extent(Idxs)) == static_cast<common_t>(other.extents().extent(Idxs))) /* || ... */);
      }
      template <class OtherExtents>
      MDSPAN_INLINE_FUNCTION
      static constexpr bool _not_eq_impl(mapping const& self, mapping<OtherExtents> const& other) noexcept {
        using common_t = std::common_type_t<index_type, typename OtherExtents::index_type>;
        return    MDSPAN_IMPL_FOLD_OR((static_cast<common_t>(self.stride(Idxs)) != static_cast<common_t>(other.stride(Idxs))) /* || ... */)
               || MDSPAN_IMPL_FOLD_OR((static_cast<common_t>(self.extents().extent(Idxs)) != static_cast<common_t>(other.extents().extent(Idxs))) /* || ... */);
      }

      template <class... Integral>
      MDSPAN_FORCE_INLINE_FUNCTION
      static constexpr size_t _call_op_impl(mapping const& self, Integral... idxs) noexcept {
        return MDSPAN_IMPL_FOLD_PLUS_RIGHT((idxs * self.stride(Idxs)), /* + ... + */ 0);
      }

      MDSPAN_INLINE_FUNCTION
      static constexpr size_t _req_span_size_impl(mapping const& self) noexcept {
        // assumes no negative strides; not sure if I'm allowed to assume that or not
        return deduction_workaround_impl::_call_op_impl(self, (self.extents().template extent<Idxs>() - 1)...) + 1;
      }

      template<class OtherMapping>
      MDSPAN_INLINE_FUNCTION
      static constexpr const strides_storage_t fill_strides(const OtherMapping& map) {
        return strides_storage_t{static_cast<index_type>(map.stride(Idxs))...};
      }

      MDSPAN_INLINE_FUNCTION
      static constexpr const strides_storage_t& fill_strides(const strides_storage_t& s) {
        return s;
      }

      template<class IntegralType>
      static constexpr const strides_storage_t fill_strides(const std::array<IntegralType,extents_type::rank()>& s) {
        return strides_storage_t{static_cast<index_type>(s[Idxs])...};
      }

      MDSPAN_TEMPLATE_REQUIRES(
        class IntegralType,
        (std::is_convertible<IntegralType, typename extents_type::index_type>::value)
      )
      MDSPAN_INLINE_FUNCTION
      // Need to avoid zero length c-array
      static constexpr const strides_storage_t fill_strides(mdspan_non_standard_tag, const IntegralType (&s)[extents_type::rank()>0?extents_type::rank():1]) {
        return strides_storage_t{static_cast<index_type>(s[Idxs])...};
      }

#ifdef __cpp_lib_span
      template<class IntegralType>
      static constexpr const strides_storage_t fill_strides(const std::span<IntegralType,extents_type::rank()>& s) {
        return strides_storage_t{static_cast<index_type>(s[Idxs])...};
      }
#endif

      MDSPAN_INLINE_FUNCTION
      static constexpr std::array<index_type, extents_type::rank()> return_strides(const strides_storage_t& s) {
        return std::array<index_type, extents_type::rank()>{s[Idxs]...};
      }

      template<size_t K>
      MDSPAN_INLINE_FUNCTION
      static constexpr size_t return_zero() { return 0; }

      template<class Mapping>
      MDSPAN_INLINE_FUNCTION
      static constexpr typename Mapping::index_type
        offset(const Mapping& m) { return m(return_zero<Idxs>()...); }
    };

    // Can't use defaulted parameter in the deduction_workaround template because of a bug in MSVC warning C4348.
    using deduction_workaround_impl = deduction_workaround<std::make_index_sequence<Extents::rank()>>;

    MDSPAN_FUNCTION
    static constexpr strides_storage_t strides_storage(detail::with_rank<0>) {
      return {};
    }

    template <std::size_t N>
    MDSPAN_FUNCTION
    static constexpr strides_storage_t strides_storage(detail::with_rank<N>) {
      strides_storage_t s{};

      extents_type e;
      index_type stride = 1;
      for(int r = static_cast<int>(extents_type::rank() - 1); r >= 0; r--) {
        s[r] = stride;
        stride *= e.extent(r);
      }

      return s;
    }

    //----------------------------------------------------------------------------

#if defined(MDSPAN_IMPL_USE_ATTRIBUTE_NO_UNIQUE_ADDRESS)
    MDSPAN_INLINE_FUNCTION constexpr explicit
    mapping(member_pair_t&& m) : m_members(::std::move(m)) {}
#else
    MDSPAN_INLINE_FUNCTION constexpr explicit
    mapping(base_t&& __b) : base_t(::std::move(__b)) {}
#endif

  public:

    //--------------------------------------------------------------------------------

    MDSPAN_INLINE_FUNCTION constexpr mapping() noexcept
#if defined(MDSPAN_IMPL_USE_ATTRIBUTE_NO_UNIQUE_ADDRESS)
      : m_members{
#else
      : base_t(base_t{member_pair_t(
#endif
          extents_type(),
          strides_storage_t(strides_storage(detail::with_rank<extents_type::rank()>{}))
#if defined(MDSPAN_IMPL_USE_ATTRIBUTE_NO_UNIQUE_ADDRESS)
        }
#else
        )})
#endif
    {}

    MDSPAN_INLINE_FUNCTION_DEFAULTED constexpr mapping(mapping const&) noexcept = default;

    MDSPAN_TEMPLATE_REQUIRES(
      class IntegralTypes,
      /* requires */ (
        // MSVC 19.32 does not like using index_type here, requires the typename Extents::index_type
        // error C2641: cannot deduce template arguments for 'MDSPAN_IMPL_STANDARD_NAMESPACE::layout_stride::mapping'
        MDSPAN_IMPL_TRAIT(std::is_convertible, const std::remove_const_t<IntegralTypes>&, typename Extents::index_type) &&
        MDSPAN_IMPL_TRAIT(std::is_nothrow_constructible, typename Extents::index_type, const std::remove_const_t<IntegralTypes>&)
      )
    )
    constexpr
    mapping(
      extents_type const& e,
      std::array<IntegralTypes, extents_type::rank()> const& s
    ) noexcept
#if defined(MDSPAN_IMPL_USE_ATTRIBUTE_NO_UNIQUE_ADDRESS)
      : m_members{
#else
      : base_t(base_t{member_pair_t(
#endif
          e, strides_storage_t(deduction_workaround_impl::fill_strides(s))
#if defined(MDSPAN_IMPL_USE_ATTRIBUTE_NO_UNIQUE_ADDRESS)
        }
#else
        )})
#endif
    {
      /*
       * TODO: check preconditions
       * - s[i] > 0 is true for all i in the range [0, rank_ ).
       * - REQUIRED-SPAN-SIZE(e, s) is a representable value of type index_type ([basic.fundamental]).
       * - If rank_ is greater than 0, then there exists a permutation P of the integers in the
       *   range [0, rank_), such that s[ pi ] >= s[ pi − 1 ] * e.extent( pi − 1 ) is true for
       *   all i in the range [1, rank_ ), where pi is the ith element of P.
       */
    }

    MDSPAN_TEMPLATE_REQUIRES(
      class IntegralTypes,
      /* requires */ (
        // MSVC 19.32 does not like using index_type here, requires the typename Extents::index_type
        // error C2641: cannot deduce template arguments for 'MDSPAN_IMPL_STANDARD_NAMESPACE::layout_stride::mapping'
        MDSPAN_IMPL_TRAIT(std::is_convertible, const std::remove_const_t<IntegralTypes>&, typename Extents::index_type) &&
        MDSPAN_IMPL_TRAIT(std::is_nothrow_constructible, typename Extents::index_type, const std::remove_const_t<IntegralTypes>&)
      )
    )
    MDSPAN_INLINE_FUNCTION
    constexpr
    mapping(
      mdspan_non_standard_tag,
      extents_type const& e,
      // Need to avoid zero-length c-array
      const IntegralTypes (&s)[extents_type::rank()>0?extents_type::rank():1]
    ) noexcept
#if defined(MDSPAN_IMPL_USE_ATTRIBUTE_NO_UNIQUE_ADDRESS)
      : m_members{
#else
      : base_t(base_t{member_pair_t(
#endif
          e, strides_storage_t(deduction_workaround_impl::fill_strides(mdspan_non_standard, s))
#if defined(MDSPAN_IMPL_USE_ATTRIBUTE_NO_UNIQUE_ADDRESS)
        }
#else
        )})
#endif
    {
      /*
       * TODO: check preconditions
       * - s[i] > 0 is true for all i in the range [0, rank_ ).
       * - REQUIRED-SPAN-SIZE(e, s) is a representable value of type index_type ([basic.fundamental]).
       * - If rank_ is greater than 0, then there exists a permutation P of the integers in the
       *   range [0, rank_), such that s[ pi ] >= s[ pi − 1 ] * e.extent( pi − 1 ) is true for
       *   all i in the range [1, rank_ ), where pi is the ith element of P.
       */
    }

#ifdef __cpp_lib_span
    MDSPAN_TEMPLATE_REQUIRES(
      class IntegralTypes,
      /* requires */ (
        // MSVC 19.32 does not like using index_type here, requires the typename Extents::index_type
        // error C2641: cannot deduce template arguments for 'MDSPAN_IMPL_STANDARD_NAMESPACE::layout_stride::mapping'
        MDSPAN_IMPL_TRAIT(std::is_convertible, const std::remove_const_t<IntegralTypes>&, typename Extents::index_type) &&
        MDSPAN_IMPL_TRAIT(std::is_nothrow_constructible, typename Extents::index_type, const std::remove_const_t<IntegralTypes>&)
      )
    )
    constexpr
    mapping(
      extents_type const& e,
      std::span<IntegralTypes, extents_type::rank()> const& s
    ) noexcept
#if defined(MDSPAN_IMPL_USE_ATTRIBUTE_NO_UNIQUE_ADDRESS)
      : m_members{
#else
      : base_t(base_t{member_pair_t(
#endif
          e, strides_storage_t(deduction_workaround_impl::fill_strides(s))
#if defined(MDSPAN_IMPL_USE_ATTRIBUTE_NO_UNIQUE_ADDRESS)
        }
#else
        )})
#endif
    {
      /*
       * TODO: check preconditions
       * - s[i] > 0 is true for all i in the range [0, rank_ ).
       * - REQUIRED-SPAN-SIZE(e, s) is a representable value of type index_type ([basic.fundamental]).
       * - If rank_ is greater than 0, then there exists a permutation P of the integers in the
       *   range [0, rank_), such that s[ pi ] >= s[ pi − 1 ] * e.extent( pi − 1 ) is true for
       *   all i in the range [1, rank_ ), where pi is the ith element of P.
       */
    }
#endif // __cpp_lib_span

#if !(defined(MDSPAN_IMPL_USE_CONCEPTS) && MDSPAN_HAS_CXX_20)
    MDSPAN_TEMPLATE_REQUIRES(
      class StridedLayoutMapping,
      /* requires */ (
        MDSPAN_IMPL_TRAIT(std::is_constructible, extents_type, typename StridedLayoutMapping::extents_type) &&
        detail::is_mapping_of<typename StridedLayoutMapping::layout_type, StridedLayoutMapping> &&
        StridedLayoutMapping::is_always_unique() &&
        StridedLayoutMapping::is_always_strided()
      )
    )
#else
    template<class StridedLayoutMapping>
    requires(
         detail::layout_mapping_alike<StridedLayoutMapping> &&
         MDSPAN_IMPL_TRAIT(std::is_constructible, extents_type, typename StridedLayoutMapping::extents_type) &&
         StridedLayoutMapping::is_always_unique() &&
         StridedLayoutMapping::is_always_strided()
    )
#endif
    MDSPAN_CONDITIONAL_EXPLICIT(
      !(std::is_convertible<typename StridedLayoutMapping::extents_type, extents_type>::value &&
       (detail::is_mapping_of<layout_left, StridedLayoutMapping> ||
        detail::is_mapping_of<layout_right, StridedLayoutMapping> ||
        detail::is_mapping_of<layout_stride, StridedLayoutMapping>))
    ) // needs two () due to comma
    MDSPAN_INLINE_FUNCTION MDSPAN_IMPL_CONSTEXPR_14
    mapping(StridedLayoutMapping const& other) noexcept // NOLINT(google-explicit-constructor)
#if defined(MDSPAN_IMPL_USE_ATTRIBUTE_NO_UNIQUE_ADDRESS)
      : m_members{
#else
      : base_t(base_t{member_pair_t(
#endif
          other.extents(), strides_storage_t(deduction_workaround_impl::fill_strides(other))
#if defined(MDSPAN_IMPL_USE_ATTRIBUTE_NO_UNIQUE_ADDRESS)
        }
#else
        )})
#endif
    {
      /*
       * TODO: check preconditions
       * - other.stride(i) > 0 is true for all i in the range [0, rank_ ).
       * - other.required_span_size() is a representable value of type index_type ([basic.fundamental]).
       * - OFFSET(other) == 0
       */
    }

    //--------------------------------------------------------------------------------

    MDSPAN_INLINE_FUNCTION_DEFAULTED MDSPAN_IMPL_CONSTEXPR_14_DEFAULTED
    mapping& operator=(mapping const&) noexcept = default;

    MDSPAN_INLINE_FUNCTION constexpr const extents_type& extents() const noexcept {
#if defined(MDSPAN_IMPL_USE_ATTRIBUTE_NO_UNIQUE_ADDRESS)
      return m_members.first();
#else
      return this->base_t::ref().first();
#endif
    };

    MDSPAN_INLINE_FUNCTION
    constexpr std::array< index_type, extents_type::rank() > strides() const noexcept {
      return deduction_workaround_impl::return_strides(strides_storage());
    }

    MDSPAN_INLINE_FUNCTION
    constexpr index_type required_span_size() const noexcept {
      index_type span_size = 1;
      // using int here to avoid warning about pointless comparison to 0
      for(int r = 0; r < static_cast<int>(extents_type::rank()); r++) {
        // Return early if any of the extents are zero
        if(extents().extent(r)==0) return 0;
        span_size += ( static_cast<index_type>(extents().extent(r) - 1 ) * strides_storage()[r]);
      }
      return span_size;
    }


    MDSPAN_TEMPLATE_REQUIRES(
      class... Indices,
      /* requires */ (
        sizeof...(Indices) == Extents::rank() &&
        (detail::are_valid_indices<index_type, Indices...>())
      )
    )
    MDSPAN_FORCE_INLINE_FUNCTION
    constexpr index_type operator()(Indices... idxs) const noexcept {
#if ! defined(NDEBUG)
      detail::check_all_indices(this->extents(), idxs...);
#endif // ! NDEBUG
      return static_cast<index_type>(deduction_workaround_impl::_call_op_impl(*this, static_cast<index_type>(idxs)...));
    }

    MDSPAN_INLINE_FUNCTION static constexpr bool is_always_unique() noexcept { return true; }
    MDSPAN_INLINE_FUNCTION static constexpr bool is_always_exhaustive() noexcept {
      return false;
    }
    MDSPAN_INLINE_FUNCTION static constexpr bool is_always_strided() noexcept { return true; }

    MDSPAN_INLINE_FUNCTION static constexpr bool is_unique() noexcept { return true; }

  private:
    MDSPAN_INLINE_FUNCTION
    constexpr bool exhaustive_for_nonzero_span_size() const
    {
      return required_span_size() == get_size(extents(), std::make_index_sequence<extents_type::rank()>());
    }

    MDSPAN_INLINE_FUNCTION
    constexpr bool is_exhaustive_impl(detail::with_rank<0>) const
    {
      return true;
    }
    MDSPAN_INLINE_FUNCTION
    constexpr bool is_exhaustive_impl(detail::with_rank<1>) const
    {
      if (required_span_size() != static_cast<index_type>(0)) {
        return exhaustive_for_nonzero_span_size();
      }
      return stride(0) == 1;
    }
    template <std::size_t N>
    MDSPAN_INLINE_FUNCTION
    constexpr bool is_exhaustive_impl(detail::with_rank<N>) const
    {
      if (required_span_size() != static_cast<index_type>(0)) {
        return exhaustive_for_nonzero_span_size();
      }

      rank_type r_largest = 0;
      for (rank_type r = 1; r < extents_type::rank(); r++) {
        if (stride(r) > stride(r_largest)) {
          r_largest = r;
        }
      }
      for (rank_type r = 0; r < extents_type::rank(); r++) {
        if (extents().extent(r) == 0 && r != r_largest) {
          return false;
        }
      }
      return true;
    }

  public:
    MDSPAN_INLINE_FUNCTION MDSPAN_IMPL_CONSTEXPR_14 bool is_exhaustive() const noexcept {
      return is_exhaustive_impl(detail::with_rank<extents_type::rank()>{});
    }
    MDSPAN_INLINE_FUNCTION static constexpr bool is_strided() noexcept { return true; }


    MDSPAN_INLINE_FUNCTION
    constexpr index_type stride(rank_type r) const noexcept {
      return strides_storage()[r];
    }

#if !(defined(MDSPAN_IMPL_USE_CONCEPTS) && MDSPAN_HAS_CXX_20)
    MDSPAN_TEMPLATE_REQUIRES(
      class StridedLayoutMapping,
      /* requires */ (
        detail::is_mapping_of<typename StridedLayoutMapping::layout_type, StridedLayoutMapping> &&
        (extents_type::rank() == StridedLayoutMapping::extents_type::rank()) &&
        StridedLayoutMapping::is_always_strided()
      )
    )
#else
    template<class StridedLayoutMapping>
    requires(
         detail::layout_mapping_alike<StridedLayoutMapping> &&
         (extents_type::rank() == StridedLayoutMapping::extents_type::rank()) &&
         StridedLayoutMapping::is_always_strided()
    )
#endif
    MDSPAN_INLINE_FUNCTION
    friend constexpr bool operator==(const mapping& x, const StridedLayoutMapping& y) noexcept {
      return (x.extents() == y.extents()) &&
             (deduction_workaround_impl::offset(y) == static_cast<typename StridedLayoutMapping::index_type>(0)) &&
             detail::rankwise_equal(detail::with_rank<extents_type::rank()>{}, x, y, detail::stride);
    }

    // This one is not technically part of the proposal. Just here to make implementation a bit more optimal hopefully
    MDSPAN_TEMPLATE_REQUIRES(
      class OtherExtents,
      /* requires */ (
        (extents_type::rank() == OtherExtents::rank())
      )
    )
    MDSPAN_INLINE_FUNCTION
    friend constexpr bool operator==(mapping const& lhs, mapping<OtherExtents> const& rhs) noexcept {
      return deduction_workaround_impl::_eq_impl(lhs, rhs);
    }

#if !MDSPAN_HAS_CXX_20
    MDSPAN_TEMPLATE_REQUIRES(
      class StridedLayoutMapping,
      /* requires */ (
        detail::is_mapping_of<typename StridedLayoutMapping::layout_type, StridedLayoutMapping> &&
        (extents_type::rank() == StridedLayoutMapping::extents_type::rank()) &&
        StridedLayoutMapping::is_always_strided()
      )
    )
    MDSPAN_INLINE_FUNCTION
    friend constexpr bool operator!=(const mapping& x, const StridedLayoutMapping& y) noexcept {
      return !(x == y);
    }

    MDSPAN_TEMPLATE_REQUIRES(
      class OtherExtents,
      /* requires */ (
        (extents_type::rank() == OtherExtents::rank())
      )
    )
    MDSPAN_INLINE_FUNCTION
    friend constexpr bool operator!=(mapping const& lhs, mapping<OtherExtents> const& rhs) noexcept {
      return deduction_workaround_impl::_not_eq_impl(lhs, rhs);
    }
#endif

   // [mdspan.submdspan.mapping], submdspan mapping specialization
   template<class... SliceSpecifiers>
   MDSPAN_INLINE_FUNCTION
   constexpr auto submdspan_mapping_impl(
       SliceSpecifiers... slices) const;

   template<class... SliceSpecifiers>
    MDSPAN_INLINE_FUNCTION
     friend constexpr auto submdspan_mapping(
       const mapping& src, SliceSpecifiers... slices) {
      return src.submdspan_mapping_impl(slices...);
    }
  };
};

namespace detail {

template <class Layout, class Extents, class Mapping>
MDSPAN_INLINE_FUNCTION
constexpr void validate_strides(with_rank<0>, Layout, const Extents&, const Mapping&)
{}

template <std::size_t N, class Layout, class Extents, class Mapping>
MDSPAN_INLINE_FUNCTION
constexpr void validate_strides(with_rank<N>, Layout, const Extents& ext, const Mapping& other)
{
  static_assert(std::is_same<typename Mapping::layout_type, layout_stride>::value &&
                (std::is_same<Layout, layout_left>::value ||
                 std::is_same<Layout, layout_right>::value)
                , "This function is only intended to validate construction of "
                  "a layout_left or layout_right mapping from a layout_stride mapping.");

  constexpr auto is_left = std::is_same<Layout, layout_left>::value;

  typename Extents::index_type expected_stride = 1;

  for (std::size_t r = 0; r < N; r++) {
    const std::size_t s = is_left ? r : N - 1 - r;

    MDSPAN_IMPL_PRECONDITION(common_integral_compare(expected_stride, other.stride(s))
                             && "invalid strides for layout_{left,right}");

    expected_stride *= ext.extent(s);
  }
}

} // namespace detail
} // end namespace MDSPAN_IMPL_STANDARD_NAMESPACE
