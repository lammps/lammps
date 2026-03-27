// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#include <Kokkos_Macros.hpp>
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#endif

#ifndef KOKKOS_VIEW_HPP
#define KOKKOS_VIEW_HPP

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_IMPL_MDSPAN
#include <View/Kokkos_BasicView.hpp>
#endif
#ifdef KOKKOS_ENABLE_IMPL_VIEW_LEGACY
#include <View/Kokkos_ViewLegacy.hpp>
#else

#include <View/Kokkos_ViewTraits.hpp>
#include <Kokkos_MemoryTraits.hpp>

// FIXME: This will eventually be removed
namespace Kokkos::Impl {
template <class, class...>
class ViewMapping;
}
#include <View/Kokkos_ViewMapping.hpp>
#include <Kokkos_MinMax.hpp>

// Class to provide a uniform type
namespace Kokkos {
namespace Impl {
template <class ViewType, int Traits>
struct ViewUniformType;

template <class ParentView>
struct ViewTracker;
} /* namespace Impl */

template <class T1, class T2>
struct is_always_assignable_impl;

template <class... ViewTDst, class... ViewTSrc>
struct is_always_assignable_impl<Kokkos::View<ViewTDst...>,
                                 Kokkos::View<ViewTSrc...> > {
  using dst_mdspan = typename Kokkos::View<ViewTDst...>::mdspan_type;
  using src_mdspan = typename Kokkos::View<ViewTSrc...>::mdspan_type;

  constexpr static bool value =
      std::is_constructible_v<dst_mdspan, src_mdspan> &&
      static_cast<int>(Kokkos::View<ViewTDst...>::rank_dynamic) >=
          static_cast<int>(Kokkos::View<ViewTSrc...>::rank_dynamic);
};

template <class View1, class View2>
using is_always_assignable = is_always_assignable_impl<
    std::remove_reference_t<View1>,
    std::remove_const_t<std::remove_reference_t<View2> > >;

template <class T1, class T2>
inline constexpr bool is_always_assignable_v =
    is_always_assignable<T1, T2>::value;

template <class... ViewTDst, class... ViewTSrc>
constexpr bool is_assignable(const Kokkos::View<ViewTDst...>& dst,
                             const Kokkos::View<ViewTSrc...>& src) {
  using dst_mdspan = typename Kokkos::View<ViewTDst...>::mdspan_type;
  using src_mdspan = typename Kokkos::View<ViewTSrc...>::mdspan_type;

  return is_always_assignable_v<Kokkos::View<ViewTDst...>,
                                Kokkos::View<ViewTSrc...> > ||
         (std::is_constructible_v<dst_mdspan, src_mdspan> &&
          ((dst_mdspan::rank_dynamic() >= 1) ||
           (dst.static_extent(0) == src.extent(0))) &&
          ((dst_mdspan::rank_dynamic() >= 2) ||
           (dst.static_extent(1) == src.extent(1))) &&
          ((dst_mdspan::rank_dynamic() >= 3) ||
           (dst.static_extent(2) == src.extent(2))) &&
          ((dst_mdspan::rank_dynamic() >= 4) ||
           (dst.static_extent(3) == src.extent(3))) &&
          ((dst_mdspan::rank_dynamic() >= 5) ||
           (dst.static_extent(4) == src.extent(4))) &&
          ((dst_mdspan::rank_dynamic() >= 6) ||
           (dst.static_extent(5) == src.extent(5))) &&
          ((dst_mdspan::rank_dynamic() >= 7) ||
           (dst.static_extent(6) == src.extent(6))) &&
          ((dst_mdspan::rank_dynamic() == 8) ||
           (dst.static_extent(7) == src.extent(7))));
}

namespace Impl {
template <class... Properties>
struct BasicViewFromTraits {
  using view_traits        = ViewTraits<Properties...>;
  using mdspan_view_traits = MDSpanViewTraits<view_traits>;
  using element_type       = typename view_traits::value_type;
  using extents_type       = typename mdspan_view_traits::extents_type;
  using layout_type        = typename mdspan_view_traits::mdspan_layout_type;
  using accessor_type      = typename mdspan_view_traits::accessor_type;

  using type =
      BV::BasicView<element_type, extents_type, layout_type, accessor_type>;
};

// Helper function to deal with cases where the data handle is
// not convertible to element_type* such as in Sacado.
// An overload for our reference counted data handle is next to its
// implementation. This one covers Unmanaged views with raw pointers.
template <class HandleType>
KOKKOS_INLINE_FUNCTION constexpr auto ptr_from_data_handle(
    const HandleType& handle) {
  // This should only be internally invoked in Kokkos with raw pointers.
  static_assert(std::is_pointer_v<HandleType>);
  return handle;
}
}  // namespace Impl

template <class DataType, class... Properties>
struct ViewTraits;

template <class DataType, class... Properties>
class View;

template <class>
struct is_view : public std::false_type {};

template <class D, class... P>
struct is_view<View<D, P...> > : public std::true_type {};

template <class D, class... P>
struct is_view<const View<D, P...> > : public std::true_type {};

template <class T>
inline constexpr bool is_view_v = is_view<T>::value;

// FIXME spurious warnings like
// error: 'SR.14123' may be used uninitialized [-Werror=maybe-uninitialized]
#if defined(KOKKOS_COMPILER_GNU) && KOKKOS_COMPILER_GNU >= 1500
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#pragma GCC diagnostic ignored "-Wuninitialized"
#endif

template <class DataType, class... Properties>
class View : public Impl::BasicViewFromTraits<DataType, Properties...>::type {
  // We are deriving from BasicView, but need a helper to translate
  // View template parameters to BasicView template parameters
 private:
  template <class, class...>
  friend class View;
  template <typename V>
  friend struct Kokkos::Impl::ViewTracker;

  using base_t =
      typename Impl::BasicViewFromTraits<DataType, Properties...>::type;

  using base_t::m_acc;
  using base_t::m_map;
  using base_t::m_ptr;

 public:
  using base_t::base_t;

  // typedefs originally from ViewTraits
  using traits               = ViewTraits<DataType, Properties...>;
  using const_value_type     = typename traits::const_value_type;
  using non_const_value_type = typename traits::non_const_value_type;
  using data_type            = DataType;
  using const_data_type      = typename traits::const_data_type;
  using non_const_data_type  = typename traits::non_const_data_type;
  using view_tracker_type    = Impl::ViewTracker<View>;
  using array_layout         = typename traits::array_layout;
  using device_type          = typename traits::device_type;
  using execution_space      = typename traits::execution_space;
  using memory_space         = typename traits::memory_space;
  using memory_traits        = typename traits::memory_traits;
  using host_mirror_space    = typename traits::host_mirror_space;
  using typename base_t::index_type;

  // aliases from BasicView

  // FIXME: Should be unsigned
  // FIXME: these are overriden so that their types are identical when using
  // BasicView or Legacy we will need to obtain these from base_t in the future
  // and deprecate old behavior
  using size_type  = typename memory_space::size_type;
  using value_type = typename traits::value_type;
  // pointer_type can be different from element_type*
  using pointer_type = decltype(Impl::ptr_from_data_handle(
      std::declval<typename base_t::data_handle_type>()));

 private:
  using raw_allocation_value_type = std::remove_pointer_t<pointer_type>;
  using hooks_policy =
      typename Impl::ViewHooksFromTraits<DataType, Properties...>::type;
  static constexpr bool has_hooks_policy = !std::is_void_v<hooks_policy>;

 public:
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_5
  using scalar_array_type KOKKOS_DEPRECATED_WITH_COMMENT(
      "Use data_type instead.") = data_type;
  using const_scalar_array_type KOKKOS_DEPRECATED_WITH_COMMENT(
      "Use const_data_type instead.") = const_data_type;
  using non_const_scalar_array_type KOKKOS_DEPRECATED_WITH_COMMENT(
      "Use non_const_data_type instead.") = non_const_data_type;
#endif

  // typedefs from BasicView
  using typename base_t::mdspan_type;
  using reference_type = typename base_t::reference;
  using typename base_t::data_handle_type;

  //----------------------------------------
  // Compatible view of a data type
  using type = std::conditional_t<
      has_hooks_policy,
      View<typename traits::data_type, typename traits::array_layout,
           typename traits::device_type, typename traits::hooks_policy,
           typename traits::memory_traits>,
      View<typename traits::data_type, typename traits::array_layout,
           typename traits::device_type, typename traits::memory_traits> >;

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_5
  //----------------------------------------
  // Compatible view of array of scalar types
  using array_type KOKKOS_DEPRECATED_WITH_COMMENT("Use type instead.") = type;
#endif

  // Compatible view of const data type
  using const_type = std::conditional_t<
      has_hooks_policy,
      View<typename traits::const_data_type, typename traits::array_layout,
           typename traits::device_type, typename traits::hooks_policy,
           typename traits::memory_traits>,
      View<typename traits::const_data_type, typename traits::array_layout,
           typename traits::device_type, typename traits::memory_traits> >;

  // Compatible view of non-const data type
  using non_const_type = std::conditional_t<
      has_hooks_policy,
      View<typename traits::non_const_data_type, typename traits::array_layout,
           typename traits::device_type, typename traits::hooks_policy,
           typename traits::memory_traits>,
      View<typename traits::non_const_data_type, typename traits::array_layout,
           typename traits::device_type, typename traits::memory_traits> >;

  // Compatible host mirror view
  using host_mirror_type = std::conditional_t<
      has_hooks_policy,
      View<typename traits::non_const_data_type, typename traits::array_layout,
           Device<DefaultHostExecutionSpace,
                  typename traits::host_mirror_space::memory_space>,
           typename traits::hooks_policy>,
      View<typename traits::non_const_data_type, typename traits::array_layout,
           Device<DefaultHostExecutionSpace,
                  typename traits::host_mirror_space::memory_space> > >;

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  /** \brief  Compatible HostMirror view */
  using HostMirror KOKKOS_DEPRECATED_WITH_COMMENT(
      "Use host_mirror_type instead.") = host_mirror_type;
#endif

  // Unified types
  using uniform_type = typename Impl::ViewUniformType<View, 0>::type;
  using uniform_const_type =
      typename Impl::ViewUniformType<View, 0>::const_type;
  using uniform_runtime_type =
      typename Impl::ViewUniformType<View, 0>::runtime_type;
  using uniform_runtime_const_type =
      typename Impl::ViewUniformType<View, 0>::runtime_const_type;
  using uniform_nomemspace_type =
      typename Impl::ViewUniformType<View, 0>::nomemspace_type;
  using uniform_const_nomemspace_type =
      typename Impl::ViewUniformType<View, 0>::const_nomemspace_type;
  using uniform_runtime_nomemspace_type =
      typename Impl::ViewUniformType<View, 0>::runtime_nomemspace_type;
  using uniform_runtime_const_nomemspace_type =
      typename Impl::ViewUniformType<View, 0>::runtime_const_nomemspace_type;

  //----------------------------------------
  // Domain rank and extents

  static constexpr Impl::integral_constant<size_t, base_t::rank()> rank = {};
  static constexpr Impl::integral_constant<size_t, base_t::rank_dynamic()>
      rank_dynamic = {};
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  enum {Rank KOKKOS_DEPRECATED_WITH_COMMENT("Use rank instead.") = rank()};
#endif

  KOKKOS_INLINE_FUNCTION constexpr array_layout layout() const {
    return Impl::array_layout_from_mapping<array_layout, mdspan_type>(
        base_t::mapping());
  }

  // clang-format off
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
  KOKKOS_DEPRECATED_WITH_COMMENT("Use stride(0) instead") KOKKOS_FUNCTION constexpr size_t stride_0() const { return stride(0); }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use stride(1) instead") KOKKOS_FUNCTION constexpr size_t stride_1() const { return stride(1); }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use stride(2) instead") KOKKOS_FUNCTION constexpr size_t stride_2() const { return stride(2); }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use stride(3) instead") KOKKOS_FUNCTION constexpr size_t stride_3() const { return stride(3); }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use stride(4) instead") KOKKOS_FUNCTION constexpr size_t stride_4() const { return stride(4); }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use stride(5) instead") KOKKOS_FUNCTION constexpr size_t stride_5() const { return stride(5); }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use stride(6) instead") KOKKOS_FUNCTION constexpr size_t stride_6() const { return stride(6); }
  KOKKOS_DEPRECATED_WITH_COMMENT("Use stride(7) instead") KOKKOS_FUNCTION constexpr size_t stride_7() const { return stride(7); }
#endif
  // clang-format on

  template <typename iType>
  KOKKOS_INLINE_FUNCTION constexpr std::enable_if_t<std::is_integral_v<iType>,
                                                    size_t>
  stride(iType r) const {
    // base class doesn't have constraint
    // FIXME: Eventually we need to deprecate this behavior and just use
    // BasicView implementation
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
    using LayoutType = typename mdspan_type::layout_type;
    if (r >= static_cast<iType>(rank())) {
      if constexpr (rank() == 0) return 1;
      if constexpr (std::is_same_v<LayoutType, layout_right> ||
                    Impl::IsLayoutRightPadded<LayoutType>::value) {
        return 1;
      }
      if constexpr (std::is_same_v<LayoutType, layout_left> ||
                    Impl::IsLayoutLeftPadded<LayoutType>::value) {
        return base_t::stride(rank() - 1) * extent(rank() - 1);
      }
      if constexpr (std::is_same_v<LayoutType, layout_stride>) {
        return 0;
      }
    }
#else
    KOKKOS_ASSERT(r < static_cast<iType>(rank()));
#endif
    return base_t::stride(r);
  }

  template <typename iType>
  KOKKOS_INLINE_FUNCTION void stride([[maybe_unused]] iType* const s) const {
    if constexpr (rank() > 0) {
      size_t max_stride     = 0;
      size_t max_stride_idx = 0;
      for (size_t r = 0; r < rank(); r++) {
        s[r] = base_t::stride(r);
        if (s[r] > static_cast<iType>(max_stride)) {
          max_stride     = s[r];
          max_stride_idx = r;
        }
      }
      s[rank()] = max_stride * base_t::extent(max_stride_idx);
    }
  }

  //----------------------------------------
  // Range span is the span which contains all members.

  static constexpr auto reference_type_is_lvalue_reference =
      std::is_lvalue_reference_v<reference_type>;

  KOKKOS_INLINE_FUNCTION constexpr size_t span() const {
    return base_t::mapping().required_span_size();
  }
  KOKKOS_INLINE_FUNCTION bool span_is_contiguous() const {
    return base_t::is_exhaustive();
  }
  KOKKOS_INLINE_FUNCTION constexpr bool is_allocated() const {
    return data() != nullptr;
  }
  KOKKOS_INLINE_FUNCTION constexpr pointer_type data() const {
    return Impl::ptr_from_data_handle(base_t::data_handle());
  }

  KOKKOS_INLINE_FUNCTION constexpr int extent_int(size_t r) const {
    return static_cast<int>(base_t::extent(r));
  }
  //----------------------------------------
  // Allow specializations to query their specialized map

  KOKKOS_INLINE_FUNCTION
  auto impl_map() const {
    using map_type =
        Kokkos::Impl::ViewMapping<traits, typename traits::specialize>;
    return map_type(Kokkos::view_wrap(data()), layout());
  }

  KOKKOS_INLINE_FUNCTION
  const Kokkos::Impl::SharedAllocationTracker& impl_track() const {
    if constexpr (traits::memory_traits::is_unmanaged) {
      static const Kokkos::Impl::SharedAllocationTracker empty_tracker = {};
      return empty_tracker;
    } else {
      return base_t::data_handle().tracker();
    }
  }
  //----------------------------------------
  // Operators always provided by View

  template <class OtherIndexType>
  KOKKOS_FUNCTION constexpr reference_type operator[](
      const OtherIndexType& idx) const {
    return base_t::operator()(idx);
  }

 private:

#ifdef KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK
  template <typename... Is>
  static KOKKOS_FUNCTION void check_access_member_function_valid_args(
      Is... is) {
    // cast to int to work around pointless comparison of unsigned to 0 warning
    static_assert(static_cast<int>(sizeof...(Is)) <=
                  static_cast<int>(8 - rank));
    static_assert(Kokkos::Impl::are_integral<Is...>::value);
    if (!((is == static_cast<Is>(0)) && ... && true))
      Kokkos::abort("Extra arguments to Kokkos::access must be zero");
  }
#else
  template <typename... Is>
  static KOKKOS_FUNCTION void check_access_member_function_valid_args(Is...) {
    // cast to int to work around pointless comparison of unsigned to 0 warning
    static_assert(static_cast<int>(sizeof...(Is)) <=
                  static_cast<int>(8 - rank));
    static_assert(Kokkos::Impl::are_integral<Is...>::value);
  }
#endif

  // The following are shortcuts to allow implicit integer precision
  // in offset calculations - meaning index calculation happens in the common
  // type of the used indices. If and when these are not needed, the data access
  // operator should come from BasicView only. BasicView will not support this
  // directly, since BasicView anyway requires you to be explicit about the
  // desired index_type

  // ROCM: The below code segfaults the compiler with ROCM 6.3 and ROCM 6.2
  // We will simply avoid the performance optimization code path for those.
  // SYCL: The below code segfaults the compiler with Intel OneAPI 2024
  // Cuda+Clang segfaults for clang-17 andt clang-18
#if !(defined(HIP_VERSION) && HIP_VERSION_MAJOR == 6 &&                     \
      HIP_VERSION_MINOR <= 3) &&                                            \
    !(defined(KOKKOS_ENABLE_SYCL) && defined(KOKKOS_COMPILER_INTEL_LLVM) && \
      KOKKOS_COMPILER_INTEL_LLVM < 20250000) &&                             \
    !(defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOS_COMPILER_CLANG) &&      \
      KOKKOS_COMPILER_CLANG >= 1700 && KOKKOS_COMPILER_CLANG < 1900)
  // Rank 0
  KOKKOS_FUNCTION constexpr auto compute_offset(std::index_sequence<>) const {
    return 0;
  }

  // Rank 1
  template <class IndexOffset>
  KOKKOS_FUNCTION constexpr auto compute_offset(
      std::index_sequence<0>, IndexOffset index_offset) const {
    if constexpr (std::is_same_v<typename base_t::layout_type,
                                 Kokkos::layout_stride>)
      return index_offset * static_cast<IndexOffset>(m_map.stride(0));
    else
      return index_offset;
  }

  // Rank > 1
  // CUDA: using requires clauses instead of if constexpr ran into issues
  // with CUDA 12.2 + GCC 11.5, hence the if constexpr approach
  template <size_t... I, class... IndexOffsets>
  KOKKOS_FUNCTION constexpr auto compute_offset(
      std::index_sequence<I...>, IndexOffsets... index_offsets) const {
    using idx_type = std::common_type_t<IndexOffsets...>;

    if constexpr (Kokkos::Impl::IsLayoutLeftPadded<
                      typename base_t::layout_type>::value) {
      idx_type indices[] = {static_cast<idx_type>(index_offsets)...};
      // self-recursive fold trick from
      // https://github.com/llvm/llvm-project/blob/96e1914aa2e6d8966acbfbe2f4d184201f1aa318/libcxx/include/__mdspan/layout_left.h#L144
      idx_type res = 0;
      ((res = indices[rank() - 1 - I] +
              static_cast<idx_type>((rank() - 1 - I) == 0u /* extent_to_pad */
                                        ? m_map.stride(1)
                                        : extent(rank() - 1 - I)) *
                  res),
       ...);
      return res;
    } else if constexpr (Kokkos::Impl::IsLayoutRightPadded<
                             typename base_t::layout_type>::value) {
      // self-recursive fold trick from
      // https://github.com/llvm/llvm-project/blob/4d9771741d40cc9cfcccb6b033f43689d36b705a/libcxx/include/__mdspan/layout_right.h#L141
      idx_type res = 0;
      ((res = static_cast<idx_type>(index_offsets) +
              static_cast<idx_type>(I == rank() - 1 ? m_map.stride(rank() - 2)
                                                    : extent(I)) *
                  res),
       ...);
      return res;
    } else {
      return ((static_cast<idx_type>(index_offsets) *
               static_cast<idx_type>(m_map.stride(I))) +
              ...);
    }
  }

#if defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK)

#define KOKKOS_IMPL_BASICVIEW_OPERATOR_VERIFY(...)                             \
  if constexpr (Impl::IsReferenceCountedDataHandle<data_handle_type>::value) { \
    Kokkos::Impl::runtime_check_memory_access_violation<memory_space>(         \
        m_ptr.tracker());                                                      \
    Kokkos::Impl::view_verify_operator_bounds(                                 \
        m_ptr.tracker(), m_map.extents(), m_ptr.get(), __VA_ARGS__);           \
  } else {                                                                     \
    Kokkos::Impl::runtime_check_memory_access_violation<memory_space>(         \
        Kokkos::Impl::SharedAllocationTracker());                              \
    Kokkos::Impl::view_verify_operator_bounds(                                 \
        Kokkos::Impl::SharedAllocationTracker(), m_map.extents(), m_ptr,       \
        __VA_ARGS__);                                                          \
  }

#else

#define KOKKOS_IMPL_BASICVIEW_OPERATOR_VERIFY(...)                             \
  if constexpr (Impl::IsReferenceCountedDataHandle<data_handle_type>::value) { \
    Kokkos::Impl::runtime_check_memory_access_violation<memory_space>(         \
        m_ptr.tracker());                                                      \
  } else {                                                                     \
    Kokkos::Impl::runtime_check_memory_access_violation<memory_space>(         \
        Kokkos::Impl::SharedAllocationTracker());                              \
  }
#endif

 public:
  template <class... OtherIndexTypes>
    requires(
        (std::is_convertible_v<OtherIndexTypes, index_type> && ...) &&
        (std::is_nothrow_constructible_v<index_type, OtherIndexTypes> && ...) &&
        (sizeof...(OtherIndexTypes) == rank()) &&
        (Kokkos::Impl::IsLayoutLeftPadded<
             typename base_t::layout_type>::value ||
         Kokkos::Impl::IsLayoutRightPadded<
             typename base_t::layout_type>::value ||
         std::is_same_v<typename base_t::layout_type, Kokkos::layout_stride>))
  KOKKOS_FUNCTION constexpr reference_type operator()(
      OtherIndexTypes... idx) const {
    KOKKOS_IMPL_BASICVIEW_OPERATOR_VERIFY(idx...);
    return m_acc.access(
        m_ptr, static_cast<size_t>(compute_offset(
                   std::index_sequence_for<OtherIndexTypes...>{}, idx...)));
  }

  template <class... OtherIndexTypes>
    requires(
        (std::is_convertible_v<OtherIndexTypes, index_type> && ...) &&
        (std::is_nothrow_constructible_v<index_type, OtherIndexTypes> && ...) &&
        (sizeof...(OtherIndexTypes) == rank()) &&
        !(Kokkos::Impl::IsLayoutLeftPadded<
              typename base_t::layout_type>::value ||
          Kokkos::Impl::IsLayoutRightPadded<
              typename base_t::layout_type>::value ||
          std::is_same_v<typename base_t::layout_type, Kokkos::layout_stride>))
  KOKKOS_FUNCTION constexpr reference_type operator()(
      OtherIndexTypes... indices) const {
    KOKKOS_IMPL_BASICVIEW_OPERATOR_VERIFY(indices...);
    return m_acc.access(m_ptr,
                        m_map(static_cast<index_type>(std::move(indices))...));
  }
#undef KOKKOS_IMPL_BASICVIEW_OPERATOR_VERIFY
#endif

 public:
  //------------------------------
  // Rank 0

  template <typename... Is>
  KOKKOS_FORCEINLINE_FUNCTION std::enable_if_t<
      (Kokkos::Impl::always_true<Is...>::value && (0 == rank)), reference_type>
  access(Is... extra) const {
    check_access_member_function_valid_args(extra...);
    return base_t::operator()();
  }

  //------------------------------
  // Rank 1

  template <typename I0, typename... Is>
  KOKKOS_FORCEINLINE_FUNCTION std::enable_if_t<
      (Kokkos::Impl::always_true<I0, Is...>::value && (1 == rank)),
      reference_type>
  access(I0 i0, Is... extra) const {
    check_access_member_function_valid_args(extra...);
    return base_t::operator()(i0);
  }

  //------------------------------
  // Rank 2

  template <typename I0, typename I1, typename... Is>
  KOKKOS_FORCEINLINE_FUNCTION std::enable_if_t<
      (Kokkos::Impl::always_true<I0, I1, Is...>::value && (2 == rank)),
      reference_type>
  access(I0 i0, I1 i1, Is... extra) const {
    check_access_member_function_valid_args(extra...);
    return base_t::operator()(i0, i1);
  }

  //------------------------------
  // Rank 3

  template <typename I0, typename I1, typename I2, typename... Is>
  KOKKOS_FORCEINLINE_FUNCTION std::enable_if_t<
      (Kokkos::Impl::always_true<I0, I1, I2, Is...>::value && (3 == rank)),
      reference_type>
  access(I0 i0, I1 i1, I2 i2, Is... extra) const {
    check_access_member_function_valid_args(extra...);
    return base_t::operator()(i0, i1, i2);
  }

  //------------------------------
  // Rank 4

  template <typename I0, typename I1, typename I2, typename I3, typename... Is>
  KOKKOS_FORCEINLINE_FUNCTION std::enable_if_t<
      (Kokkos::Impl::always_true<I0, I1, I2, I3, Is...>::value && (4 == rank)),
      reference_type>
  access(I0 i0, I1 i1, I2 i2, I3 i3, Is... extra) const {
    check_access_member_function_valid_args(extra...);
    return base_t::operator()(i0, i1, i2, i3);
  }

  //------------------------------
  // Rank 5

  template <typename I0, typename I1, typename I2, typename I3, typename I4,
            typename... Is>
  KOKKOS_FORCEINLINE_FUNCTION std::enable_if_t<
      (Kokkos::Impl::always_true<I0, I1, I2, I3, I4, Is...>::value &&
       (5 == rank)),
      reference_type>
  access(I0 i0, I1 i1, I2 i2, I3 i3, I4 i4, Is... extra) const {
    check_access_member_function_valid_args(extra...);
    return base_t::operator()(i0, i1, i2, i3, i4);
  }

  //------------------------------
  // Rank 6

  template <typename I0, typename I1, typename I2, typename I3, typename I4,
            typename I5, typename... Is>
  KOKKOS_FORCEINLINE_FUNCTION std::enable_if_t<
      (Kokkos::Impl::always_true<I0, I1, I2, I3, I4, I5, Is...>::value &&
       (6 == rank)),
      reference_type>
  access(I0 i0, I1 i1, I2 i2, I3 i3, I4 i4, I5 i5, Is... extra) const {
    check_access_member_function_valid_args(extra...);
    return base_t::operator()(i0, i1, i2, i3, i4, i5);
  }

  //------------------------------
  // Rank 7

  template <typename I0, typename I1, typename I2, typename I3, typename I4,
            typename I5, typename I6, typename... Is>
  KOKKOS_FORCEINLINE_FUNCTION std::enable_if_t<
      (Kokkos::Impl::always_true<I0, I1, I2, I3, I4, I5, I6, Is...>::value &&
       (7 == rank)),
      reference_type>
  access(I0 i0, I1 i1, I2 i2, I3 i3, I4 i4, I5 i5, I6 i6, Is... extra) const {
    check_access_member_function_valid_args(extra...);
    return base_t::operator()(i0, i1, i2, i3, i4, i5, i6);
  }

  //------------------------------
  // Rank 8

  template <typename I0, typename I1, typename I2, typename I3, typename I4,
            typename I5, typename I6, typename I7, typename... Is>
  KOKKOS_FORCEINLINE_FUNCTION
      std::enable_if_t<(Kokkos::Impl::always_true<I0, I1, I2, I3, I4, I5, I6,
                                                  I7, Is...>::value &&
                        (8 == rank)),
                       reference_type>
      access(I0 i0, I1 i1, I2 i2, I3 i3, I4 i4, I5 i5, I6 i6, I7 i7,
             Is... extra) const {
    check_access_member_function_valid_args(extra...);
    return base_t::operator()(i0, i1, i2, i3, i4, i5, i6, i7);
  }

  //----------------------------------------
  // Standard destructor, constructors, and assignment operators

  KOKKOS_DEFAULTED_FUNCTION
  ~View() = default;

  KOKKOS_DEFAULTED_FUNCTION
  View() = default;

// FIXME_NVCC: nvcc 12.2 and 12.3 view these as ambiguous even though they have
// exclusive requirements clauses. 12.6 Also has some issues though it manifests
// differently
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_COMPILER_NVHPC)
#define KOKKOS_IMPL_VIEW_HOOKS_NVCC_WORKAROUND 1
#endif
#ifdef KOKKOS_IMPL_VIEW_HOOKS_NVCC_WORKAROUND
  KOKKOS_FUNCTION
  View(const View& other) : base_t{other} {
    if constexpr (has_hooks_policy) {
      KOKKOS_IF_ON_HOST((hooks_policy::copy_construct(*this, other);))
    }
  }
#else
  KOKKOS_DEFAULTED_FUNCTION
  View(const View&)
    requires(!has_hooks_policy)
  = default;

  KOKKOS_FUNCTION
  View(const View& other)
    requires(has_hooks_policy)
      : base_t{other} {
    KOKKOS_IF_ON_HOST((hooks_policy::copy_construct(*this, other);))
  }
#endif

#ifdef KOKKOS_IMPL_VIEW_HOOKS_NVCC_WORKAROUND
  KOKKOS_FUNCTION
  View(View&& other) : base_t{std::move(static_cast<base_t&&>(other))} {
    if constexpr (has_hooks_policy) {
      KOKKOS_IF_ON_HOST((hooks_policy::move_construct(*this, other);))
    }
  }
#else
  KOKKOS_DEFAULTED_FUNCTION
  View(View&&)
    requires(!has_hooks_policy)
  = default;

  KOKKOS_FUNCTION
  View(View&& other)
    requires(has_hooks_policy)
      : base_t{std::move(static_cast<base_t&&>(other))} {
    KOKKOS_IF_ON_HOST((hooks_policy::move_construct(*this, other);))
  }
#endif

#ifdef KOKKOS_IMPL_VIEW_HOOKS_NVCC_WORKAROUND
  KOKKOS_FUNCTION
  View& operator=(const View& other) {
    base_t::operator=(other);

    if constexpr (has_hooks_policy) {
      KOKKOS_IF_ON_HOST(
          (if (&other != this) { hooks_policy::copy_assign(*this, other); }))
    }

    return *this;
  }
#else
  KOKKOS_DEFAULTED_FUNCTION
  View& operator=(const View&)
    requires(!has_hooks_policy)
  = default;

  KOKKOS_FUNCTION
  View& operator=(const View& other)
    requires(has_hooks_policy)
  {
    base_t::operator=(other);
    KOKKOS_IF_ON_HOST(
        (if (&other != this) { hooks_policy::copy_assign(*this, other); }))

    return *this;
  }
#endif

// FIXME_NVCC: nvcc 12.2 and 12.3 view these as ambiguous even though they have
// exclusive requirements clauses. 12.6 Also has some issues though it manifests
// differently
#ifdef KOKKOS_IMPL_VIEW_HOOKS_NVCC_WORKAROUND
  KOKKOS_FUNCTION
  View& operator=(View&& other) {
    base_t::operator=(std::move(static_cast<base_t&&>(other)));

    if constexpr (has_hooks_policy) {
      KOKKOS_IF_ON_HOST(
          (if (&other != this) { hooks_policy::move_assign(*this, other); }))
    }

    return *this;
  }
#else
  KOKKOS_DEFAULTED_FUNCTION
  View& operator=(View&&)
    requires(!has_hooks_policy)
  = default;

  KOKKOS_FUNCTION
  View& operator=(View&& other)
    requires(has_hooks_policy)
  {
    base_t::operator=(std::move(static_cast<base_t&&>(other)));
    KOKKOS_IF_ON_HOST(
        (if (&other != this) { hooks_policy::move_assign(*this, other); }))

    return *this;
  }
#endif
#undef KOKKOS_IMPL_VIEW_HOOKS_NVCC_WORKAROUND

  KOKKOS_FUNCTION
  View(typename base_t::data_handle_type p,
       const typename base_t::mapping_type& m)
      : base_t(p, m) {}

  //----------------------------------------
  // Compatible view copy constructor and assignment
  // may assign unmanaged from managed.

  template <class OtherT, class... OtherArgs>
  //    requires(std::is_constructible_v<
  //             mdspan_type, typename View<OtherT, OtherArgs...>::mdspan_type>)
  KOKKOS_INLINE_FUNCTION View(
      const View<OtherT, OtherArgs...>& other,
      std::enable_if_t<
          std::is_constructible_v<
              mdspan_type, typename View<OtherT, OtherArgs...>::mdspan_type>,
          void*> = nullptr)
      : base_t([&] {
          // use an immediately invoked lambda so we can run our own checks with
          // better error messages before mdspan diagnoses problems
          base_t::check_basic_view_constructibility(other.mapping());
          return base_t(
              static_cast<typename mdspan_type::data_handle_type>(
                  other.data_handle()),
              static_cast<typename mdspan_type::mapping_type>(other.mapping()),
              static_cast<typename mdspan_type::accessor_type>(
                  other.accessor()));
        }()) {}

  //----------------------------------------
  // Compatible subview constructor
  // may assign unmanaged from managed.

  template <class RT, class... RP, class Arg0, class... Args>
  KOKKOS_INLINE_FUNCTION View(const View<RT, RP...>& src_view, const Arg0 arg0,
                              Args... args)
      : base_t(Impl::subview_ctor_tag, src_view, arg0, args...) {}

  //----------------------------------------
  // Allocation according to allocation properties and array layout

  template <class... P>
  explicit View(const Impl::ViewCtorProp<P...>& arg_prop,
                std::enable_if_t<!Impl::ViewCtorProp<P...>::has_pointer,
                                 const typename traits::array_layout&>
                    arg_layout)
      : base_t(
            arg_prop,
            Impl::mapping_from_array_layout<typename mdspan_type::mapping_type>(
                arg_layout)) {
    static_assert(!traits::memory_traits::is_unmanaged,
                  "Can't construct managed View with unmanaged memory trait!");
  }

  template <class... P>
  KOKKOS_FUNCTION explicit View(
      const Impl::ViewCtorProp<P...>& arg_prop,
      std::enable_if_t<Impl::ViewCtorProp<P...>::has_pointer,
                       const typename traits::array_layout&>
          arg_layout)
      : base_t(
            arg_prop,
            Impl::mapping_from_array_layout<typename mdspan_type::mapping_type>(
                arg_layout)) {}

  // Constructors from legacy layouts when using Views of the new layouts
  // LayoutLeft -> layout_left, layout_left_padded
  // LayoutRight -> layout_right, layout_right_padded
  // LayoutStride -> layout_stride
  KOKKOS_FUNCTION
  explicit View(const typename base_t::data_handle_type& handle,
                const LayoutStride& arg_layout)
    requires(std::is_same_v<typename base_t::layout_type, layout_stride>)
      : base_t(
            handle,
            Impl::mapping_from_array_layout<typename mdspan_type::mapping_type>(
                arg_layout)) {}

  KOKKOS_FUNCTION
  explicit View(const typename base_t::data_handle_type& handle,
                const LayoutLeft& arg_layout)
    requires(std::is_same_v<typename base_t::layout_type,
                            Experimental::layout_left_padded<> >)
      : base_t(
            handle,
            Impl::mapping_from_array_layout<typename mdspan_type::mapping_type>(
                arg_layout)) {}

  KOKKOS_FUNCTION
  explicit View(const typename base_t::data_handle_type& handle,
                const LayoutRight& arg_layout)
    requires(std::is_same_v<typename base_t::layout_type,
                            Experimental::layout_right_padded<> >)
      : base_t(
            handle,
            Impl::mapping_from_array_layout<typename mdspan_type::mapping_type>(
                arg_layout)) {}

  KOKKOS_FUNCTION
  explicit View(const typename base_t::data_handle_type& handle,
                const LayoutLeft& arg_layout)
    requires(std::is_same_v<typename base_t::layout_type, layout_left>)
      : base_t(
            handle,
            Impl::mapping_from_array_layout<typename mdspan_type::mapping_type>(
                arg_layout)) {}

  KOKKOS_FUNCTION
  explicit View(const typename base_t::data_handle_type& handle,
                const LayoutRight& arg_layout)
    requires(std::is_same_v<typename base_t::layout_type, layout_right>)
      : base_t(
            handle,
            Impl::mapping_from_array_layout<typename mdspan_type::mapping_type>(
                arg_layout)) {}

  template <class P, class... Args>
    requires(!std::is_null_pointer_v<P> &&
             std::is_constructible_v<typename base_t::data_handle_type, P> &&
             sizeof...(Args) != rank() + 1)
  KOKKOS_FUNCTION explicit View(P ptr_, Args... args)
      : View(Kokkos::view_wrap(static_cast<pointer_type>(ptr_)), args...) {}

  // Special function to be preferred over the above for string literals
  // when pointer type is char*
  template <class L, class... Args>
    requires(std::is_same_v<pointer_type, char*> &&
             std::is_same_v<const char*, L>)
  explicit View(L label, Args... args)
      : View(Kokkos::view_alloc(std::string(label)), args...) {}

  // Special function to be preferred over the above for passing in 0, NULL or
  // nullptr when pointer type is char*
  template <class... Args>
  KOKKOS_FUNCTION explicit View(std::nullptr_t, Args... args)
      : View(Kokkos::view_wrap(pointer_type(nullptr)), args...) {}

  // FIXME: Constructor which allows always 8 sizes should be deprecated
  template <class... P>
  explicit View(
      const Impl::ViewCtorProp<P...>& arg_prop,
      std::enable_if_t<!Impl::ViewCtorProp<P...>::has_pointer, const size_t>
          arg_N0          = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      const size_t arg_N1 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      const size_t arg_N2 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      const size_t arg_N3 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      const size_t arg_N4 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      const size_t arg_N5 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      const size_t arg_N6 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      const size_t arg_N7 = KOKKOS_IMPL_CTOR_DEFAULT_ARG)
      : base_t([&] {
  // use an immediately invoked lambda so we can run our own checks with
  // better error messages before mdspan diagnoses problems
#ifdef KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK
          if constexpr (std::is_same_v<typename traits::array_layout,
                                       Kokkos::LayoutLeft> ||
                        std::is_same_v<typename traits::array_layout,
                                       Kokkos::LayoutRight> ||
                        std::is_same_v<typename traits::array_layout,
                                       Kokkos::LayoutStride>) {
            auto prop_copy =
                Impl::with_properties_if_unset(arg_prop, std::string{});
            const std::string& alloc_name =
                Impl::get_property<Impl::LabelTag>(prop_copy);

            Impl::runtime_check_rank(*this, !traits::impl_is_customized, arg_N0,
                                     arg_N1, arg_N2, arg_N3, arg_N4, arg_N5,
                                     arg_N6, arg_N7, alloc_name.c_str());
          }
#endif
          return base_t(
              arg_prop,
              Impl::mapping_from_ctor_and_8sizes<
                  typename mdspan_type::mapping_type, sizeof(value_type)>(
                  arg_prop, arg_N0, arg_N1, arg_N2, arg_N3, arg_N4, arg_N5,
                  arg_N6, arg_N7));
        }()) {
    static_assert(traits::array_layout::is_extent_constructible,
                  "Layout is not constructible from extent arguments. Use "
                  "overload taking a layout object instead.");
    static_assert(!traits::memory_traits::is_unmanaged,
                  "Can't construct managed View with unmanaged memory trait!");
  }

  template <class... P>
  KOKKOS_FUNCTION explicit View(
      const Impl::ViewCtorProp<P...>& arg_prop,
      std::enable_if_t<Impl::ViewCtorProp<P...>::has_pointer, const size_t>
          arg_N0          = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      const size_t arg_N1 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      const size_t arg_N2 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      const size_t arg_N3 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      const size_t arg_N4 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      const size_t arg_N5 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      const size_t arg_N6 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      const size_t arg_N7 = KOKKOS_IMPL_CTOR_DEFAULT_ARG)
      : base_t([&] {
  // use an immediately invoked lambda so we can run our own checks with
  // better error messages before mdspan diagnoses problems
#ifdef KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK
          if constexpr (std::is_same_v<typename traits::array_layout,
                                       Kokkos::LayoutLeft> ||
                        std::is_same_v<typename traits::array_layout,
                                       Kokkos::LayoutRight> ||
                        std::is_same_v<typename traits::array_layout,
                                       Kokkos::LayoutStride>) {
            Impl::runtime_check_rank(*this, !traits::impl_is_customized, arg_N0,
                                     arg_N1, arg_N2, arg_N3, arg_N4, arg_N5,
                                     arg_N6, arg_N7, "UNMANAGED");
          }
#endif
          return base_t(
              arg_prop,
              Impl::mapping_from_ctor_and_8sizes<
                  typename mdspan_type::mapping_type, sizeof(value_type)>(
                  arg_prop, arg_N0, arg_N1, arg_N2, arg_N3, arg_N4, arg_N5,
                  arg_N6, arg_N7));
        }()) {
    static_assert(traits::array_layout::is_extent_constructible,
                  "Layout is not constructible from extent arguments. Use "
                  "overload taking a layout object instead.");
  }

  // Allocate with label and layout
  template <typename Label>
  explicit View(
      const Label& arg_label,
      std::enable_if_t<Kokkos::Impl::is_view_label<Label>::value,
                       typename traits::array_layout> const& arg_layout)
      : View(Impl::ViewCtorProp<std::string>(arg_label), arg_layout) {}

#ifdef KOKKOS_COMPILER_MSVC  // FIXME_MSVC
  // MSVC had pack expansion issues with the condition inside the enable_if
 private:
  template <class... Args>
  static constexpr bool msvc_workaround_ctor_condition_1() {
    size_t num_args = sizeof...(Args);
    bool are_constructible =
        (std::is_constructible_v<size_t, Args> && ... && true);
    return (num_args != rank() + 1) && are_constructible;
  }

 public:
#endif

  template <class... Args>
  View(std::enable_if_t<
#ifndef KOKKOS_COMPILER_MSVC
           ((sizeof...(Args)) != rank() + 1) &&
               (std::is_constructible_v<size_t, Args> && ... && true),
#else
           msvc_workaround_ctor_condition_1<Args...>(),
#endif
           const std::string&>
           arg_label,
       const Args... args)
#ifdef KOKKOS_COMPILER_INTEL_LLVM  // FIXME_INTEL
      // Eventually we want to get rid of the array_layout thing entirely.
      // For now this avoids a bug in the intel compiler 2024.2, and 2025 tested
      // that only happens with O2 or higher and makes some extents not being
      // set See https://github.com/kokkos/kokkos/pull/8202
      : View(Impl::ViewCtorProp<std::string>(arg_label),
             typename traits::array_layout(args...)) {
#else
      : View(Impl::ViewCtorProp<std::string>(arg_label), args...) {
#endif
  }

 private:
  // Special thing for Sacado taking rank()+1 integers, where the last integer
  // is the FAD dimension
  template <class... Args, size_t... Idx>
  static auto view_alloc_from_label_and_integrals(std::true_type,
                                                  const std::string& arg_label,
                                                  std::index_sequence<Idx...>,
                                                  Args... args) {
    return view_alloc(arg_label, Impl::AccessorArg_t{static_cast<size_t>(
                                     ((Idx == rank() ? args : 0) + ... + 0))});
  }

  template <class... Args, size_t... Idx>
  static auto view_alloc_from_label_and_integrals(std::false_type,
                                                  const std::string& arg_label,
                                                  std::index_sequence<Idx...>,
                                                  Args...) {
    return view_alloc(arg_label);
  }

#ifdef KOKKOS_COMPILER_MSVC  // FIXME_MSVC
  // Same as above but checking for num_args equal to rank()+1
  template <class... Args>
  static constexpr bool msvc_workaround_ctor_condition_2() {
    size_t num_args = sizeof...(Args);
    bool are_constructible =
        (std::is_constructible_v<size_t, Args> && ... && true);
    return (num_args == rank() + 1) && are_constructible;
  }
#endif

 public:
  template <class... Args>
  View(std::enable_if_t<
#ifndef KOKKOS_COMPILER_MSVC
           ((sizeof...(Args)) == rank() + 1) &&
               (std::is_constructible_v<size_t, Args> && ... && true),
#else
           msvc_workaround_ctor_condition_2<Args...>(),
#endif
           const std::string&>
           arg_label,
       const Args... args)
      : View(
            view_alloc_from_label_and_integrals(
                std::bool_constant<traits::impl_is_customized>(), arg_label,
                std::make_index_sequence<sizeof...(Args)>(), args...),
#ifdef KOKKOS_COMPILER_INTEL_LLVM  // FIXME_INTEL
            // Eventually we want to get rid of the array_layout thing entirely.
            // For now this avoids a bug in the intel compiler 2024.2, and 2025
            // tested that only happens with O2 or higher and makes some extents
            // not being set See https://github.com/kokkos/kokkos/pull/8202
            typename traits::array_layout(args...)) {
#else
            args...) {
#endif
  }

  template <class... Args>
  KOKKOS_FUNCTION View(
      std::enable_if_t<
#ifndef KOKKOS_COMPILER_MSVC
          ((sizeof...(Args)) == rank() + 1) &&
              (std::is_constructible_v<size_t, Args> && ... && true),
#else
           msvc_workaround_ctor_condition_2<Args...>(),
#endif
          const pointer_type&>
          arg_ptr,
      const Args... args)
      : View(
            Kokkos::view_wrap(arg_ptr,
                              Kokkos::Impl::AccessorArg_t{
                                  Kokkos::Array<size_t, sizeof...(Args)>{
                                      static_cast<size_t>(args)...}[rank()]}),
#ifdef KOKKOS_COMPILER_INTEL_LLVM  // FIXME_INTEL
            // Eventually we want to get rid of the array_layout thing entirely.
            // For now this avoids a bug in the intel compiler 2024.2, and 2025
            // tested that only happens with O2 or higher and makes some extents
            // not being set See https://github.com/kokkos/kokkos/pull/8202
            typename traits::array_layout(args...)) {
#else
             args...) {
#endif
  }

  //----------------------------------------
  // Memory span required to wrap these dimensions.
  KOKKOS_FUNCTION
  static constexpr size_t required_allocation_size(
      typename traits::array_layout const& layout) {
    return Impl::mapping_from_array_layout<typename base_t::mapping_type>(
               layout)
               .required_span_size() *
           sizeof(raw_allocation_value_type);
  }

 private:
  template <size_t... RankIdx, std::integral... Args>
  KOKKOS_FUNCTION static constexpr size_t impl_required_allocation_size(
      std::index_sequence<RankIdx...>, Args... args) {
    constexpr size_t num_passed_args = sizeof...(Args);
    // Deal with customized view with extra args first.
    // Secondly, handle case where the number of arguments is valid.
    // Thirdly, deal with the case where the number of arguments is
    // invalid, which the old impl allowed.
    if constexpr (traits::impl_is_customized && num_passed_args == rank() + 1) {
      size_t args_array[num_passed_args] = {static_cast<size_t>(args)...};
      size_t req_span_size =
          typename base_t::mapping_type(
              typename base_t::extents_type{args_array[RankIdx]...})
              .required_span_size();
      return req_span_size * args_array[rank()] *
             sizeof(raw_allocation_value_type);
    } else if constexpr (num_passed_args == rank_dynamic ||
                         num_passed_args == rank()) {
      size_t req_span_size =
          typename base_t::mapping_type(typename base_t::extents_type{args...})
              .required_span_size();
      return req_span_size * sizeof(typename base_t::element_type);
    }
#ifndef KOKKOS_ENABLE_DEPRECATED_CODE_5
    static_assert(
        (traits::impl_is_customized && num_passed_args == rank() + 1) ||
            num_passed_args == rank_dynamic || num_passed_args == rank(),
        "Kokkos::View::required_span_size(...) - invalid number of arguments");
#else
    else {
      size_t args_array[num_passed_args] = {static_cast<size_t>(args)...};
      size_t req_span_size =
          typename base_t::mapping_type(
              typename base_t::extents_type{args_array[RankIdx]...})
              .required_span_size();
      return req_span_size * sizeof(typename base_t::element_type);
    }
#endif
  }

 public:
  template <std::integral... Args>
  KOKKOS_FUNCTION static constexpr size_t required_allocation_size(
      Args... args) {
    static_assert(traits::array_layout::is_extent_constructible,
                  "Layout is not constructible from extent arguments. Use "
                  "overload taking a layout object instead.");
    static_assert(sizeof...(Args) == rank_dynamic || sizeof...(Args) >= rank(),
                  "Number of extents is invalid");
    return impl_required_allocation_size(std::make_index_sequence<rank()>(),
                                         args...);
  }

  //----------------------------------------
  // Shared scratch memory constructor

  static KOKKOS_INLINE_FUNCTION size_t
  shmem_size(const size_t arg_N0 = KOKKOS_INVALID_INDEX,
             const size_t arg_N1 = KOKKOS_INVALID_INDEX,
             const size_t arg_N2 = KOKKOS_INVALID_INDEX,
             const size_t arg_N3 = KOKKOS_INVALID_INDEX,
             const size_t arg_N4 = KOKKOS_INVALID_INDEX,
             const size_t arg_N5 = KOKKOS_INVALID_INDEX,
             const size_t arg_N6 = KOKKOS_INVALID_INDEX,
             const size_t arg_N7 = KOKKOS_INVALID_INDEX) {
    static_assert(traits::array_layout::is_extent_constructible,
                  "Layout is not constructible from extent arguments. Use "
                  "overload taking a layout object instead.");
    const size_t num_passed_args = Impl::count_valid_integers(
        arg_N0, arg_N1, arg_N2, arg_N3, arg_N4, arg_N5, arg_N6, arg_N7);

    // Special case to cover sacado which passes in an extra integer
    if (traits::impl_is_customized && num_passed_args == rank_dynamic + 1) {
      size_t extra_dim = 0;
      switch (rank_dynamic) {
        case 0: extra_dim = arg_N0; break;
        case 1: extra_dim = arg_N1; break;
        case 2: extra_dim = arg_N2; break;
        case 3: extra_dim = arg_N3; break;
        case 4: extra_dim = arg_N4; break;
        case 5: extra_dim = arg_N5; break;
        case 6: extra_dim = arg_N6; break;
        case 7: extra_dim = arg_N7; break;
        default:
          Kokkos::abort("This can't happen: rank_dynamic is smaller than 8");
      }
      return View::shmem_size(
          typename traits::array_layout(arg_N0, arg_N1, arg_N2, arg_N3, arg_N4,
                                        arg_N5, arg_N6, arg_N7),
          extra_dim);
    } else {
      if (num_passed_args != rank_dynamic) {
        Kokkos::abort(
            "Kokkos::View::shmem_size() rank_dynamic != number of "
            "arguments.\n");
      }
      return View::shmem_size(
          typename traits::array_layout(arg_N0, arg_N1, arg_N2, arg_N3, arg_N4,
                                        arg_N5, arg_N6, arg_N7),
          1);
    }
  }

  static KOKKOS_INLINE_FUNCTION size_t
  shmem_size(typename traits::array_layout const& arg_layout) {
    return shmem_size(arg_layout, 1);
  }

 private:
  // Want to be able to align to minimum scratch alignment or sizeof or alignof
  // elements
  static constexpr size_t scratch_value_alignment = max(
      {sizeof(raw_allocation_value_type), alignof(raw_allocation_value_type),
       static_cast<size_t>(
           traits::execution_space::scratch_memory_space::ALIGN)});

  static KOKKOS_INLINE_FUNCTION size_t shmem_size(
      typename traits::array_layout const& arg_layout, size_t extra_dim) {
    return Impl::mapping_from_array_layout<typename base_t::mapping_type>(
               arg_layout)
                   .required_span_size() *
               sizeof(raw_allocation_value_type) * extra_dim +
           scratch_value_alignment;
  }

 public:
  explicit KOKKOS_INLINE_FUNCTION View(
      const typename traits::execution_space::scratch_memory_space& arg_space,
      const typename traits::array_layout& arg_layout)
      : View(Impl::ViewCtorProp<pointer_type>(
                 static_cast<pointer_type>(arg_space.get_shmem_aligned(
                     Kokkos::Impl::ViewMapping<
                         traits,
                         typename traits::specialize>::memory_span(arg_layout),
                     scratch_value_alignment))),
             arg_layout) {}

 private:
  // Function to support use case of Trilinos Sacado where an extra dimension
  // is handed in to pass to the accessor (i.e. ensemble dimension)
  template <size_t... Idx, class... Sizes>
  KOKKOS_FUNCTION auto construct_scratch_view_from_extra_dim(
      std::index_sequence<Idx...>,
      const typename traits::execution_space::scratch_memory_space& arg_space,
      Sizes... sizes_in) const {
    size_t sizes[rank() + 1] = {static_cast<size_t>(sizes_in)...};
    const auto map =
        Impl::mapping_from_ctor_and_sizes<typename mdspan_type::mapping_type,
                                          sizeof(value_type)>(
            Kokkos::view_wrap(static_cast<pointer_type>(nullptr)),
            sizes[Idx]...);

    size_t extra_dim = sizes[rank()];
    const auto acc   = accessor_from_mapping_and_accessor_arg(
        Kokkos::Impl::AccessorTypeTag<typename base_t::accessor_type>(), map,
        Kokkos::Impl::AccessorArg_t{extra_dim});

    const size_t allocation_size =
        allocation_size_from_mapping_and_accessor(map, acc) *
        sizeof(raw_allocation_value_type);
    return base_t(static_cast<pointer_type>(arg_space.get_shmem_aligned(
                      allocation_size, scratch_value_alignment)),
                  std::move(map), std::move(acc));
  }

 public:
  // Constructor to support use case of Trilinos Sacado where an extra dimension
  // is handed in to pass to the accessor (i.e. ensemble dimension)
  // Only eligible of View customization points exists.
  template <class... Sizes>
    requires((sizeof...(Sizes) == rank() + 1) && traits::impl_is_customized)
  explicit KOKKOS_INLINE_FUNCTION View(
      const typename traits::execution_space::scratch_memory_space& arg_space,
      Sizes... sizes)
      : base_t(construct_scratch_view_from_extra_dim(
            std::make_index_sequence<rank()>(), arg_space, sizes...)) {}

  // Constructor to support cases where View is customized but no extra argument
  // is passed in. In this case the accessor is default constructed, but the
  // customization point for allocation size still needs to be used.
  template <std::integral... Sizes>
    requires((sizeof...(Sizes) == rank()) && traits::impl_is_customized)
  explicit KOKKOS_INLINE_FUNCTION View(
      const typename traits::execution_space::scratch_memory_space& arg_space,
      Sizes... sizes)
      : base_t([&] {
          const auto map = Impl::mapping_from_ctor_and_sizes<
              typename mdspan_type::mapping_type, sizeof(value_type)>(
              Kokkos::view_wrap(static_cast<pointer_type>(nullptr)), sizes...);

          const auto acc = typename base_t::accessor_type();

          const size_t allocation_size =
              allocation_size_from_mapping_and_accessor(map, acc) *
              sizeof(raw_allocation_value_type);
          return base_t(static_cast<pointer_type>(arg_space.get_shmem_aligned(
                            allocation_size, scratch_value_alignment)),
                        std::move(map), std::move(acc));
        }()) {}

  // Constructor supporting scratch view construction without customization.
  // don't need to use customization point for allocation size, and use
  // default constructor for accessor.
  explicit KOKKOS_INLINE_FUNCTION View(
      const typename traits::execution_space::scratch_memory_space& arg_space,
      const size_t arg_N0 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      const size_t arg_N1 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      const size_t arg_N2 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      const size_t arg_N3 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      const size_t arg_N4 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      const size_t arg_N5 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      const size_t arg_N6 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
      const size_t arg_N7 = KOKKOS_IMPL_CTOR_DEFAULT_ARG)
    requires(!traits::impl_is_customized)
      : base_t([&] {
          const auto map = Impl::mapping_from_ctor_and_8sizes<
              typename mdspan_type::mapping_type, sizeof(value_type)>(
              Kokkos::view_wrap(static_cast<pointer_type>(nullptr)), arg_N0,
              arg_N1, arg_N2, arg_N3, arg_N4, arg_N5, arg_N6, arg_N7);

          size_t allocation_size =
              map.required_span_size() * sizeof(value_type);
          return base_t(static_cast<pointer_type>(arg_space.get_shmem_aligned(
                            allocation_size, scratch_value_alignment)),
                        std::move(map), typename base_t::accessor_type());
        }()) {}

 public:
  //----------------------------------------
  // Allocation tracking properties
  std::string label() const {
    if constexpr (traits::memory_traits::is_unmanaged) {
      return "";
    } else {
      return this->data_handle().get_label();
    }
  }

  int use_count() const {
    if constexpr (traits::memory_traits::is_unmanaged) {
      return 0;
    } else {
      return this->data_handle().use_count();
    }
  }

  KOKKOS_FUNCTION
  constexpr typename base_t::index_type extent(size_t r) const noexcept {
    // casting to int to avoid warning for pointless comparison of unsigned
    // with 0
    if (static_cast<int>(r) >= static_cast<int>(base_t::extents_type::rank()))
      return 1;
    return base_t::extent(r);
  }

  KOKKOS_FUNCTION
  static constexpr size_t static_extent(size_t r) noexcept {
    // casting to int to avoid warning for pointless comparison of unsigned
    // with 0
    if (static_cast<int>(r) >= static_cast<int>(base_t::extents_type::rank()))
      return 1;
    size_t value = base_t::extents_type::static_extent(r);
    return value == Kokkos::dynamic_extent ? 0 : value;
  }
};

#if defined(KOKKOS_COMPILER_GNU) && KOKKOS_COMPILER_GNU >= 1500
#pragma GCC diagnostic pop
#endif

template <typename D, class... P>
KOKKOS_INLINE_FUNCTION constexpr unsigned rank(const View<D, P...>&) {
  return View<D, P...>::rank();
}

namespace Impl {

template <typename ValueType, unsigned int Rank>
struct RankDataType {
  using type = typename RankDataType<ValueType, Rank - 1>::type*;
};

template <typename ValueType>
struct RankDataType<ValueType, 0> {
  using type = ValueType;
};

template <unsigned N, typename... Args>
KOKKOS_FUNCTION std::enable_if_t<
    N == View<Args...>::rank() &&
        std::is_same_v<typename ViewTraits<Args...>::specialize, void>,
    View<Args...> >
as_view_of_rank_n(View<Args...> v) {
  return v;
}

// Placeholder implementation to compile generic code for DynRankView; should
// never be called
template <unsigned N, typename T, typename... Args>
KOKKOS_FUNCTION std::enable_if_t<
    N != View<T, Args...>::rank() &&
        std::is_same_v<typename ViewTraits<T, Args...>::specialize, void>,
    View<typename RankDataType<typename View<T, Args...>::value_type, N>::type,
         Args...> >
as_view_of_rank_n(View<T, Args...>) {
  Kokkos::abort("Trying to get at a View of the wrong rank");
  return {};
}

template <typename ViewType>
struct ApplyToViewOfStaticRank {
  template <typename Function>
  static void apply(Function&& f, ViewType a) {
    f(a);
  }
};

}  // namespace Impl
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template <class D, class... P, class... Args>
KOKKOS_INLINE_FUNCTION auto subview(const View<D, P...>& src, Args... args) {
  static_assert(View<D, P...>::rank == sizeof...(Args),
                "subview requires one argument for each source View rank");

  return typename Kokkos::Impl::ViewMapping<
      void /* deduce subview type from source view traits */
      ,
      typename Impl::RemoveAlignedMemoryTrait<D, P...>::type,
      Args...>::type(src, args...);
}

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
template <class MemoryTraits, class D, class... P, class... Args>
KOKKOS_DEPRECATED KOKKOS_INLINE_FUNCTION auto subview(const View<D, P...>& src,
                                                      Args... args) {
  static_assert(View<D, P...>::rank == sizeof...(Args),
                "subview requires one argument for each source View rank");
  static_assert(Kokkos::is_memory_traits<MemoryTraits>::value);

  return typename Kokkos::Impl::ViewMapping<
      void /* deduce subview type from source view traits */
      ,
      typename Impl::RemoveAlignedMemoryTrait<D, P..., MemoryTraits>::type,
      Args...>::type(src, args...);
}
#endif

template <class V, class... Args>
using Subview = decltype(subview(std::declval<V>(), std::declval<Args>()...));

} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

template <class LT, class... LP, class RT, class... RP>
KOKKOS_INLINE_FUNCTION bool operator==(const View<LT, LP...>& lhs,
                                       const View<RT, RP...>& rhs) {
  // Same data, layout, dimensions
  using lhs_traits = ViewTraits<LT, LP...>;
  using rhs_traits = ViewTraits<RT, RP...>;

  return std::is_same_v<typename lhs_traits::const_value_type,
                        typename rhs_traits::const_value_type> &&
         std::is_same_v<typename lhs_traits::array_layout,
                        typename rhs_traits::array_layout> &&
         std::is_same_v<typename lhs_traits::memory_space,
                        typename rhs_traits::memory_space> &&
         View<LT, LP...>::rank() == View<RT, RP...>::rank() &&
         lhs.data() == rhs.data() && lhs.span() == rhs.span() &&
         lhs.extent(0) == rhs.extent(0) && lhs.extent(1) == rhs.extent(1) &&
         lhs.extent(2) == rhs.extent(2) && lhs.extent(3) == rhs.extent(3) &&
         lhs.extent(4) == rhs.extent(4) && lhs.extent(5) == rhs.extent(5) &&
         lhs.extent(6) == rhs.extent(6) && lhs.extent(7) == rhs.extent(7);
}

template <class LT, class... LP, class RT, class... RP>
KOKKOS_INLINE_FUNCTION bool operator!=(const View<LT, LP...>& lhs,
                                       const View<RT, RP...>& rhs) {
  return !(operator==(lhs, rhs));
}

} /* namespace Kokkos */

// FIXME: https://github.com/kokkos/kokkos/issues/7736 We may want to move these
// out
#include <View/Kokkos_ViewCommonType.hpp>
#include <View/Kokkos_ViewUniformType.hpp>
#include <View/Kokkos_ViewAtomic.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* KOKKOS_ENABLE_IMPL_VIEW_LEGACY */
#endif /* #ifndef KOKKOS_VIEW_HPP */
