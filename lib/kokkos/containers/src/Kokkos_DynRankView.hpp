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

/// \file Kokkos_DynRankView.hpp
/// \brief Declaration and definition of Kokkos::DynRankView.
///
/// This header file declares and defines Kokkos::DynRankView and its
/// related nonmember functions.

#ifndef KOKKOS_DYNRANKVIEW_HPP
#define KOKKOS_DYNRANKVIEW_HPP
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_DYNRANKVIEW
#endif

#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Error.hpp>
#include <type_traits>

namespace Kokkos {

template <typename DataType, class... Properties>
class DynRankView;  // forward declare

namespace Impl {

template <typename Specialize>
struct DynRankDimTraits {
  enum : size_t { unspecified = KOKKOS_INVALID_INDEX };

  // Compute the rank of the view from the nonzero dimension arguments.
  KOKKOS_INLINE_FUNCTION
  static size_t computeRank(const size_t N0, const size_t N1, const size_t N2,
                            const size_t N3, const size_t N4, const size_t N5,
                            const size_t N6, const size_t /* N7 */) {
    return (
        (N6 == unspecified && N5 == unspecified && N4 == unspecified &&
         N3 == unspecified && N2 == unspecified && N1 == unspecified &&
         N0 == unspecified)
            ? 0
            : ((N6 == unspecified && N5 == unspecified && N4 == unspecified &&
                N3 == unspecified && N2 == unspecified && N1 == unspecified)
                   ? 1
                   : ((N6 == unspecified && N5 == unspecified &&
                       N4 == unspecified && N3 == unspecified &&
                       N2 == unspecified)
                          ? 2
                          : ((N6 == unspecified && N5 == unspecified &&
                              N4 == unspecified && N3 == unspecified)
                                 ? 3
                                 : ((N6 == unspecified && N5 == unspecified &&
                                     N4 == unspecified)
                                        ? 4
                                        : ((N6 == unspecified &&
                                            N5 == unspecified)
                                               ? 5
                                               : ((N6 == unspecified)
                                                      ? 6
                                                      : 7)))))));
  }

  // Compute the rank of the view from the nonzero layout arguments.
  template <typename Layout>
  KOKKOS_INLINE_FUNCTION static size_t computeRank(const Layout& layout) {
    return computeRank(layout.dimension[0], layout.dimension[1],
                       layout.dimension[2], layout.dimension[3],
                       layout.dimension[4], layout.dimension[5],
                       layout.dimension[6], layout.dimension[7]);
  }

  // Extra overload to match that for specialize types v2
  template <typename Layout, typename... P>
  KOKKOS_INLINE_FUNCTION static size_t computeRank(
      const Kokkos::Impl::ViewCtorProp<P...>& /* prop */,
      const Layout& layout) {
    return computeRank(layout);
  }

  // Create the layout for the rank-7 view.
  // Non-strided Layout
  template <typename Layout>
  KOKKOS_INLINE_FUNCTION static std::enable_if_t<
      (std::is_same<Layout, Kokkos::LayoutRight>::value ||
       std::is_same<Layout, Kokkos::LayoutLeft>::value),
      Layout>
  createLayout(const Layout& layout) {
    return Layout(layout.dimension[0] != unspecified ? layout.dimension[0] : 1,
                  layout.dimension[1] != unspecified ? layout.dimension[1] : 1,
                  layout.dimension[2] != unspecified ? layout.dimension[2] : 1,
                  layout.dimension[3] != unspecified ? layout.dimension[3] : 1,
                  layout.dimension[4] != unspecified ? layout.dimension[4] : 1,
                  layout.dimension[5] != unspecified ? layout.dimension[5] : 1,
                  layout.dimension[6] != unspecified ? layout.dimension[6] : 1,
                  layout.dimension[7] != unspecified ? layout.dimension[7] : 1);
  }

  // LayoutStride
  template <typename Layout>
  KOKKOS_INLINE_FUNCTION static std::enable_if_t<
      (std::is_same<Layout, Kokkos::LayoutStride>::value), Layout>
  createLayout(const Layout& layout) {
    return Layout(layout.dimension[0] != unspecified ? layout.dimension[0] : 1,
                  layout.stride[0],
                  layout.dimension[1] != unspecified ? layout.dimension[1] : 1,
                  layout.stride[1],
                  layout.dimension[2] != unspecified ? layout.dimension[2] : 1,
                  layout.stride[2],
                  layout.dimension[3] != unspecified ? layout.dimension[3] : 1,
                  layout.stride[3],
                  layout.dimension[4] != unspecified ? layout.dimension[4] : 1,
                  layout.stride[4],
                  layout.dimension[5] != unspecified ? layout.dimension[5] : 1,
                  layout.stride[5],
                  layout.dimension[6] != unspecified ? layout.dimension[6] : 1,
                  layout.stride[6],
                  layout.dimension[7] != unspecified ? layout.dimension[7] : 1,
                  layout.stride[7]);
  }

  // Extra overload to match that for specialize types
  template <typename Traits, typename... P>
  KOKKOS_INLINE_FUNCTION static std::enable_if_t<
      (std::is_same<typename Traits::array_layout,
                    Kokkos::LayoutRight>::value ||
       std::is_same<typename Traits::array_layout, Kokkos::LayoutLeft>::value ||
       std::is_same<typename Traits::array_layout,
                    Kokkos::LayoutStride>::value),
      typename Traits::array_layout>
  createLayout(const Kokkos::Impl::ViewCtorProp<P...>& /* prop */,
               const typename Traits::array_layout& layout) {
    return createLayout(layout);
  }

  // Create a view from the given dimension arguments.
  // This is only necessary because the shmem constructor doesn't take a layout.
  //   NDE shmem View's are not compatible with the added view_alloc value_type
  //   / fad_dim deduction functionality
  template <typename ViewType, typename ViewArg>
  static ViewType createView(const ViewArg& arg, const size_t N0,
                             const size_t N1, const size_t N2, const size_t N3,
                             const size_t N4, const size_t N5, const size_t N6,
                             const size_t N7) {
    return ViewType(arg, N0 != unspecified ? N0 : 1, N1 != unspecified ? N1 : 1,
                    N2 != unspecified ? N2 : 1, N3 != unspecified ? N3 : 1,
                    N4 != unspecified ? N4 : 1, N5 != unspecified ? N5 : 1,
                    N6 != unspecified ? N6 : 1, N7 != unspecified ? N7 : 1);
  }
};

// Non-strided Layout
template <typename Layout, typename iType>
KOKKOS_INLINE_FUNCTION static std::enable_if_t<
    (std::is_same<Layout, Kokkos::LayoutRight>::value ||
     std::is_same<Layout, Kokkos::LayoutLeft>::value) &&
        std::is_integral<iType>::value,
    Layout>
reconstructLayout(const Layout& layout, iType dynrank) {
  return Layout(dynrank > 0 ? layout.dimension[0] : KOKKOS_INVALID_INDEX,
                dynrank > 1 ? layout.dimension[1] : KOKKOS_INVALID_INDEX,
                dynrank > 2 ? layout.dimension[2] : KOKKOS_INVALID_INDEX,
                dynrank > 3 ? layout.dimension[3] : KOKKOS_INVALID_INDEX,
                dynrank > 4 ? layout.dimension[4] : KOKKOS_INVALID_INDEX,
                dynrank > 5 ? layout.dimension[5] : KOKKOS_INVALID_INDEX,
                dynrank > 6 ? layout.dimension[6] : KOKKOS_INVALID_INDEX,
                dynrank > 7 ? layout.dimension[7] : KOKKOS_INVALID_INDEX);
}

// LayoutStride
template <typename Layout, typename iType>
KOKKOS_INLINE_FUNCTION static std::enable_if_t<
    (std::is_same<Layout, Kokkos::LayoutStride>::value) &&
        std::is_integral<iType>::value,
    Layout>
reconstructLayout(const Layout& layout, iType dynrank) {
  return Layout(dynrank > 0 ? layout.dimension[0] : KOKKOS_INVALID_INDEX,
                dynrank > 0 ? layout.stride[0] : (0),
                dynrank > 1 ? layout.dimension[1] : KOKKOS_INVALID_INDEX,
                dynrank > 1 ? layout.stride[1] : (0),
                dynrank > 2 ? layout.dimension[2] : KOKKOS_INVALID_INDEX,
                dynrank > 2 ? layout.stride[2] : (0),
                dynrank > 3 ? layout.dimension[3] : KOKKOS_INVALID_INDEX,
                dynrank > 3 ? layout.stride[3] : (0),
                dynrank > 4 ? layout.dimension[4] : KOKKOS_INVALID_INDEX,
                dynrank > 4 ? layout.stride[4] : (0),
                dynrank > 5 ? layout.dimension[5] : KOKKOS_INVALID_INDEX,
                dynrank > 5 ? layout.stride[5] : (0),
                dynrank > 6 ? layout.dimension[6] : KOKKOS_INVALID_INDEX,
                dynrank > 6 ? layout.stride[6] : (0),
                dynrank > 7 ? layout.dimension[7] : KOKKOS_INVALID_INDEX,
                dynrank > 7 ? layout.stride[7] : (0));
}

/** \brief  Debug bounds-checking routines */
// Enhanced debug checking - most infrastructure matches that of functions in
// Kokkos_ViewMapping; additional checks for extra arguments beyond rank are 0
template <unsigned, typename iType0, class MapType>
KOKKOS_INLINE_FUNCTION bool dyn_rank_view_verify_operator_bounds(
    const iType0&, const MapType&) {
  return true;
}

template <unsigned R, typename iType0, class MapType, typename iType1,
          class... Args>
KOKKOS_INLINE_FUNCTION bool dyn_rank_view_verify_operator_bounds(
    const iType0& rank, const MapType& map, const iType1& i, Args... args) {
  if (static_cast<iType0>(R) < rank) {
    return (size_t(i) < map.extent(R)) &&
           dyn_rank_view_verify_operator_bounds<R + 1>(rank, map, args...);
  } else if (i != 0) {
    Kokkos::printf(
        "DynRankView Debug Bounds Checking Error: at rank %u\n  Extra "
        "arguments beyond the rank must be zero \n",
        R);
    return (false) &&
           dyn_rank_view_verify_operator_bounds<R + 1>(rank, map, args...);
  } else {
    return (true) &&
           dyn_rank_view_verify_operator_bounds<R + 1>(rank, map, args...);
  }
}

template <unsigned, class MapType>
inline void dyn_rank_view_error_operator_bounds(char*, int, const MapType&) {}

template <unsigned R, class MapType, class iType, class... Args>
inline void dyn_rank_view_error_operator_bounds(char* buf, int len,
                                                const MapType& map,
                                                const iType& i, Args... args) {
  const int n = snprintf(
      buf, len, " %ld < %ld %c", static_cast<unsigned long>(i),
      static_cast<unsigned long>(map.extent(R)), (sizeof...(Args) ? ',' : ')'));
  dyn_rank_view_error_operator_bounds<R + 1>(buf + n, len - n, map, args...);
}

// op_rank = rank of the operator version that was called
template <typename MemorySpace, typename iType0, typename iType1, class MapType,
          class... Args>
KOKKOS_INLINE_FUNCTION void dyn_rank_view_verify_operator_bounds(
    const iType0& op_rank, const iType1& rank,
    const Kokkos::Impl::SharedAllocationTracker& tracker, const MapType& map,
    Args... args) {
  if (static_cast<iType0>(rank) > op_rank) {
    Kokkos::abort(
        "DynRankView Bounds Checking Error: Need at least rank arguments to "
        "the operator()");
  }

  if (!dyn_rank_view_verify_operator_bounds<0>(rank, map, args...)) {
    KOKKOS_IF_ON_HOST(
        (enum {LEN = 1024}; char buffer[LEN];
         const std::string label = tracker.template get_label<MemorySpace>();
         int n = snprintf(buffer, LEN, "DynRankView bounds error of view %s (",
                          label.c_str());
         dyn_rank_view_error_operator_bounds<0>(buffer + n, LEN - n, map,
                                                args...);
         Kokkos::Impl::throw_runtime_exception(std::string(buffer));))

    KOKKOS_IF_ON_DEVICE(
        ((void)tracker; Kokkos::abort("DynRankView bounds error");))
  }
}

/** \brief  Assign compatible default mappings */
struct ViewToDynRankViewTag {};

}  // namespace Impl

namespace Impl {

template <class DstTraits, class SrcTraits>
class ViewMapping<
    DstTraits, SrcTraits,
    std::enable_if_t<(std::is_same<typename DstTraits::memory_space,
                                   typename SrcTraits::memory_space>::value &&
                      std::is_void<typename DstTraits::specialize>::value &&
                      std::is_void<typename SrcTraits::specialize>::value &&
                      (std::is_same<typename DstTraits::array_layout,
                                    typename SrcTraits::array_layout>::value ||
                       ((std::is_same<typename DstTraits::array_layout,
                                      Kokkos::LayoutLeft>::value ||
                         std::is_same<typename DstTraits::array_layout,
                                      Kokkos::LayoutRight>::value ||
                         std::is_same<typename DstTraits::array_layout,
                                      Kokkos::LayoutStride>::value) &&
                        (std::is_same<typename SrcTraits::array_layout,
                                      Kokkos::LayoutLeft>::value ||
                         std::is_same<typename SrcTraits::array_layout,
                                      Kokkos::LayoutRight>::value ||
                         std::is_same<typename SrcTraits::array_layout,
                                      Kokkos::LayoutStride>::value)))),
                     Kokkos::Impl::ViewToDynRankViewTag>> {
 private:
  enum {
    is_assignable_value_type =
        std::is_same<typename DstTraits::value_type,
                     typename SrcTraits::value_type>::value ||
        std::is_same<typename DstTraits::value_type,
                     typename SrcTraits::const_value_type>::value
  };

  enum {
    is_assignable_layout =
        std::is_same<typename DstTraits::array_layout,
                     typename SrcTraits::array_layout>::value ||
        std::is_same<typename DstTraits::array_layout,
                     Kokkos::LayoutStride>::value
  };

 public:
  enum { is_assignable = is_assignable_value_type && is_assignable_layout };

  using DstType = ViewMapping<DstTraits, typename DstTraits::specialize>;
  using SrcType = ViewMapping<SrcTraits, typename SrcTraits::specialize>;

  template <typename DT, typename... DP, typename ST, typename... SP>
  KOKKOS_INLINE_FUNCTION static void assign(
      Kokkos::DynRankView<DT, DP...>& dst, const Kokkos::View<ST, SP...>& src) {
    static_assert(
        is_assignable_value_type,
        "View assignment must have same value type or const = non-const");

    static_assert(
        is_assignable_layout,
        "View assignment must have compatible layout or have rank <= 1");

    // Removed dimension checks...

    using dst_offset_type   = typename DstType::offset_type;
    dst.m_map.m_impl_offset = dst_offset_type(
        std::integral_constant<unsigned, 0>(),
        src.layout());  // Check this for integer input1 for padding, etc
    dst.m_map.m_impl_handle = Kokkos::Impl::ViewDataHandle<DstTraits>::assign(
        src.m_map.m_impl_handle, src.m_track.m_tracker);
    dst.m_track.assign(src.m_track.m_tracker, DstTraits::is_managed);
    dst.m_rank = Kokkos::View<ST, SP...>::rank();
  }
};

}  // namespace Impl

/* \class DynRankView
 * \brief Container that creates a Kokkos view with rank determined at runtime.
 *   Essentially this is a rank 7 view
 *
 *   Changes from View
 *   1. The rank of the DynRankView is returned by the method rank()
 *   2. Max rank of a DynRankView is 7
 *   3. subview called with 'subview(...)' or 'subdynrankview(...)' (backward
 * compatibility)
 *   4. Every subview is returned with LayoutStride
 *   5. Copy and Copy-Assign View to DynRankView
 *   6. deep_copy between Views and DynRankViews
 *   7. rank( view ); returns the rank of View or DynRankView
 *
 */

template <class>
struct is_dyn_rank_view : public std::false_type {};

template <class D, class... P>
struct is_dyn_rank_view<Kokkos::DynRankView<D, P...>> : public std::true_type {
};

template <class T>
inline constexpr bool is_dyn_rank_view_v = is_dyn_rank_view<T>::value;

template <typename DataType, class... Properties>
class DynRankView : public ViewTraits<DataType, Properties...> {
  static_assert(!std::is_array<DataType>::value &&
                    !std::is_pointer<DataType>::value,
                "Cannot template DynRankView with array or pointer datatype - "
                "must be pod");

 private:
  template <class, class...>
  friend class DynRankView;
  template <class, class...>
  friend class Kokkos::Impl::ViewMapping;

 public:
  using drvtraits = ViewTraits<DataType, Properties...>;

  using view_type = View<DataType*******, Properties...>;

  using traits = ViewTraits<DataType*******, Properties...>;

 private:
  using map_type =
      Kokkos::Impl::ViewMapping<traits, typename traits::specialize>;
  using track_type = Kokkos::Impl::SharedAllocationTracker;

  track_type m_track;
  map_type m_map;
  unsigned m_rank;

 public:
  KOKKOS_INLINE_FUNCTION
  view_type& DownCast() const { return (view_type&)(*this); }
  KOKKOS_INLINE_FUNCTION
  const view_type& ConstDownCast() const { return (const view_type&)(*this); }

  // Types below - at least the HostMirror requires the value_type, NOT the rank
  // 7 data_type of the traits

  /** \brief  Compatible view of array of scalar types */
  using array_type = DynRankView<
      typename drvtraits::scalar_array_type, typename drvtraits::array_layout,
      typename drvtraits::device_type, typename drvtraits::memory_traits>;

  /** \brief  Compatible view of const data type */
  using const_type = DynRankView<
      typename drvtraits::const_data_type, typename drvtraits::array_layout,
      typename drvtraits::device_type, typename drvtraits::memory_traits>;

  /** \brief  Compatible view of non-const data type */
  using non_const_type = DynRankView<
      typename drvtraits::non_const_data_type, typename drvtraits::array_layout,
      typename drvtraits::device_type, typename drvtraits::memory_traits>;

  /** \brief  Compatible HostMirror view */
  using HostMirror = DynRankView<typename drvtraits::non_const_data_type,
                                 typename drvtraits::array_layout,
                                 typename drvtraits::host_mirror_space>;

  //----------------------------------------
  // Domain rank and extents

  //  enum { Rank = map_type::Rank }; //Will be dyn rank of 7 always, keep the
  //  enum?

  template <typename iType>
  KOKKOS_INLINE_FUNCTION constexpr std::enable_if_t<
      std::is_integral<iType>::value, size_t>
  extent(const iType& r) const {
    return m_map.extent(r);
  }

  template <typename iType>
  KOKKOS_INLINE_FUNCTION constexpr std::enable_if_t<
      std::is_integral<iType>::value, int>
  extent_int(const iType& r) const {
    return static_cast<int>(m_map.extent(r));
  }

  KOKKOS_INLINE_FUNCTION constexpr typename traits::array_layout layout() const;

  //----------------------------------------
  /*  Deprecate all 'dimension' functions in favor of
   *  ISO/C++ vocabulary 'extent'.
   */

  KOKKOS_INLINE_FUNCTION constexpr size_t size() const {
    return m_map.extent(0) * m_map.extent(1) * m_map.extent(2) *
           m_map.extent(3) * m_map.extent(4) * m_map.extent(5) *
           m_map.extent(6) * m_map.extent(7);
  }

  KOKKOS_INLINE_FUNCTION constexpr size_t stride_0() const {
    return m_map.stride_0();
  }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_1() const {
    return m_map.stride_1();
  }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_2() const {
    return m_map.stride_2();
  }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_3() const {
    return m_map.stride_3();
  }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_4() const {
    return m_map.stride_4();
  }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_5() const {
    return m_map.stride_5();
  }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_6() const {
    return m_map.stride_6();
  }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_7() const {
    return m_map.stride_7();
  }

  template <typename iType>
  KOKKOS_INLINE_FUNCTION void stride(iType* const s) const {
    m_map.stride(s);
  }

  //----------------------------------------
  // Range span is the span which contains all members.

  using reference_type = typename map_type::reference_type;
  using pointer_type   = typename map_type::pointer_type;

  enum {
    reference_type_is_lvalue_reference =
        std::is_lvalue_reference<reference_type>::value
  };

  KOKKOS_INLINE_FUNCTION constexpr size_t span() const { return m_map.span(); }
  KOKKOS_INLINE_FUNCTION constexpr bool span_is_contiguous() const {
    return m_map.span_is_contiguous();
  }
  KOKKOS_INLINE_FUNCTION constexpr pointer_type data() const {
    return m_map.data();
  }
  KOKKOS_INLINE_FUNCTION constexpr bool is_allocated() const {
    return (m_map.data() != nullptr);
  }

  //----------------------------------------
  // Allow specializations to query their specialized map
  KOKKOS_INLINE_FUNCTION
  const Kokkos::Impl::ViewMapping<traits, typename traits::specialize>&
  impl_map() const {
    return m_map;
  }

  //----------------------------------------

 private:
  enum {
    is_layout_left =
        std::is_same<typename traits::array_layout, Kokkos::LayoutLeft>::value,

    is_layout_right =
        std::is_same<typename traits::array_layout, Kokkos::LayoutRight>::value,

    is_layout_stride = std::is_same<typename traits::array_layout,
                                    Kokkos::LayoutStride>::value,

    is_default_map = std::is_void<typename traits::specialize>::value &&
                     (is_layout_left || is_layout_right || is_layout_stride)
  };

// Bounds checking macros
#if defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK)

// rank of the calling operator - included as first argument in ARG
#define KOKKOS_IMPL_VIEW_OPERATOR_VERIFY(ARG)                             \
  Kokkos::Impl::runtime_check_memory_access_violation<                    \
      typename traits::memory_space>(                                     \
      "Kokkos::DynRankView ERROR: attempt to access inaccessible memory " \
      "space");                                                           \
  Kokkos::Impl::dyn_rank_view_verify_operator_bounds<                     \
      typename traits::memory_space>                                      \
      ARG;

#else

#define KOKKOS_IMPL_VIEW_OPERATOR_VERIFY(ARG)                             \
  Kokkos::Impl::runtime_check_memory_access_violation<                    \
      typename traits::memory_space>(                                     \
      "Kokkos::DynRankView ERROR: attempt to access inaccessible memory " \
      "space");

#endif

 public:
  KOKKOS_INLINE_FUNCTION
  constexpr unsigned rank() const { return m_rank; }

  // operators ()
  // Rank 0
  KOKKOS_INLINE_FUNCTION
  reference_type operator()() const {
    KOKKOS_IMPL_VIEW_OPERATOR_VERIFY((0, this->rank(), m_track, m_map))
    return impl_map().reference();
    // return m_map.reference(0,0,0,0,0,0,0);
  }

  // Rank 1
  // This assumes a contiguous underlying memory (i.e. no padding, no
  // striding...)
  template <typename iType>
  KOKKOS_INLINE_FUNCTION std::enable_if_t<
      std::is_same<typename drvtraits::value_type,
                   typename drvtraits::scalar_array_type>::value &&
          std::is_integral<iType>::value,
      reference_type>
  operator[](const iType& i0) const {
    // Phalanx is violating this, since they use the operator to access ALL
    // elements in the allocation KOKKOS_IMPL_VIEW_OPERATOR_VERIFY( (1 ,
    // this->rank(), m_track, m_map) )
    return data()[i0];
  }

  // This assumes a contiguous underlying memory (i.e. no padding, no
  // striding... AND a Trilinos/Sacado scalar type )
  template <typename iType>
  KOKKOS_INLINE_FUNCTION std::enable_if_t<
      !std::is_same<typename drvtraits::value_type,
                    typename drvtraits::scalar_array_type>::value &&
          std::is_integral<iType>::value,
      reference_type>
  operator[](const iType& i0) const {
    //      auto map = impl_map();
    const size_t dim_scalar = m_map.dimension_scalar();
    const size_t bytes      = this->span() / dim_scalar;

    using tmp_view_type = Kokkos::View<
        DataType*, typename traits::array_layout, typename traits::device_type,
        Kokkos::MemoryTraits<traits::memory_traits::is_unmanaged |
                             traits::memory_traits::is_random_access |
                             traits::memory_traits::is_atomic>>;
    tmp_view_type rankone_view(this->data(), bytes, dim_scalar);
    return rankone_view(i0);
  }

  // Rank 1 parenthesis
  template <typename iType>
  KOKKOS_INLINE_FUNCTION
      std::enable_if_t<(std::is_void<typename traits::specialize>::value &&
                        std::is_integral<iType>::value),
                       reference_type>
      operator()(const iType& i0) const {
    KOKKOS_IMPL_VIEW_OPERATOR_VERIFY((1, this->rank(), m_track, m_map, i0))
    return m_map.reference(i0);
  }

  template <typename iType>
  KOKKOS_INLINE_FUNCTION
      std::enable_if_t<!(std::is_void<typename traits::specialize>::value &&
                         std::is_integral<iType>::value),
                       reference_type>
      operator()(const iType& i0) const {
    KOKKOS_IMPL_VIEW_OPERATOR_VERIFY((1, this->rank(), m_track, m_map, i0))
    return m_map.reference(i0, 0, 0, 0, 0, 0, 0);
  }

  // Rank 2
  template <typename iType0, typename iType1>
  KOKKOS_INLINE_FUNCTION std::enable_if_t<
      (std::is_void<typename traits::specialize>::value &&
       std::is_integral<iType0>::value && std::is_integral<iType1>::value),
      reference_type>
  operator()(const iType0& i0, const iType1& i1) const {
    KOKKOS_IMPL_VIEW_OPERATOR_VERIFY((2, this->rank(), m_track, m_map, i0, i1))
    return m_map.reference(i0, i1);
  }

  template <typename iType0, typename iType1>
  KOKKOS_INLINE_FUNCTION
      std::enable_if_t<!(std::is_void<typename drvtraits::specialize>::value &&
                         std::is_integral<iType0>::value),
                       reference_type>
      operator()(const iType0& i0, const iType1& i1) const {
    KOKKOS_IMPL_VIEW_OPERATOR_VERIFY((2, this->rank(), m_track, m_map, i0, i1))
    return m_map.reference(i0, i1, 0, 0, 0, 0, 0);
  }

  // Rank 3
  template <typename iType0, typename iType1, typename iType2>
  KOKKOS_INLINE_FUNCTION std::enable_if_t<
      (std::is_void<typename traits::specialize>::value &&
       std::is_integral<iType0>::value && std::is_integral<iType1>::value &&
       std::is_integral<iType2>::value),
      reference_type>
  operator()(const iType0& i0, const iType1& i1, const iType2& i2) const {
    KOKKOS_IMPL_VIEW_OPERATOR_VERIFY(
        (3, this->rank(), m_track, m_map, i0, i1, i2))
    return m_map.reference(i0, i1, i2);
  }

  template <typename iType0, typename iType1, typename iType2>
  KOKKOS_INLINE_FUNCTION
      std::enable_if_t<!(std::is_void<typename drvtraits::specialize>::value &&
                         std::is_integral<iType0>::value),
                       reference_type>
      operator()(const iType0& i0, const iType1& i1, const iType2& i2) const {
    KOKKOS_IMPL_VIEW_OPERATOR_VERIFY(
        (3, this->rank(), m_track, m_map, i0, i1, i2))
    return m_map.reference(i0, i1, i2, 0, 0, 0, 0);
  }

  // Rank 4
  template <typename iType0, typename iType1, typename iType2, typename iType3>
  KOKKOS_INLINE_FUNCTION std::enable_if_t<
      (std::is_void<typename traits::specialize>::value &&
       std::is_integral<iType0>::value && std::is_integral<iType1>::value &&
       std::is_integral<iType2>::value && std::is_integral<iType3>::value),
      reference_type>
  operator()(const iType0& i0, const iType1& i1, const iType2& i2,
             const iType3& i3) const {
    KOKKOS_IMPL_VIEW_OPERATOR_VERIFY(
        (4, this->rank(), m_track, m_map, i0, i1, i2, i3))
    return m_map.reference(i0, i1, i2, i3);
  }

  template <typename iType0, typename iType1, typename iType2, typename iType3>
  KOKKOS_INLINE_FUNCTION
      std::enable_if_t<!(std::is_void<typename drvtraits::specialize>::value &&
                         std::is_integral<iType0>::value),
                       reference_type>
      operator()(const iType0& i0, const iType1& i1, const iType2& i2,
                 const iType3& i3) const {
    KOKKOS_IMPL_VIEW_OPERATOR_VERIFY(
        (4, this->rank(), m_track, m_map, i0, i1, i2, i3))
    return m_map.reference(i0, i1, i2, i3, 0, 0, 0);
  }

  // Rank 5
  template <typename iType0, typename iType1, typename iType2, typename iType3,
            typename iType4>
  KOKKOS_INLINE_FUNCTION std::enable_if_t<
      (std::is_void<typename traits::specialize>::value &&
       std::is_integral<iType0>::value && std::is_integral<iType1>::value &&
       std::is_integral<iType2>::value && std::is_integral<iType3>::value &&
       std::is_integral<iType4>::value),
      reference_type>
  operator()(const iType0& i0, const iType1& i1, const iType2& i2,
             const iType3& i3, const iType4& i4) const {
    KOKKOS_IMPL_VIEW_OPERATOR_VERIFY(
        (5, this->rank(), m_track, m_map, i0, i1, i2, i3, i4))
    return m_map.reference(i0, i1, i2, i3, i4);
  }

  template <typename iType0, typename iType1, typename iType2, typename iType3,
            typename iType4>
  KOKKOS_INLINE_FUNCTION
      std::enable_if_t<!(std::is_void<typename drvtraits::specialize>::value &&
                         std::is_integral<iType0>::value),
                       reference_type>
      operator()(const iType0& i0, const iType1& i1, const iType2& i2,
                 const iType3& i3, const iType4& i4) const {
    KOKKOS_IMPL_VIEW_OPERATOR_VERIFY(
        (5, this->rank(), m_track, m_map, i0, i1, i2, i3, i4))
    return m_map.reference(i0, i1, i2, i3, i4, 0, 0);
  }

  // Rank 6
  template <typename iType0, typename iType1, typename iType2, typename iType3,
            typename iType4, typename iType5>
  KOKKOS_INLINE_FUNCTION std::enable_if_t<
      (std::is_void<typename traits::specialize>::value &&
       std::is_integral<iType0>::value && std::is_integral<iType1>::value &&
       std::is_integral<iType2>::value && std::is_integral<iType3>::value &&
       std::is_integral<iType4>::value && std::is_integral<iType5>::value),
      reference_type>
  operator()(const iType0& i0, const iType1& i1, const iType2& i2,
             const iType3& i3, const iType4& i4, const iType5& i5) const {
    KOKKOS_IMPL_VIEW_OPERATOR_VERIFY(
        (6, this->rank(), m_track, m_map, i0, i1, i2, i3, i4, i5))
    return m_map.reference(i0, i1, i2, i3, i4, i5);
  }

  template <typename iType0, typename iType1, typename iType2, typename iType3,
            typename iType4, typename iType5>
  KOKKOS_INLINE_FUNCTION
      std::enable_if_t<!(std::is_void<typename drvtraits::specialize>::value &&
                         std::is_integral<iType0>::value),
                       reference_type>
      operator()(const iType0& i0, const iType1& i1, const iType2& i2,
                 const iType3& i3, const iType4& i4, const iType5& i5) const {
    KOKKOS_IMPL_VIEW_OPERATOR_VERIFY(
        (6, this->rank(), m_track, m_map, i0, i1, i2, i3, i4, i5))
    return m_map.reference(i0, i1, i2, i3, i4, i5, 0);
  }

  // Rank 7
  template <typename iType0, typename iType1, typename iType2, typename iType3,
            typename iType4, typename iType5, typename iType6>
  KOKKOS_INLINE_FUNCTION std::enable_if_t<
      (std::is_integral<iType0>::value && std::is_integral<iType1>::value &&
       std::is_integral<iType2>::value && std::is_integral<iType3>::value &&
       std::is_integral<iType4>::value && std::is_integral<iType5>::value &&
       std::is_integral<iType6>::value),
      reference_type>
  operator()(const iType0& i0, const iType1& i1, const iType2& i2,
             const iType3& i3, const iType4& i4, const iType5& i5,
             const iType6& i6) const {
    KOKKOS_IMPL_VIEW_OPERATOR_VERIFY(
        (7, this->rank(), m_track, m_map, i0, i1, i2, i3, i4, i5, i6))
    return m_map.reference(i0, i1, i2, i3, i4, i5, i6);
  }

  // Rank 0
  KOKKOS_INLINE_FUNCTION
  reference_type access() const {
    KOKKOS_IMPL_VIEW_OPERATOR_VERIFY((0, this->rank(), m_track, m_map))
    return impl_map().reference();
    // return m_map.reference(0,0,0,0,0,0,0);
  }

  // Rank 1
  // Rank 1 parenthesis
  template <typename iType>
  KOKKOS_INLINE_FUNCTION
      std::enable_if_t<(std::is_void<typename traits::specialize>::value &&
                        std::is_integral<iType>::value),
                       reference_type>
      access(const iType& i0) const {
    KOKKOS_IMPL_VIEW_OPERATOR_VERIFY((1, this->rank(), m_track, m_map, i0))
    return m_map.reference(i0);
  }

  template <typename iType>
  KOKKOS_INLINE_FUNCTION
      std::enable_if_t<!(std::is_void<typename traits::specialize>::value &&
                         std::is_integral<iType>::value),
                       reference_type>
      access(const iType& i0) const {
    KOKKOS_IMPL_VIEW_OPERATOR_VERIFY((1, this->rank(), m_track, m_map, i0))
    return m_map.reference(i0, 0, 0, 0, 0, 0, 0);
  }

  // Rank 2
  template <typename iType0, typename iType1>
  KOKKOS_INLINE_FUNCTION std::enable_if_t<
      (std::is_void<typename traits::specialize>::value &&
       std::is_integral<iType0>::value && std::is_integral<iType1>::value),
      reference_type>
  access(const iType0& i0, const iType1& i1) const {
    KOKKOS_IMPL_VIEW_OPERATOR_VERIFY((2, this->rank(), m_track, m_map, i0, i1))
    return m_map.reference(i0, i1);
  }

  template <typename iType0, typename iType1>
  KOKKOS_INLINE_FUNCTION
      std::enable_if_t<!(std::is_void<typename drvtraits::specialize>::value &&
                         std::is_integral<iType0>::value),
                       reference_type>
      access(const iType0& i0, const iType1& i1) const {
    KOKKOS_IMPL_VIEW_OPERATOR_VERIFY((2, this->rank(), m_track, m_map, i0, i1))
    return m_map.reference(i0, i1, 0, 0, 0, 0, 0);
  }

  // Rank 3
  template <typename iType0, typename iType1, typename iType2>
  KOKKOS_INLINE_FUNCTION std::enable_if_t<
      (std::is_void<typename traits::specialize>::value &&
       std::is_integral<iType0>::value && std::is_integral<iType1>::value &&
       std::is_integral<iType2>::value),
      reference_type>
  access(const iType0& i0, const iType1& i1, const iType2& i2) const {
    KOKKOS_IMPL_VIEW_OPERATOR_VERIFY(
        (3, this->rank(), m_track, m_map, i0, i1, i2))
    return m_map.reference(i0, i1, i2);
  }

  template <typename iType0, typename iType1, typename iType2>
  KOKKOS_INLINE_FUNCTION
      std::enable_if_t<!(std::is_void<typename drvtraits::specialize>::value &&
                         std::is_integral<iType0>::value),
                       reference_type>
      access(const iType0& i0, const iType1& i1, const iType2& i2) const {
    KOKKOS_IMPL_VIEW_OPERATOR_VERIFY(
        (3, this->rank(), m_track, m_map, i0, i1, i2))
    return m_map.reference(i0, i1, i2, 0, 0, 0, 0);
  }

  // Rank 4
  template <typename iType0, typename iType1, typename iType2, typename iType3>
  KOKKOS_INLINE_FUNCTION std::enable_if_t<
      (std::is_void<typename traits::specialize>::value &&
       std::is_integral<iType0>::value && std::is_integral<iType1>::value &&
       std::is_integral<iType2>::value && std::is_integral<iType3>::value),
      reference_type>
  access(const iType0& i0, const iType1& i1, const iType2& i2,
         const iType3& i3) const {
    KOKKOS_IMPL_VIEW_OPERATOR_VERIFY(
        (4, this->rank(), m_track, m_map, i0, i1, i2, i3))
    return m_map.reference(i0, i1, i2, i3);
  }

  template <typename iType0, typename iType1, typename iType2, typename iType3>
  KOKKOS_INLINE_FUNCTION
      std::enable_if_t<!(std::is_void<typename drvtraits::specialize>::value &&
                         std::is_integral<iType0>::value),
                       reference_type>
      access(const iType0& i0, const iType1& i1, const iType2& i2,
             const iType3& i3) const {
    KOKKOS_IMPL_VIEW_OPERATOR_VERIFY(
        (4, this->rank(), m_track, m_map, i0, i1, i2, i3))
    return m_map.reference(i0, i1, i2, i3, 0, 0, 0);
  }

  // Rank 5
  template <typename iType0, typename iType1, typename iType2, typename iType3,
            typename iType4>
  KOKKOS_INLINE_FUNCTION std::enable_if_t<
      (std::is_void<typename traits::specialize>::value &&
       std::is_integral<iType0>::value && std::is_integral<iType1>::value &&
       std::is_integral<iType2>::value && std::is_integral<iType3>::value &&
       std::is_integral<iType4>::value),
      reference_type>
  access(const iType0& i0, const iType1& i1, const iType2& i2, const iType3& i3,
         const iType4& i4) const {
    KOKKOS_IMPL_VIEW_OPERATOR_VERIFY(
        (5, this->rank(), m_track, m_map, i0, i1, i2, i3, i4))
    return m_map.reference(i0, i1, i2, i3, i4);
  }

  template <typename iType0, typename iType1, typename iType2, typename iType3,
            typename iType4>
  KOKKOS_INLINE_FUNCTION
      std::enable_if_t<!(std::is_void<typename drvtraits::specialize>::value &&
                         std::is_integral<iType0>::value),
                       reference_type>
      access(const iType0& i0, const iType1& i1, const iType2& i2,
             const iType3& i3, const iType4& i4) const {
    KOKKOS_IMPL_VIEW_OPERATOR_VERIFY(
        (5, this->rank(), m_track, m_map, i0, i1, i2, i3, i4))
    return m_map.reference(i0, i1, i2, i3, i4, 0, 0);
  }

  // Rank 6
  template <typename iType0, typename iType1, typename iType2, typename iType3,
            typename iType4, typename iType5>
  KOKKOS_INLINE_FUNCTION std::enable_if_t<
      (std::is_void<typename traits::specialize>::value &&
       std::is_integral<iType0>::value && std::is_integral<iType1>::value &&
       std::is_integral<iType2>::value && std::is_integral<iType3>::value &&
       std::is_integral<iType4>::value && std::is_integral<iType5>::value),
      reference_type>
  access(const iType0& i0, const iType1& i1, const iType2& i2, const iType3& i3,
         const iType4& i4, const iType5& i5) const {
    KOKKOS_IMPL_VIEW_OPERATOR_VERIFY(
        (6, this->rank(), m_track, m_map, i0, i1, i2, i3, i4, i5))
    return m_map.reference(i0, i1, i2, i3, i4, i5);
  }

  template <typename iType0, typename iType1, typename iType2, typename iType3,
            typename iType4, typename iType5>
  KOKKOS_INLINE_FUNCTION
      std::enable_if_t<!(std::is_void<typename drvtraits::specialize>::value &&
                         std::is_integral<iType0>::value),
                       reference_type>
      access(const iType0& i0, const iType1& i1, const iType2& i2,
             const iType3& i3, const iType4& i4, const iType5& i5) const {
    KOKKOS_IMPL_VIEW_OPERATOR_VERIFY(
        (6, this->rank(), m_track, m_map, i0, i1, i2, i3, i4, i5))
    return m_map.reference(i0, i1, i2, i3, i4, i5, 0);
  }

  // Rank 7
  template <typename iType0, typename iType1, typename iType2, typename iType3,
            typename iType4, typename iType5, typename iType6>
  KOKKOS_INLINE_FUNCTION std::enable_if_t<
      (std::is_integral<iType0>::value && std::is_integral<iType1>::value &&
       std::is_integral<iType2>::value && std::is_integral<iType3>::value &&
       std::is_integral<iType4>::value && std::is_integral<iType5>::value &&
       std::is_integral<iType6>::value),
      reference_type>
  access(const iType0& i0, const iType1& i1, const iType2& i2, const iType3& i3,
         const iType4& i4, const iType5& i5, const iType6& i6) const {
    KOKKOS_IMPL_VIEW_OPERATOR_VERIFY(
        (7, this->rank(), m_track, m_map, i0, i1, i2, i3, i4, i5, i6))
    return m_map.reference(i0, i1, i2, i3, i4, i5, i6);
  }

#undef KOKKOS_IMPL_VIEW_OPERATOR_VERIFY

  //----------------------------------------
  // Standard constructor, destructor, and assignment operators...

  KOKKOS_DEFAULTED_FUNCTION
  ~DynRankView() = default;

  KOKKOS_INLINE_FUNCTION
  DynRankView() : m_track(), m_map(), m_rank() {}  // Default ctor

  KOKKOS_INLINE_FUNCTION
  DynRankView(const DynRankView& rhs)
      : m_track(rhs.m_track), m_map(rhs.m_map), m_rank(rhs.m_rank) {}

  KOKKOS_INLINE_FUNCTION
  DynRankView(DynRankView&& rhs)
      : m_track(rhs.m_track), m_map(rhs.m_map), m_rank(rhs.m_rank) {}

  KOKKOS_INLINE_FUNCTION
  DynRankView& operator=(const DynRankView& rhs) {
    m_track = rhs.m_track;
    m_map   = rhs.m_map;
    m_rank  = rhs.m_rank;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  DynRankView& operator=(DynRankView&& rhs) {
    m_track = rhs.m_track;
    m_map   = rhs.m_map;
    m_rank  = rhs.m_rank;
    return *this;
  }

  //----------------------------------------
  // Compatible view copy constructor and assignment
  // may assign unmanaged from managed.
  template <class RT, class... RP>
  KOKKOS_INLINE_FUNCTION DynRankView(const DynRankView<RT, RP...>& rhs)
      : m_track(rhs.m_track, traits::is_managed), m_map(), m_rank(rhs.m_rank) {
    using SrcTraits = typename DynRankView<RT, RP...>::traits;
    using Mapping   = Kokkos::Impl::ViewMapping<traits, SrcTraits,
                                              typename traits::specialize>;
    static_assert(Mapping::is_assignable,
                  "Incompatible DynRankView copy construction");
    Mapping::assign(m_map, rhs.m_map, rhs.m_track);
  }

  template <class RT, class... RP>
  KOKKOS_INLINE_FUNCTION DynRankView& operator=(
      const DynRankView<RT, RP...>& rhs) {
    using SrcTraits = typename DynRankView<RT, RP...>::traits;
    using Mapping   = Kokkos::Impl::ViewMapping<traits, SrcTraits,
                                              typename traits::specialize>;
    static_assert(Mapping::is_assignable,
                  "Incompatible DynRankView copy construction");
    Mapping::assign(m_map, rhs.m_map, rhs.m_track);
    m_track.assign(rhs.m_track, traits::is_managed);
    m_rank = rhs.rank();
    return *this;
  }

  // Copy/Assign View to DynRankView
  template <class RT, class... RP>
  KOKKOS_INLINE_FUNCTION DynRankView(const View<RT, RP...>& rhs)
      : m_track(), m_map(), m_rank(View<RT, RP...>::rank()) {
    using SrcTraits = typename View<RT, RP...>::traits;
    using Mapping =
        Kokkos::Impl::ViewMapping<traits, SrcTraits,
                                  Kokkos::Impl::ViewToDynRankViewTag>;
    static_assert(Mapping::is_assignable,
                  "Incompatible View to DynRankView copy construction");
    Mapping::assign(*this, rhs);
  }

  template <class RT, class... RP>
  KOKKOS_INLINE_FUNCTION DynRankView& operator=(const View<RT, RP...>& rhs) {
    using SrcTraits = typename View<RT, RP...>::traits;
    using Mapping =
        Kokkos::Impl::ViewMapping<traits, SrcTraits,
                                  Kokkos::Impl::ViewToDynRankViewTag>;
    static_assert(Mapping::is_assignable,
                  "Incompatible View to DynRankView copy assignment");
    Mapping::assign(*this, rhs);
    return *this;
  }

  //----------------------------------------
  // Allocation tracking properties

  KOKKOS_INLINE_FUNCTION
  int use_count() const { return m_track.use_count(); }

  inline const std::string label() const {
    return m_track.template get_label<typename traits::memory_space>();
  }

  //----------------------------------------
  // Allocation according to allocation properties and array layout
  // unused arg_layout dimensions must be set to KOKKOS_INVALID_INDEX so that
  // rank deduction can properly take place
  template <class... P>
  explicit inline DynRankView(
      const Kokkos::Impl::ViewCtorProp<P...>& arg_prop,
      std::enable_if_t<!Kokkos::Impl::ViewCtorProp<P...>::has_pointer,
                       typename traits::array_layout> const& arg_layout)
      : m_track(),
        m_map(),
        m_rank(Impl::DynRankDimTraits<typename traits::specialize>::
                   template computeRank<typename traits::array_layout, P...>(
                       arg_prop, arg_layout)) {
    // Copy the input allocation properties with possibly defaulted properties
    auto prop_copy = Impl::with_properties_if_unset(
        arg_prop, std::string{}, typename traits::device_type::memory_space{},
        typename traits::device_type::execution_space{});
    using alloc_prop = decltype(prop_copy);

    static_assert(traits::is_managed,
                  "View allocation constructor requires managed memory");

    if (alloc_prop::initialize &&
        !alloc_prop::execution_space::impl_is_initialized()) {
      // If initializing view data then
      // the execution space must be initialized.
      Kokkos::Impl::throw_runtime_exception(
          "Constructing DynRankView and initializing data with uninitialized "
          "execution space");
    }

    Kokkos::Impl::SharedAllocationRecord<>* record = m_map.allocate_shared(
        prop_copy,
        Impl::DynRankDimTraits<typename traits::specialize>::
            template createLayout<traits, P...>(arg_prop, arg_layout),
        Impl::ViewCtorProp<P...>::has_execution_space);

    // Setup and initialization complete, start tracking
    m_track.assign_allocated_record_to_uninitialized(record);
  }

  // Wrappers
  template <class... P>
  explicit KOKKOS_INLINE_FUNCTION DynRankView(
      const Kokkos::Impl::ViewCtorProp<P...>& arg_prop,
      std::enable_if_t<Kokkos::Impl::ViewCtorProp<P...>::has_pointer,
                       typename traits::array_layout> const& arg_layout)
      : m_track()  // No memory tracking
        ,
        m_map(arg_prop,
              Impl::DynRankDimTraits<typename traits::specialize>::
                  template createLayout<traits, P...>(arg_prop, arg_layout)),
        m_rank(Impl::DynRankDimTraits<typename traits::specialize>::
                   template computeRank<typename traits::array_layout, P...>(
                       arg_prop, arg_layout)) {
    static_assert(
        std::is_same<pointer_type,
                     typename Impl::ViewCtorProp<P...>::pointer_type>::value,
        "Constructing DynRankView to wrap user memory must supply matching "
        "pointer type");
  }

  //----------------------------------------
  // Constructor(s)

  // Simple dimension-only layout
  template <class... P>
  explicit inline DynRankView(
      const Kokkos::Impl::ViewCtorProp<P...>& arg_prop,
      std::enable_if_t<!Kokkos::Impl::ViewCtorProp<P...>::has_pointer,
                       size_t> const arg_N0 = KOKKOS_INVALID_INDEX,
      const size_t arg_N1                   = KOKKOS_INVALID_INDEX,
      const size_t arg_N2                   = KOKKOS_INVALID_INDEX,
      const size_t arg_N3                   = KOKKOS_INVALID_INDEX,
      const size_t arg_N4                   = KOKKOS_INVALID_INDEX,
      const size_t arg_N5                   = KOKKOS_INVALID_INDEX,
      const size_t arg_N6                   = KOKKOS_INVALID_INDEX,
      const size_t arg_N7                   = KOKKOS_INVALID_INDEX)
      : DynRankView(arg_prop, typename traits::array_layout(
                                  arg_N0, arg_N1, arg_N2, arg_N3, arg_N4,
                                  arg_N5, arg_N6, arg_N7)) {}

  template <class... P>
  explicit KOKKOS_INLINE_FUNCTION DynRankView(
      const Kokkos::Impl::ViewCtorProp<P...>& arg_prop,
      std::enable_if_t<Kokkos::Impl::ViewCtorProp<P...>::has_pointer,
                       size_t> const arg_N0 = KOKKOS_INVALID_INDEX,
      const size_t arg_N1                   = KOKKOS_INVALID_INDEX,
      const size_t arg_N2                   = KOKKOS_INVALID_INDEX,
      const size_t arg_N3                   = KOKKOS_INVALID_INDEX,
      const size_t arg_N4                   = KOKKOS_INVALID_INDEX,
      const size_t arg_N5                   = KOKKOS_INVALID_INDEX,
      const size_t arg_N6                   = KOKKOS_INVALID_INDEX,
      const size_t arg_N7                   = KOKKOS_INVALID_INDEX)
      : DynRankView(arg_prop, typename traits::array_layout(
                                  arg_N0, arg_N1, arg_N2, arg_N3, arg_N4,
                                  arg_N5, arg_N6, arg_N7)) {}

  // Allocate with label and layout
  template <typename Label>
  explicit inline DynRankView(
      const Label& arg_label,
      std::enable_if_t<Kokkos::Impl::is_view_label<Label>::value,
                       typename traits::array_layout> const& arg_layout)
      : DynRankView(Kokkos::Impl::ViewCtorProp<std::string>(arg_label),
                    arg_layout) {}

  // Allocate label and layout, must disambiguate from subview constructor
  template <typename Label>
  explicit inline DynRankView(
      const Label& arg_label,
      std::enable_if_t<Kokkos::Impl::is_view_label<Label>::value, const size_t>
          arg_N0          = KOKKOS_INVALID_INDEX,
      const size_t arg_N1 = KOKKOS_INVALID_INDEX,
      const size_t arg_N2 = KOKKOS_INVALID_INDEX,
      const size_t arg_N3 = KOKKOS_INVALID_INDEX,
      const size_t arg_N4 = KOKKOS_INVALID_INDEX,
      const size_t arg_N5 = KOKKOS_INVALID_INDEX,
      const size_t arg_N6 = KOKKOS_INVALID_INDEX,
      const size_t arg_N7 = KOKKOS_INVALID_INDEX)
      : DynRankView(
            Kokkos::Impl::ViewCtorProp<std::string>(arg_label),
            typename traits::array_layout(arg_N0, arg_N1, arg_N2, arg_N3,
                                          arg_N4, arg_N5, arg_N6, arg_N7)) {}

  //----------------------------------------
  // Memory span required to wrap these dimensions.
  static constexpr size_t required_allocation_size(
      const size_t arg_N0 = 0, const size_t arg_N1 = 0, const size_t arg_N2 = 0,
      const size_t arg_N3 = 0, const size_t arg_N4 = 0, const size_t arg_N5 = 0,
      const size_t arg_N6 = 0, const size_t arg_N7 = 0) {
    return map_type::memory_span(typename traits::array_layout(
        arg_N0, arg_N1, arg_N2, arg_N3, arg_N4, arg_N5, arg_N6, arg_N7));
  }

  explicit KOKKOS_INLINE_FUNCTION DynRankView(
      pointer_type arg_ptr, const size_t arg_N0 = KOKKOS_INVALID_INDEX,
      const size_t arg_N1 = KOKKOS_INVALID_INDEX,
      const size_t arg_N2 = KOKKOS_INVALID_INDEX,
      const size_t arg_N3 = KOKKOS_INVALID_INDEX,
      const size_t arg_N4 = KOKKOS_INVALID_INDEX,
      const size_t arg_N5 = KOKKOS_INVALID_INDEX,
      const size_t arg_N6 = KOKKOS_INVALID_INDEX,
      const size_t arg_N7 = KOKKOS_INVALID_INDEX)
      : DynRankView(Kokkos::Impl::ViewCtorProp<pointer_type>(arg_ptr), arg_N0,
                    arg_N1, arg_N2, arg_N3, arg_N4, arg_N5, arg_N6, arg_N7) {}

  explicit KOKKOS_INLINE_FUNCTION DynRankView(
      pointer_type arg_ptr, typename traits::array_layout& arg_layout)
      : DynRankView(Kokkos::Impl::ViewCtorProp<pointer_type>(arg_ptr),
                    arg_layout) {}

  //----------------------------------------
  // Shared scratch memory constructor

  static inline size_t shmem_size(const size_t arg_N0 = KOKKOS_INVALID_INDEX,
                                  const size_t arg_N1 = KOKKOS_INVALID_INDEX,
                                  const size_t arg_N2 = KOKKOS_INVALID_INDEX,
                                  const size_t arg_N3 = KOKKOS_INVALID_INDEX,
                                  const size_t arg_N4 = KOKKOS_INVALID_INDEX,
                                  const size_t arg_N5 = KOKKOS_INVALID_INDEX,
                                  const size_t arg_N6 = KOKKOS_INVALID_INDEX,
                                  const size_t arg_N7 = KOKKOS_INVALID_INDEX) {
    const size_t num_passed_args =
        (arg_N0 != KOKKOS_INVALID_INDEX) + (arg_N1 != KOKKOS_INVALID_INDEX) +
        (arg_N2 != KOKKOS_INVALID_INDEX) + (arg_N3 != KOKKOS_INVALID_INDEX) +
        (arg_N4 != KOKKOS_INVALID_INDEX) + (arg_N5 != KOKKOS_INVALID_INDEX) +
        (arg_N6 != KOKKOS_INVALID_INDEX) + (arg_N7 != KOKKOS_INVALID_INDEX);

    if (std::is_void<typename traits::specialize>::value &&
        num_passed_args != traits::rank_dynamic) {
      Kokkos::abort(
          "Kokkos::View::shmem_size() rank_dynamic != number of arguments.\n");
    }
    {}

    return map_type::memory_span(typename traits::array_layout(
        arg_N0, arg_N1, arg_N2, arg_N3, arg_N4, arg_N5, arg_N6, arg_N7));
  }

  explicit KOKKOS_INLINE_FUNCTION DynRankView(
      const typename traits::execution_space::scratch_memory_space& arg_space,
      const typename traits::array_layout& arg_layout)
      : DynRankView(
            Kokkos::Impl::ViewCtorProp<pointer_type>(
                reinterpret_cast<pointer_type>(
                    arg_space.get_shmem(map_type::memory_span(
                        Impl::DynRankDimTraits<typename traits::specialize>::
                            createLayout(arg_layout)  // is this correct?
                        )))),
            arg_layout) {}

  explicit KOKKOS_INLINE_FUNCTION DynRankView(
      const typename traits::execution_space::scratch_memory_space& arg_space,
      const size_t arg_N0 = KOKKOS_INVALID_INDEX,
      const size_t arg_N1 = KOKKOS_INVALID_INDEX,
      const size_t arg_N2 = KOKKOS_INVALID_INDEX,
      const size_t arg_N3 = KOKKOS_INVALID_INDEX,
      const size_t arg_N4 = KOKKOS_INVALID_INDEX,
      const size_t arg_N5 = KOKKOS_INVALID_INDEX,
      const size_t arg_N6 = KOKKOS_INVALID_INDEX,
      const size_t arg_N7 = KOKKOS_INVALID_INDEX)

      : DynRankView(
            Kokkos::Impl::ViewCtorProp<pointer_type>(
                reinterpret_cast<pointer_type>(
                    arg_space.get_shmem(map_type::memory_span(
                        Impl::DynRankDimTraits<typename traits::specialize>::
                            createLayout(typename traits::array_layout(
                                arg_N0, arg_N1, arg_N2, arg_N3, arg_N4, arg_N5,
                                arg_N6, arg_N7)))))),
            typename traits::array_layout(arg_N0, arg_N1, arg_N2, arg_N3,
                                          arg_N4, arg_N5, arg_N6, arg_N7)) {}
};

template <typename D, class... P>
KOKKOS_INLINE_FUNCTION constexpr unsigned rank(
    const DynRankView<D, P...>& DRV) {
  return DRV.rank();
}  // needed for transition to common constexpr method in view and dynrankview
   // to return rank

//----------------------------------------------------------------------------
// Subview mapping.
// Deduce destination view type from source view traits and subview arguments

namespace Impl {

struct DynRankSubviewTag {};

}  // namespace Impl

namespace Impl {

template <class SrcTraits, class... Args>
class ViewMapping<
    std::enable_if_t<(std::is_void<typename SrcTraits::specialize>::value &&
                      (std::is_same<typename SrcTraits::array_layout,
                                    Kokkos::LayoutLeft>::value ||
                       std::is_same<typename SrcTraits::array_layout,
                                    Kokkos::LayoutRight>::value ||
                       std::is_same<typename SrcTraits::array_layout,
                                    Kokkos::LayoutStride>::value)),
                     Kokkos::Impl::DynRankSubviewTag>,
    SrcTraits, Args...> {
 private:
  enum {
    RZ = false,
    R0 = bool(is_integral_extent<0, Args...>::value),
    R1 = bool(is_integral_extent<1, Args...>::value),
    R2 = bool(is_integral_extent<2, Args...>::value),
    R3 = bool(is_integral_extent<3, Args...>::value),
    R4 = bool(is_integral_extent<4, Args...>::value),
    R5 = bool(is_integral_extent<5, Args...>::value),
    R6 = bool(is_integral_extent<6, Args...>::value)
  };

  enum {
    rank = unsigned(R0) + unsigned(R1) + unsigned(R2) + unsigned(R3) +
           unsigned(R4) + unsigned(R5) + unsigned(R6)
  };

  using array_layout = Kokkos::LayoutStride;

  using value_type = typename SrcTraits::value_type;

  using data_type = value_type*******;

 public:
  using traits_type = Kokkos::ViewTraits<data_type, array_layout,
                                         typename SrcTraits::device_type,
                                         typename SrcTraits::memory_traits>;

  using type =
      Kokkos::View<data_type, array_layout, typename SrcTraits::device_type,
                   typename SrcTraits::memory_traits>;

  template <class MemoryTraits>
  struct apply {
    static_assert(Kokkos::is_memory_traits<MemoryTraits>::value);

    using traits_type =
        Kokkos::ViewTraits<data_type, array_layout,
                           typename SrcTraits::device_type, MemoryTraits>;

    using type = Kokkos::View<data_type, array_layout,
                              typename SrcTraits::device_type, MemoryTraits>;
  };

  using dimension = typename SrcTraits::dimension;

  template <class Arg0 = int, class Arg1 = int, class Arg2 = int,
            class Arg3 = int, class Arg4 = int, class Arg5 = int,
            class Arg6 = int>
  struct ExtentGenerator {
    KOKKOS_INLINE_FUNCTION
    static SubviewExtents<7, rank> generator(
        const dimension& dim, Arg0 arg0 = Arg0(), Arg1 arg1 = Arg1(),
        Arg2 arg2 = Arg2(), Arg3 arg3 = Arg3(), Arg4 arg4 = Arg4(),
        Arg5 arg5 = Arg5(), Arg6 arg6 = Arg6()) {
      return SubviewExtents<7, rank>(dim, arg0, arg1, arg2, arg3, arg4, arg5,
                                     arg6);
    }
  };

  using ret_type = Kokkos::DynRankView<value_type, array_layout,
                                       typename SrcTraits::device_type,
                                       typename SrcTraits::memory_traits>;

  template <typename T, class... P>
  KOKKOS_INLINE_FUNCTION static ret_type subview(
      const unsigned src_rank, Kokkos::DynRankView<T, P...> const& src,
      Args... args) {
    using DstType = ViewMapping<traits_type, typename traits_type::specialize>;

    using DstDimType = std::conditional_t<
        (rank == 0), ViewDimension<>,
        std::conditional_t<
            (rank == 1), ViewDimension<0>,
            std::conditional_t<
                (rank == 2), ViewDimension<0, 0>,
                std::conditional_t<
                    (rank == 3), ViewDimension<0, 0, 0>,
                    std::conditional_t<
                        (rank == 4), ViewDimension<0, 0, 0, 0>,
                        std::conditional_t<
                            (rank == 5), ViewDimension<0, 0, 0, 0, 0>,
                            std::conditional_t<
                                (rank == 6), ViewDimension<0, 0, 0, 0, 0, 0>,
                                ViewDimension<0, 0, 0, 0, 0, 0, 0>>>>>>>>;

    using dst_offset_type = ViewOffset<DstDimType, Kokkos::LayoutStride>;
    using dst_handle_type = typename DstType::handle_type;

    ret_type dst;

    const SubviewExtents<7, rank> extents = ExtentGenerator<Args...>::generator(
        src.m_map.m_impl_offset.m_dim, args...);

    dst_offset_type tempdst(src.m_map.m_impl_offset, extents);

    dst.m_track = src.m_track;

    dst.m_map.m_impl_offset.m_dim.N0 = tempdst.m_dim.N0;
    dst.m_map.m_impl_offset.m_dim.N1 = tempdst.m_dim.N1;
    dst.m_map.m_impl_offset.m_dim.N2 = tempdst.m_dim.N2;
    dst.m_map.m_impl_offset.m_dim.N3 = tempdst.m_dim.N3;
    dst.m_map.m_impl_offset.m_dim.N4 = tempdst.m_dim.N4;
    dst.m_map.m_impl_offset.m_dim.N5 = tempdst.m_dim.N5;
    dst.m_map.m_impl_offset.m_dim.N6 = tempdst.m_dim.N6;

    dst.m_map.m_impl_offset.m_stride.S0 = tempdst.m_stride.S0;
    dst.m_map.m_impl_offset.m_stride.S1 = tempdst.m_stride.S1;
    dst.m_map.m_impl_offset.m_stride.S2 = tempdst.m_stride.S2;
    dst.m_map.m_impl_offset.m_stride.S3 = tempdst.m_stride.S3;
    dst.m_map.m_impl_offset.m_stride.S4 = tempdst.m_stride.S4;
    dst.m_map.m_impl_offset.m_stride.S5 = tempdst.m_stride.S5;
    dst.m_map.m_impl_offset.m_stride.S6 = tempdst.m_stride.S6;

    dst.m_map.m_impl_handle =
        dst_handle_type(src.m_map.m_impl_handle +
                        src.m_map.m_impl_offset(
                            extents.domain_offset(0), extents.domain_offset(1),
                            extents.domain_offset(2), extents.domain_offset(3),
                            extents.domain_offset(4), extents.domain_offset(5),
                            extents.domain_offset(6)));

    dst.m_rank =
        (src_rank > 0 ? unsigned(R0) : 0) + (src_rank > 1 ? unsigned(R1) : 0) +
        (src_rank > 2 ? unsigned(R2) : 0) + (src_rank > 3 ? unsigned(R3) : 0) +
        (src_rank > 4 ? unsigned(R4) : 0) + (src_rank > 5 ? unsigned(R5) : 0) +
        (src_rank > 6 ? unsigned(R6) : 0);

    return dst;
  }
};

}  // namespace Impl

template <class V, class... Args>
using Subdynrankview =
    typename Kokkos::Impl::ViewMapping<Kokkos::Impl::DynRankSubviewTag, V,
                                       Args...>::ret_type;

template <class D, class... P, class... Args>
KOKKOS_INLINE_FUNCTION Subdynrankview<ViewTraits<D*******, P...>, Args...>
subdynrankview(const Kokkos::DynRankView<D, P...>& src, Args... args) {
  if (src.rank() > sizeof...(Args))  // allow sizeof...(Args) >= src.rank(),
                                     // ignore the remaining args
  {
    Kokkos::abort(
        "subdynrankview: num of args must be >= rank of the source "
        "DynRankView");
  }

  using metafcn =
      Kokkos::Impl::ViewMapping<Kokkos::Impl::DynRankSubviewTag,
                                Kokkos::ViewTraits<D*******, P...>, Args...>;

  return metafcn::subview(src.rank(), src, args...);
}

// Wrapper to allow subview function name
template <class D, class... P, class... Args>
KOKKOS_INLINE_FUNCTION Subdynrankview<ViewTraits<D*******, P...>, Args...>
subview(const Kokkos::DynRankView<D, P...>& src, Args... args) {
  return subdynrankview(src, args...);
}

}  // namespace Kokkos

namespace Kokkos {

// overload == and !=
template <class LT, class... LP, class RT, class... RP>
KOKKOS_INLINE_FUNCTION bool operator==(const DynRankView<LT, LP...>& lhs,
                                       const DynRankView<RT, RP...>& rhs) {
  // Same data, layout, dimensions
  using lhs_traits = ViewTraits<LT, LP...>;
  using rhs_traits = ViewTraits<RT, RP...>;

  return std::is_same<typename lhs_traits::const_value_type,
                      typename rhs_traits::const_value_type>::value &&
         std::is_same<typename lhs_traits::array_layout,
                      typename rhs_traits::array_layout>::value &&
         std::is_same<typename lhs_traits::memory_space,
                      typename rhs_traits::memory_space>::value &&
         lhs.rank() == rhs.rank() && lhs.data() == rhs.data() &&
         lhs.span() == rhs.span() && lhs.extent(0) == rhs.extent(0) &&
         lhs.extent(1) == rhs.extent(1) && lhs.extent(2) == rhs.extent(2) &&
         lhs.extent(3) == rhs.extent(3) && lhs.extent(4) == rhs.extent(4) &&
         lhs.extent(5) == rhs.extent(5) && lhs.extent(6) == rhs.extent(6) &&
         lhs.extent(7) == rhs.extent(7);
}

template <class LT, class... LP, class RT, class... RP>
KOKKOS_INLINE_FUNCTION bool operator!=(const DynRankView<LT, LP...>& lhs,
                                       const DynRankView<RT, RP...>& rhs) {
  return !(operator==(lhs, rhs));
}

}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
namespace Kokkos {
namespace Impl {

template <class OutputView, class Enable = void>
struct DynRankViewFill {
  using const_value_type = typename OutputView::traits::const_value_type;

  const OutputView output;
  const_value_type input;

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t i0) const {
    const size_t n1 = output.extent(1);
    const size_t n2 = output.extent(2);
    const size_t n3 = output.extent(3);
    const size_t n4 = output.extent(4);
    const size_t n5 = output.extent(5);
    const size_t n6 = output.extent(6);

    for (size_t i1 = 0; i1 < n1; ++i1) {
      for (size_t i2 = 0; i2 < n2; ++i2) {
        for (size_t i3 = 0; i3 < n3; ++i3) {
          for (size_t i4 = 0; i4 < n4; ++i4) {
            for (size_t i5 = 0; i5 < n5; ++i5) {
              for (size_t i6 = 0; i6 < n6; ++i6) {
                output.access(i0, i1, i2, i3, i4, i5, i6) = input;
              }
            }
          }
        }
      }
    }
  }

  DynRankViewFill(const OutputView& arg_out, const_value_type& arg_in)
      : output(arg_out), input(arg_in) {
    using execution_space = typename OutputView::execution_space;
    using Policy          = Kokkos::RangePolicy<execution_space>;

    Kokkos::parallel_for("Kokkos::DynRankViewFill", Policy(0, output.extent(0)),
                         *this);
  }
};

template <class OutputView>
struct DynRankViewFill<OutputView, std::enable_if_t<OutputView::rank == 0>> {
  DynRankViewFill(const OutputView& dst,
                  const typename OutputView::const_value_type& src) {
    Kokkos::Impl::DeepCopy<typename OutputView::memory_space,
                           Kokkos::HostSpace>(
        dst.data(), &src, sizeof(typename OutputView::const_value_type));
  }
};

template <class OutputView, class InputView,
          class ExecSpace = typename OutputView::execution_space>
struct DynRankViewRemap {
  const OutputView output;
  const InputView input;
  const size_t n0;
  const size_t n1;
  const size_t n2;
  const size_t n3;
  const size_t n4;
  const size_t n5;
  const size_t n6;
  const size_t n7;

  DynRankViewRemap(const ExecSpace& exec_space, const OutputView& arg_out,
                   const InputView& arg_in)
      : output(arg_out),
        input(arg_in),
        n0(std::min((size_t)arg_out.extent(0), (size_t)arg_in.extent(0))),
        n1(std::min((size_t)arg_out.extent(1), (size_t)arg_in.extent(1))),
        n2(std::min((size_t)arg_out.extent(2), (size_t)arg_in.extent(2))),
        n3(std::min((size_t)arg_out.extent(3), (size_t)arg_in.extent(3))),
        n4(std::min((size_t)arg_out.extent(4), (size_t)arg_in.extent(4))),
        n5(std::min((size_t)arg_out.extent(5), (size_t)arg_in.extent(5))),
        n6(std::min((size_t)arg_out.extent(6), (size_t)arg_in.extent(6))),
        n7(std::min((size_t)arg_out.extent(7), (size_t)arg_in.extent(7))) {
    using Policy = Kokkos::RangePolicy<ExecSpace>;

    Kokkos::parallel_for("Kokkos::DynRankViewRemap", Policy(exec_space, 0, n0),
                         *this);
  }

  DynRankViewRemap(const OutputView& arg_out, const InputView& arg_in)
      : output(arg_out),
        input(arg_in),
        n0(std::min((size_t)arg_out.extent(0), (size_t)arg_in.extent(0))),
        n1(std::min((size_t)arg_out.extent(1), (size_t)arg_in.extent(1))),
        n2(std::min((size_t)arg_out.extent(2), (size_t)arg_in.extent(2))),
        n3(std::min((size_t)arg_out.extent(3), (size_t)arg_in.extent(3))),
        n4(std::min((size_t)arg_out.extent(4), (size_t)arg_in.extent(4))),
        n5(std::min((size_t)arg_out.extent(5), (size_t)arg_in.extent(5))),
        n6(std::min((size_t)arg_out.extent(6), (size_t)arg_in.extent(6))),
        n7(std::min((size_t)arg_out.extent(7), (size_t)arg_in.extent(7))) {
    using Policy = Kokkos::RangePolicy<ExecSpace>;

    Kokkos::parallel_for("Kokkos::DynRankViewRemap", Policy(0, n0), *this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const size_t i0) const {
    for (size_t i1 = 0; i1 < n1; ++i1) {
      for (size_t i2 = 0; i2 < n2; ++i2) {
        for (size_t i3 = 0; i3 < n3; ++i3) {
          for (size_t i4 = 0; i4 < n4; ++i4) {
            for (size_t i5 = 0; i5 < n5; ++i5) {
              for (size_t i6 = 0; i6 < n6; ++i6) {
                output.access(i0, i1, i2, i3, i4, i5, i6) =
                    input.access(i0, i1, i2, i3, i4, i5, i6);
              }
            }
          }
        }
      }
    }
  }
};

} /* namespace Impl */
} /* namespace Kokkos */

namespace Kokkos {

namespace Impl {

/* \brief Returns a View of the requested rank, aliasing the
   underlying memory, to facilitate implementation of deep_copy() and
   other routines that are defined on View */
template <unsigned N, typename T, typename... Args>
KOKKOS_FUNCTION auto as_view_of_rank_n(
    DynRankView<T, Args...> v,
    typename std::enable_if<std::is_same<
        typename ViewTraits<T, Args...>::specialize, void>::value>::type* =
        nullptr) {
  if (v.rank() != N) {
    KOKKOS_IF_ON_HOST(
        const std::string message =
            "Converting DynRankView of rank " + std::to_string(v.rank()) +
            " to a View of mis-matched rank " + std::to_string(N) + "!";
        Kokkos::abort(message.c_str());)
    KOKKOS_IF_ON_DEVICE(
        Kokkos::abort("Converting DynRankView to a View of mis-matched rank!");)
  }

  auto layout = v.impl_map().layout();

  if constexpr (std::is_same_v<decltype(layout), Kokkos::LayoutLeft> ||
                std::is_same_v<decltype(layout), Kokkos::LayoutRight> ||
                std::is_same_v<decltype(layout), Kokkos::LayoutStride>) {
    for (int i = N; i < 7; ++i)
      layout.dimension[i] = KOKKOS_IMPL_CTOR_DEFAULT_ARG;
  }

  return View<typename RankDataType<T, N>::type, Args...>(v.data(), layout);
}

template <typename Function, typename... Args>
void apply_to_view_of_static_rank(Function&& f, DynRankView<Args...> a) {
  switch (rank(a)) {
    case 0: f(as_view_of_rank_n<0>(a)); break;
    case 1: f(as_view_of_rank_n<1>(a)); break;
    case 2: f(as_view_of_rank_n<2>(a)); break;
    case 3: f(as_view_of_rank_n<3>(a)); break;
    case 4: f(as_view_of_rank_n<4>(a)); break;
    case 5: f(as_view_of_rank_n<5>(a)); break;
    case 6: f(as_view_of_rank_n<6>(a)); break;
    case 7: f(as_view_of_rank_n<7>(a)); break;
    default:
      KOKKOS_IF_ON_HOST(
          Kokkos::abort(
              std::string(
                  "Trying to apply a function to a view of unexpected rank " +
                  std::to_string(rank(a)))
                  .c_str());)
      KOKKOS_IF_ON_DEVICE(
          Kokkos::abort(
              "Trying to apply a function to a view of unexpected rank");)
  }
}

}  // namespace Impl

template <typename D, class... P>
KOKKOS_INLINE_FUNCTION constexpr auto DynRankView<D, P...>::layout() const ->
    typename traits::array_layout {
  switch (rank()) {
    case 0: return Impl::as_view_of_rank_n<0>(*this).layout();
    case 1: return Impl::as_view_of_rank_n<1>(*this).layout();
    case 2: return Impl::as_view_of_rank_n<2>(*this).layout();
    case 3: return Impl::as_view_of_rank_n<3>(*this).layout();
    case 4: return Impl::as_view_of_rank_n<4>(*this).layout();
    case 5: return Impl::as_view_of_rank_n<5>(*this).layout();
    case 6: return Impl::as_view_of_rank_n<6>(*this).layout();
    case 7: return Impl::as_view_of_rank_n<7>(*this).layout();
    default:
      KOKKOS_IF_ON_HOST(
          Kokkos::abort(
              std::string(
                  "Calling DynRankView::layout on DRV of unexpected rank " +
                  std::to_string(rank()))
                  .c_str());)
      KOKKOS_IF_ON_DEVICE(
          Kokkos::abort(
              "Calling DynRankView::layout on DRV of unexpected rank");)
  }
  // control flow should never reach here
  return m_map.layout();
}

/** \brief  Deep copy a value from Host memory into a view.  */
template <class ExecSpace, class DT, class... DP>
inline void deep_copy(
    const ExecSpace& e, const DynRankView<DT, DP...>& dst,
    typename ViewTraits<DT, DP...>::const_value_type& value,
    std::enable_if_t<std::is_same<typename ViewTraits<DT, DP...>::specialize,
                                  void>::value>* = nullptr) {
  static_assert(
      std::is_same<typename ViewTraits<DT, DP...>::non_const_value_type,
                   typename ViewTraits<DT, DP...>::value_type>::value,
      "deep_copy requires non-const type");

  Impl::apply_to_view_of_static_rank(
      [=](auto view) { deep_copy(e, view, value); }, dst);
}

template <class DT, class... DP>
inline void deep_copy(
    const DynRankView<DT, DP...>& dst,
    typename ViewTraits<DT, DP...>::const_value_type& value,
    std::enable_if_t<std::is_same<typename ViewTraits<DT, DP...>::specialize,
                                  void>::value>* = nullptr) {
  Impl::apply_to_view_of_static_rank([=](auto view) { deep_copy(view, value); },
                                     dst);
}

/** \brief  Deep copy into a value in Host memory from a view.  */
template <class ExecSpace, class ST, class... SP>
inline void deep_copy(
    const ExecSpace& e,
    typename ViewTraits<ST, SP...>::non_const_value_type& dst,
    const DynRankView<ST, SP...>& src,
    std::enable_if_t<std::is_same<typename ViewTraits<ST, SP...>::specialize,
                                  void>::value>* = 0) {
  deep_copy(e, dst, Impl::as_view_of_rank_n<0>(src));
}

template <class ST, class... SP>
inline void deep_copy(
    typename ViewTraits<ST, SP...>::non_const_value_type& dst,
    const DynRankView<ST, SP...>& src,
    std::enable_if_t<std::is_same<typename ViewTraits<ST, SP...>::specialize,
                                  void>::value>* = 0) {
  deep_copy(dst, Impl::as_view_of_rank_n<0>(src));
}

//----------------------------------------------------------------------------
/** \brief  A deep copy between views of the default specialization, compatible
 * type, same rank, same contiguous layout.
 *
 * A rank mismatch will error out in the attempt to convert to a View
 */
template <class ExecSpace, class DstType, class SrcType>
inline void deep_copy(
    const ExecSpace& exec_space, const DstType& dst, const SrcType& src,
    std::enable_if_t<
        (std::is_void<typename DstType::traits::specialize>::value &&
         std::is_void<typename SrcType::traits::specialize>::value &&
         (Kokkos::is_dyn_rank_view<DstType>::value ||
          Kokkos::is_dyn_rank_view<SrcType>::value))>* = nullptr) {
  static_assert(
      std::is_same<typename DstType::traits::value_type,
                   typename DstType::traits::non_const_value_type>::value,
      "deep_copy requires non-const destination type");

  switch (rank(dst)) {
    case 0:
      deep_copy(exec_space, Impl::as_view_of_rank_n<0>(dst),
                Impl::as_view_of_rank_n<0>(src));
      break;
    case 1:
      deep_copy(exec_space, Impl::as_view_of_rank_n<1>(dst),
                Impl::as_view_of_rank_n<1>(src));
      break;
    case 2:
      deep_copy(exec_space, Impl::as_view_of_rank_n<2>(dst),
                Impl::as_view_of_rank_n<2>(src));
      break;
    case 3:
      deep_copy(exec_space, Impl::as_view_of_rank_n<3>(dst),
                Impl::as_view_of_rank_n<3>(src));
      break;
    case 4:
      deep_copy(exec_space, Impl::as_view_of_rank_n<4>(dst),
                Impl::as_view_of_rank_n<4>(src));
      break;
    case 5:
      deep_copy(exec_space, Impl::as_view_of_rank_n<5>(dst),
                Impl::as_view_of_rank_n<5>(src));
      break;
    case 6:
      deep_copy(exec_space, Impl::as_view_of_rank_n<6>(dst),
                Impl::as_view_of_rank_n<6>(src));
      break;
    case 7:
      deep_copy(exec_space, Impl::as_view_of_rank_n<7>(dst),
                Impl::as_view_of_rank_n<7>(src));
      break;
    default:
      Kokkos::Impl::throw_runtime_exception(
          "Calling DynRankView deep_copy with a view of unexpected rank " +
          std::to_string(rank(dst)));
  }
}

template <class DstType, class SrcType>
inline void deep_copy(
    const DstType& dst, const SrcType& src,
    std::enable_if_t<
        (std::is_void<typename DstType::traits::specialize>::value &&
         std::is_void<typename SrcType::traits::specialize>::value &&
         (Kokkos::is_dyn_rank_view<DstType>::value ||
          Kokkos::is_dyn_rank_view<SrcType>::value))>* = nullptr) {
  static_assert(
      std::is_same<typename DstType::traits::value_type,
                   typename DstType::traits::non_const_value_type>::value,
      "deep_copy requires non-const destination type");

  switch (rank(dst)) {
    case 0:
      deep_copy(Impl::as_view_of_rank_n<0>(dst),
                Impl::as_view_of_rank_n<0>(src));
      break;
    case 1:
      deep_copy(Impl::as_view_of_rank_n<1>(dst),
                Impl::as_view_of_rank_n<1>(src));
      break;
    case 2:
      deep_copy(Impl::as_view_of_rank_n<2>(dst),
                Impl::as_view_of_rank_n<2>(src));
      break;
    case 3:
      deep_copy(Impl::as_view_of_rank_n<3>(dst),
                Impl::as_view_of_rank_n<3>(src));
      break;
    case 4:
      deep_copy(Impl::as_view_of_rank_n<4>(dst),
                Impl::as_view_of_rank_n<4>(src));
      break;
    case 5:
      deep_copy(Impl::as_view_of_rank_n<5>(dst),
                Impl::as_view_of_rank_n<5>(src));
      break;
    case 6:
      deep_copy(Impl::as_view_of_rank_n<6>(dst),
                Impl::as_view_of_rank_n<6>(src));
      break;
    case 7:
      deep_copy(Impl::as_view_of_rank_n<7>(dst),
                Impl::as_view_of_rank_n<7>(src));
      break;
    default:
      Kokkos::Impl::throw_runtime_exception(
          "Calling DynRankView deep_copy with a view of unexpected rank " +
          std::to_string(rank(dst)));
  }
}

}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

// Deduce Mirror Types
template <class Space, class T, class... P>
struct MirrorDRViewType {
  // The incoming view_type
  using src_view_type = typename Kokkos::DynRankView<T, P...>;
  // The memory space for the mirror view
  using memory_space = typename Space::memory_space;
  // Check whether it is the same memory space
  enum {
    is_same_memspace =
        std::is_same<memory_space, typename src_view_type::memory_space>::value
  };
  // The array_layout
  using array_layout = typename src_view_type::array_layout;
  // The data type (we probably want it non-const since otherwise we can't even
  // deep_copy to it.
  using data_type = typename src_view_type::non_const_data_type;
  // The destination view type if it is not the same memory space
  using dest_view_type = Kokkos::DynRankView<data_type, array_layout, Space>;
  // If it is the same memory_space return the existsing view_type
  // This will also keep the unmanaged trait if necessary
  using view_type =
      std::conditional_t<is_same_memspace, src_view_type, dest_view_type>;
};

template <class Space, class T, class... P>
struct MirrorDRVType {
  // The incoming view_type
  using src_view_type = typename Kokkos::DynRankView<T, P...>;
  // The memory space for the mirror view
  using memory_space = typename Space::memory_space;
  // Check whether it is the same memory space
  enum {
    is_same_memspace =
        std::is_same<memory_space, typename src_view_type::memory_space>::value
  };
  // The array_layout
  using array_layout = typename src_view_type::array_layout;
  // The data type (we probably want it non-const since otherwise we can't even
  // deep_copy to it.
  using data_type = typename src_view_type::non_const_data_type;
  // The destination view type if it is not the same memory space
  using view_type = Kokkos::DynRankView<data_type, array_layout, Space>;
};

}  // namespace Impl

namespace Impl {

// create a mirror
// private interface that accepts arbitrary view constructor args passed by a
// view_alloc
template <class T, class... P, class... ViewCtorArgs>
inline auto create_mirror(const DynRankView<T, P...>& src,
                          const Impl::ViewCtorProp<ViewCtorArgs...>& arg_prop) {
  check_view_ctor_args_create_mirror<ViewCtorArgs...>();

  auto prop_copy = Impl::with_properties_if_unset(
      arg_prop, std::string(src.label()).append("_mirror"));

  if constexpr (Impl::ViewCtorProp<ViewCtorArgs...>::has_memory_space) {
    using dst_type = typename Impl::MirrorDRVType<
        typename Impl::ViewCtorProp<ViewCtorArgs...>::memory_space, T,
        P...>::view_type;

    return dst_type(prop_copy,
                    Impl::reconstructLayout(src.layout(), src.rank()));
  } else {
    using src_type = DynRankView<T, P...>;
    using dst_type = typename src_type::HostMirror;

    return dst_type(prop_copy,
                    Impl::reconstructLayout(src.layout(), src.rank()));
  }
#if defined(KOKKOS_COMPILER_INTEL) ||                                 \
    (defined(KOKKOS_COMPILER_NVCC) && KOKKOS_COMPILER_NVCC >= 1130 && \
     !defined(KOKKOS_COMPILER_MSVC))
  __builtin_unreachable();
#endif
}

}  // namespace Impl

// public interface
template <class T, class... P,
          class Enable = std::enable_if_t<
              std::is_void_v<typename ViewTraits<T, P...>::specialize>>>
inline auto create_mirror(const DynRankView<T, P...>& src) {
  return Impl::create_mirror(src, Kokkos::view_alloc());
}

// public interface that accepts a without initializing flag
template <class T, class... P,
          class Enable = std::enable_if_t<
              std::is_void_v<typename ViewTraits<T, P...>::specialize>>>
inline auto create_mirror(Kokkos::Impl::WithoutInitializing_t wi,
                          const DynRankView<T, P...>& src) {
  return Impl::create_mirror(src, Kokkos::view_alloc(wi));
}

// public interface that accepts a space
template <class Space, class T, class... P,
          class Enable = std::enable_if_t<
              Kokkos::is_space<Space>::value &&
              std::is_void_v<typename ViewTraits<T, P...>::specialize>>>
auto create_mirror(const Space&, const Kokkos::DynRankView<T, P...>& src) {
  return Impl::create_mirror(
      src, Kokkos::view_alloc(typename Space::memory_space{}));
}

// public interface that accepts a space and a without initializing flag
template <class Space, class T, class... P,
          class Enable = std::enable_if_t<
              Kokkos::is_space<Space>::value &&
              std::is_void_v<typename ViewTraits<T, P...>::specialize>>>
auto create_mirror(Kokkos::Impl::WithoutInitializing_t wi, const Space&,
                   const Kokkos::DynRankView<T, P...>& src) {
  return Impl::create_mirror(
      src, Kokkos::view_alloc(wi, typename Space::memory_space{}));
}

// public interface that accepts arbitrary view constructor args passed by a
// view_alloc
template <class T, class... P, class... ViewCtorArgs,
          typename Enable = std::enable_if_t<
              std::is_void_v<typename ViewTraits<T, P...>::specialize>>>
inline auto create_mirror(const Impl::ViewCtorProp<ViewCtorArgs...>& arg_prop,
                          const DynRankView<T, P...>& src) {
  return Impl::create_mirror(src, arg_prop);
}

namespace Impl {

// create a mirror view
// private interface that accepts arbitrary view constructor args passed by a
// view_alloc
template <class T, class... P, class... ViewCtorArgs>
inline auto create_mirror_view(
    const DynRankView<T, P...>& src,
    [[maybe_unused]] const typename Impl::ViewCtorProp<ViewCtorArgs...>&
        arg_prop) {
  if constexpr (!Impl::ViewCtorProp<ViewCtorArgs...>::has_memory_space) {
    if constexpr (std::is_same<typename DynRankView<T, P...>::memory_space,
                               typename DynRankView<
                                   T, P...>::HostMirror::memory_space>::value &&
                  std::is_same<typename DynRankView<T, P...>::data_type,
                               typename DynRankView<
                                   T, P...>::HostMirror::data_type>::value) {
      return typename DynRankView<T, P...>::HostMirror(src);
    } else {
      return Kokkos::Impl::choose_create_mirror(src, arg_prop);
    }
  } else {
    if constexpr (Impl::MirrorDRViewType<typename Impl::ViewCtorProp<
                                             ViewCtorArgs...>::memory_space,
                                         T, P...>::is_same_memspace) {
      return typename Impl::MirrorDRViewType<
          typename Impl::ViewCtorProp<ViewCtorArgs...>::memory_space, T,
          P...>::view_type(src);
    } else {
      return Kokkos::Impl::choose_create_mirror(src, arg_prop);
    }
  }
#if defined(KOKKOS_COMPILER_INTEL) ||                                 \
    (defined(KOKKOS_COMPILER_NVCC) && KOKKOS_COMPILER_NVCC >= 1130 && \
     !defined(KOKKOS_COMPILER_MSVC))
  __builtin_unreachable();
#endif
}

}  // namespace Impl

// public interface
template <class T, class... P>
inline auto create_mirror_view(const Kokkos::DynRankView<T, P...>& src) {
  return Impl::create_mirror_view(src, Kokkos::view_alloc());
}

// public interface that accepts a without initializing flag
template <class T, class... P>
inline auto create_mirror_view(Kokkos::Impl::WithoutInitializing_t wi,
                               const DynRankView<T, P...>& src) {
  return Impl::create_mirror_view(src, Kokkos::view_alloc(wi));
}

// public interface that accepts a space
template <class Space, class T, class... P,
          class Enable = std::enable_if_t<Kokkos::is_space<Space>::value>>
inline auto create_mirror_view(const Space&,
                               const Kokkos::DynRankView<T, P...>& src) {
  return Impl::create_mirror_view(
      src, Kokkos::view_alloc(typename Space::memory_space()));
}

// public interface that accepts a space and a without initializing flag
template <class Space, class T, class... P,
          typename Enable = std::enable_if_t<Kokkos::is_space<Space>::value>>
inline auto create_mirror_view(Kokkos::Impl::WithoutInitializing_t wi,
                               const Space&,
                               const Kokkos::DynRankView<T, P...>& src) {
  return Impl::create_mirror_view(
      src, Kokkos::view_alloc(typename Space::memory_space{}, wi));
}

// public interface that accepts arbitrary view constructor args passed by a
// view_alloc
template <class T, class... P, class... ViewCtorArgs>
inline auto create_mirror_view(
    const typename Impl::ViewCtorProp<ViewCtorArgs...>& arg_prop,
    const Kokkos::DynRankView<T, P...>& src) {
  return Impl::create_mirror_view(src, arg_prop);
}

// create a mirror view and deep copy it
// public interface that accepts arbitrary view constructor args passed by a
// view_alloc
template <class... ViewCtorArgs, class T, class... P,
          class Enable = std::enable_if_t<
              std::is_void<typename ViewTraits<T, P...>::specialize>::value>>
auto create_mirror_view_and_copy(
    [[maybe_unused]] const Impl::ViewCtorProp<ViewCtorArgs...>& arg_prop,
    const Kokkos::DynRankView<T, P...>& src) {
  using alloc_prop_input = Impl::ViewCtorProp<ViewCtorArgs...>;

  Impl::check_view_ctor_args_create_mirror_view_and_copy<ViewCtorArgs...>();

  if constexpr (Impl::MirrorDRViewType<
                    typename Impl::ViewCtorProp<ViewCtorArgs...>::memory_space,
                    T, P...>::is_same_memspace) {
    // same behavior as deep_copy(src, src)
    if constexpr (!alloc_prop_input::has_execution_space)
      fence(
          "Kokkos::create_mirror_view_and_copy: fence before returning src "
          "view");
    return src;
  } else {
    using Space  = typename alloc_prop_input::memory_space;
    using Mirror = typename Impl::MirrorDRViewType<Space, T, P...>::view_type;

    auto arg_prop_copy = Impl::with_properties_if_unset(
        arg_prop, std::string{}, WithoutInitializing,
        typename Space::execution_space{});

    std::string& label = Impl::get_property<Impl::LabelTag>(arg_prop_copy);
    if (label.empty()) label = src.label();
    auto mirror = typename Mirror::non_const_type{
        arg_prop_copy, Impl::reconstructLayout(src.layout(), src.rank())};
    if constexpr (alloc_prop_input::has_execution_space) {
      deep_copy(Impl::get_property<Impl::ExecutionSpaceTag>(arg_prop_copy),
                mirror, src);
    } else
      deep_copy(mirror, src);
    return mirror;
  }
#if defined(KOKKOS_COMPILER_NVCC) && KOKKOS_COMPILER_NVCC >= 1130 && \
    !defined(KOKKOS_COMPILER_MSVC)
  __builtin_unreachable();
#endif
}

template <class Space, class T, class... P>
auto create_mirror_view_and_copy(const Space&,
                                 const Kokkos::DynRankView<T, P...>& src,
                                 std::string const& name = "") {
  return create_mirror_view_and_copy(
      Kokkos::view_alloc(typename Space::memory_space{}, name), src);
}

}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
/** \brief  Resize a view with copying old data to new data at the corresponding
 * indices. */
template <class... ViewCtorArgs, class T, class... P>
inline void impl_resize(const Impl::ViewCtorProp<ViewCtorArgs...>& arg_prop,
                        DynRankView<T, P...>& v, const size_t n0,
                        const size_t n1, const size_t n2, const size_t n3,
                        const size_t n4, const size_t n5, const size_t n6,
                        const size_t n7) {
  using drview_type      = DynRankView<T, P...>;
  using alloc_prop_input = Impl::ViewCtorProp<ViewCtorArgs...>;

  static_assert(Kokkos::ViewTraits<T, P...>::is_managed,
                "Can only resize managed views");
  static_assert(!alloc_prop_input::has_label,
                "The view constructor arguments passed to Kokkos::resize "
                "must not include a label!");
  static_assert(!alloc_prop_input::has_pointer,
                "The view constructor arguments passed to Kokkos::resize must "
                "not include a pointer!");
  static_assert(!alloc_prop_input::has_memory_space,
                "The view constructor arguments passed to Kokkos::resize must "
                "not include a memory space instance!");

  auto prop_copy = Impl::with_properties_if_unset(
      arg_prop, v.label(), typename drview_type::execution_space{});

  drview_type v_resized(prop_copy, n0, n1, n2, n3, n4, n5, n6, n7);

  if constexpr (alloc_prop_input::has_execution_space)
    Kokkos::Impl::DynRankViewRemap<drview_type, drview_type>(
        Impl::get_property<Impl::ExecutionSpaceTag>(prop_copy), v_resized, v);
  else {
    Kokkos::Impl::DynRankViewRemap<drview_type, drview_type>(v_resized, v);
    Kokkos::fence("Kokkos::resize(DynRankView)");
  }
  v = v_resized;
}

template <class T, class... P>
inline void resize(DynRankView<T, P...>& v,
                   const size_t n0 = KOKKOS_INVALID_INDEX,
                   const size_t n1 = KOKKOS_INVALID_INDEX,
                   const size_t n2 = KOKKOS_INVALID_INDEX,
                   const size_t n3 = KOKKOS_INVALID_INDEX,
                   const size_t n4 = KOKKOS_INVALID_INDEX,
                   const size_t n5 = KOKKOS_INVALID_INDEX,
                   const size_t n6 = KOKKOS_INVALID_INDEX,
                   const size_t n7 = KOKKOS_INVALID_INDEX) {
  impl_resize(Impl::ViewCtorProp<>{}, v, n0, n1, n2, n3, n4, n5, n6, n7);
}

template <class... ViewCtorArgs, class T, class... P>
void resize(const Impl::ViewCtorProp<ViewCtorArgs...>& arg_prop,
            DynRankView<T, P...>& v,
            const size_t n0 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
            const size_t n1 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
            const size_t n2 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
            const size_t n3 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
            const size_t n4 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
            const size_t n5 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
            const size_t n6 = KOKKOS_IMPL_CTOR_DEFAULT_ARG,
            const size_t n7 = KOKKOS_IMPL_CTOR_DEFAULT_ARG) {
  impl_resize(arg_prop, v, n0, n1, n2, n3, n4, n5, n6, n7);
}

template <class I, class T, class... P>
inline std::enable_if_t<Impl::is_view_ctor_property<I>::value> resize(
    const I& arg_prop, DynRankView<T, P...>& v,
    const size_t n0 = KOKKOS_INVALID_INDEX,
    const size_t n1 = KOKKOS_INVALID_INDEX,
    const size_t n2 = KOKKOS_INVALID_INDEX,
    const size_t n3 = KOKKOS_INVALID_INDEX,
    const size_t n4 = KOKKOS_INVALID_INDEX,
    const size_t n5 = KOKKOS_INVALID_INDEX,
    const size_t n6 = KOKKOS_INVALID_INDEX,
    const size_t n7 = KOKKOS_INVALID_INDEX) {
  impl_resize(Kokkos::view_alloc(arg_prop), v, n0, n1, n2, n3, n4, n5, n6, n7);
}

/** \brief  Resize a view with copying old data to new data at the corresponding
 * indices. */
template <class... ViewCtorArgs, class T, class... P>
inline void impl_realloc(DynRankView<T, P...>& v, const size_t n0,
                         const size_t n1, const size_t n2, const size_t n3,
                         const size_t n4, const size_t n5, const size_t n6,
                         const size_t n7,
                         const Impl::ViewCtorProp<ViewCtorArgs...>& arg_prop) {
  using drview_type      = DynRankView<T, P...>;
  using alloc_prop_input = Impl::ViewCtorProp<ViewCtorArgs...>;

  static_assert(Kokkos::ViewTraits<T, P...>::is_managed,
                "Can only realloc managed views");
  static_assert(!alloc_prop_input::has_label,
                "The view constructor arguments passed to Kokkos::realloc must "
                "not include a label!");
  static_assert(!alloc_prop_input::has_pointer,
                "The view constructor arguments passed to Kokkos::realloc must "
                "not include a pointer!");
  static_assert(!alloc_prop_input::has_memory_space,
                "The view constructor arguments passed to Kokkos::realloc must "
                "not include a memory space instance!");

  auto arg_prop_copy = Impl::with_properties_if_unset(arg_prop, v.label());

  v = drview_type();  // Deallocate first, if the only view to allocation
  v = drview_type(arg_prop_copy, n0, n1, n2, n3, n4, n5, n6, n7);
}

template <class T, class... P, class... ViewCtorArgs>
inline void realloc(const Impl::ViewCtorProp<ViewCtorArgs...>& arg_prop,
                    DynRankView<T, P...>& v,
                    const size_t n0 = KOKKOS_INVALID_INDEX,
                    const size_t n1 = KOKKOS_INVALID_INDEX,
                    const size_t n2 = KOKKOS_INVALID_INDEX,
                    const size_t n3 = KOKKOS_INVALID_INDEX,
                    const size_t n4 = KOKKOS_INVALID_INDEX,
                    const size_t n5 = KOKKOS_INVALID_INDEX,
                    const size_t n6 = KOKKOS_INVALID_INDEX,
                    const size_t n7 = KOKKOS_INVALID_INDEX) {
  impl_realloc(v, n0, n1, n2, n3, n4, n5, n6, n7, arg_prop);
}

template <class T, class... P>
inline void realloc(DynRankView<T, P...>& v,
                    const size_t n0 = KOKKOS_INVALID_INDEX,
                    const size_t n1 = KOKKOS_INVALID_INDEX,
                    const size_t n2 = KOKKOS_INVALID_INDEX,
                    const size_t n3 = KOKKOS_INVALID_INDEX,
                    const size_t n4 = KOKKOS_INVALID_INDEX,
                    const size_t n5 = KOKKOS_INVALID_INDEX,
                    const size_t n6 = KOKKOS_INVALID_INDEX,
                    const size_t n7 = KOKKOS_INVALID_INDEX) {
  impl_realloc(v, n0, n1, n2, n3, n4, n5, n6, n7, Impl::ViewCtorProp<>{});
}

template <class I, class T, class... P>
inline std::enable_if_t<Impl::is_view_ctor_property<I>::value> realloc(
    const I& arg_prop, DynRankView<T, P...>& v,
    const size_t n0 = KOKKOS_INVALID_INDEX,
    const size_t n1 = KOKKOS_INVALID_INDEX,
    const size_t n2 = KOKKOS_INVALID_INDEX,
    const size_t n3 = KOKKOS_INVALID_INDEX,
    const size_t n4 = KOKKOS_INVALID_INDEX,
    const size_t n5 = KOKKOS_INVALID_INDEX,
    const size_t n6 = KOKKOS_INVALID_INDEX,
    const size_t n7 = KOKKOS_INVALID_INDEX) {
  impl_realloc(v, n0, n1, n2, n3, n4, n5, n6, n7, Kokkos::view_alloc(arg_prop));
}

}  // namespace Kokkos

#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_DYNRANKVIEW
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_DYNRANKVIEW
#endif
#endif
