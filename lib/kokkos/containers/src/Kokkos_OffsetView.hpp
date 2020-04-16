/*
 * Kokkos_OffsetView.hpp
 *
 *  Created on: Apr 23, 2018
 *      Author: swbova
 */

#ifndef KOKKOS_OFFSETVIEW_HPP_
#define KOKKOS_OFFSETVIEW_HPP_

#include <Kokkos_Core.hpp>

#include <Kokkos_View.hpp>

namespace Kokkos {

namespace Experimental {
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template <class DataType, class... Properties>
class OffsetView;

template <class>
struct is_offset_view : public std::false_type {};

template <class D, class... P>
struct is_offset_view<OffsetView<D, P...> > : public std::true_type {};

template <class D, class... P>
struct is_offset_view<const OffsetView<D, P...> > : public std::true_type {};

#define KOKKOS_INVALID_OFFSET int64_t(0x7FFFFFFFFFFFFFFFLL)
#define KOKKOS_INVALID_INDEX_RANGE \
  { KOKKOS_INVALID_OFFSET, KOKKOS_INVALID_OFFSET }

template <typename iType,
          typename std::enable_if<std::is_integral<iType>::value &&
                                      std::is_signed<iType>::value,
                                  iType>::type = 0>
using IndexRange = Kokkos::Array<iType, 2>;

using index_list_type = std::initializer_list<int64_t>;

//  template <typename iType,
//    typename std::enable_if< std::is_integral<iType>::value &&
//      std::is_signed<iType>::value, iType >::type = 0> using min_index_type =
//      std::initializer_list<iType>;

namespace Impl {

template <class ViewType>
struct GetOffsetViewTypeFromViewType {
  typedef OffsetView<
      typename ViewType::data_type, typename ViewType::array_layout,
      typename ViewType::device_type, typename ViewType::memory_traits>
      type;
};

template <unsigned, class MapType, class BeginsType>
KOKKOS_INLINE_FUNCTION bool offsetview_verify_operator_bounds(
    const MapType&, const BeginsType&) {
  return true;
}

template <unsigned R, class MapType, class BeginsType, class iType,
          class... Args>
KOKKOS_INLINE_FUNCTION bool offsetview_verify_operator_bounds(
    const MapType& map, const BeginsType& begins, const iType& i,
    Args... args) {
  const bool legalIndex =
      (int64_t(i) >= begins[R]) &&
      (int64_t(i) <= int64_t(begins[R] + map.extent(R) - 1));
  return legalIndex &&
         offsetview_verify_operator_bounds<R + 1>(map, begins, args...);
}
template <unsigned, class MapType, class BeginsType>
inline void offsetview_error_operator_bounds(char*, int, const MapType&,
                                             const BeginsType&) {}

template <unsigned R, class MapType, class BeginsType, class iType,
          class... Args>
inline void offsetview_error_operator_bounds(char* buf, int len,
                                             const MapType& map,
                                             const BeginsType begins,
                                             const iType& i, Args... args) {
  const int64_t b = begins[R];
  const int64_t e = b + map.extent(R) - 1;
  const int n =
      snprintf(buf, len, " %ld <= %ld <= %ld %c", static_cast<unsigned long>(b),
               static_cast<unsigned long>(i), static_cast<unsigned long>(e),
               (sizeof...(Args) ? ',' : ')'));
  offsetview_error_operator_bounds<R + 1>(buf + n, len - n, map, begins,
                                          args...);
}

template <class MemorySpace, class MapType, class BeginsType, class... Args>
KOKKOS_INLINE_FUNCTION void offsetview_verify_operator_bounds(
    Kokkos::Impl::SharedAllocationTracker const& tracker, const MapType& map,
    const BeginsType& begins, Args... args) {
  if (!offsetview_verify_operator_bounds<0>(map, begins, args...)) {
#if defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST)
    enum { LEN = 1024 };
    char buffer[LEN];
    const std::string label = tracker.template get_label<MemorySpace>();
    int n =
        snprintf(buffer, LEN, "OffsetView bounds error of view labeled %s (",
                 label.c_str());
    offsetview_error_operator_bounds<0>(buffer + n, LEN - n, map, begins,
                                        args...);
    Kokkos::Impl::throw_runtime_exception(std::string(buffer));
#else
    /* Check #1: is there a SharedAllocationRecord?
      (we won't use it, but if its not there then there isn't
       a corresponding SharedAllocationHeader containing a label).
      This check should cover the case of Views that don't
      have the Unmanaged trait but were initialized by pointer. */
    if (tracker.has_record()) {
      Kokkos::Impl::operator_bounds_error_on_device<MapType>(
          map, Kokkos::Impl::has_printable_label_typedef<MapType>());
    } else {
      Kokkos::abort("OffsetView bounds error");
    }
#endif
  }
}

#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
KOKKOS_INLINE_FUNCTION
void runtime_check_rank_host(const size_t rank_dynamic, const size_t rank,
                             const index_list_type minIndices,
                             const std::string& label) {
  bool isBad = false;
  std::string message =
      "Kokkos::Experimental::OffsetView ERROR: for OffsetView labeled '" +
      label + "':";
  if (rank_dynamic != rank) {
    message +=
        "The full rank must be the same as the dynamic rank. full rank = ";
    message += std::to_string(rank) +
               " dynamic rank = " + std::to_string(rank_dynamic) + "\n";
    isBad = true;
  }

  size_t numOffsets = 0;
  for (size_t i = 0; i < minIndices.size(); ++i) {
    if (minIndices.begin()[i] != KOKKOS_INVALID_OFFSET) numOffsets++;
  }
  if (numOffsets != rank_dynamic) {
    message += "The number of offsets provided ( " +
               std::to_string(numOffsets) +
               " ) must equal the dynamic rank ( " +
               std::to_string(rank_dynamic) + " ).";
    isBad = true;
  }

  if (isBad) Kokkos::abort(message.c_str());
}
#endif

KOKKOS_INLINE_FUNCTION
void runtime_check_rank_device(const size_t rank_dynamic, const size_t rank,
                               const index_list_type minIndices) {
  if (rank_dynamic != rank) {
    Kokkos::abort(
        "The full rank of an OffsetView must be the same as the dynamic rank.");
  }
  size_t numOffsets = 0;
  for (size_t i = 0; i < minIndices.size(); ++i) {
    if (minIndices.begin()[i] != KOKKOS_INVALID_OFFSET) numOffsets++;
  }
  if (numOffsets != rank) {
    Kokkos::abort(
        "The number of offsets provided to an OffsetView constructor must "
        "equal the dynamic rank.");
  }
}
}  // namespace Impl

template <class DataType, class... Properties>
class OffsetView : public ViewTraits<DataType, Properties...> {
 public:
  typedef ViewTraits<DataType, Properties...> traits;

 private:
  template <class, class...>
  friend class OffsetView;
  template <class, class...>
  friend class View;  // FIXME delete this line
  template <class, class...>
  friend class Kokkos::Impl::ViewMapping;

  typedef Kokkos::Impl::ViewMapping<traits, void> map_type;
  typedef Kokkos::Impl::SharedAllocationTracker track_type;

 public:
  enum { Rank = map_type::Rank };
  typedef Kokkos::Array<int64_t, Rank> begins_type;

  template <
      typename iType,
      typename std::enable_if<std::is_integral<iType>::value, iType>::type = 0>
  KOKKOS_INLINE_FUNCTION int64_t begin(const iType local_dimension) const {
    return local_dimension < Rank ? m_begins[local_dimension]
                                  : KOKKOS_INVALID_OFFSET;
  }

  KOKKOS_INLINE_FUNCTION
  begins_type begins() const { return m_begins; }

  template <
      typename iType,
      typename std::enable_if<std::is_integral<iType>::value, iType>::type = 0>
  KOKKOS_INLINE_FUNCTION int64_t end(const iType local_dimension) const {
    return begin(local_dimension) + m_map.extent(local_dimension);
  }

 private:
  track_type m_track;
  map_type m_map;
  begins_type m_begins;

 public:
  //----------------------------------------
  /** \brief  Compatible view of array of scalar types */
  typedef OffsetView<
      typename traits::scalar_array_type, typename traits::array_layout,
      typename traits::device_type, typename traits::memory_traits>
      array_type;

  /** \brief  Compatible view of const data type */
  typedef OffsetView<
      typename traits::const_data_type, typename traits::array_layout,
      typename traits::device_type, typename traits::memory_traits>
      const_type;

  /** \brief  Compatible view of non-const data type */
  typedef OffsetView<
      typename traits::non_const_data_type, typename traits::array_layout,
      typename traits::device_type, typename traits::memory_traits>
      non_const_type;

  /** \brief  Compatible HostMirror view */
  typedef OffsetView<typename traits::non_const_data_type,
                     typename traits::array_layout,
                     typename traits::host_mirror_space>
      HostMirror;

  //----------------------------------------
  // Domain rank and extents

  /** \brief rank() to be implemented
   */
  // KOKKOS_INLINE_FUNCTION
  // static
  // constexpr unsigned rank() { return map_type::Rank; }

  template <typename iType>
  KOKKOS_INLINE_FUNCTION constexpr
      typename std::enable_if<std::is_integral<iType>::value, size_t>::type
      extent(const iType& r) const {
    return m_map.extent(r);
  }

  template <typename iType>
  KOKKOS_INLINE_FUNCTION constexpr
      typename std::enable_if<std::is_integral<iType>::value, int>::type
      extent_int(const iType& r) const {
    return static_cast<int>(m_map.extent(r));
  }

  KOKKOS_INLINE_FUNCTION constexpr typename traits::array_layout layout()
      const {
    return m_map.layout();
  }

  KOKKOS_INLINE_FUNCTION constexpr size_t size() const {
    return m_map.dimension_0() * m_map.dimension_1() * m_map.dimension_2() *
           m_map.dimension_3() * m_map.dimension_4() * m_map.dimension_5() *
           m_map.dimension_6() * m_map.dimension_7();
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
  KOKKOS_INLINE_FUNCTION constexpr
      typename std::enable_if<std::is_integral<iType>::value, size_t>::type
      stride(iType r) const {
    return (
        r == 0
            ? m_map.stride_0()
            : (r == 1
                   ? m_map.stride_1()
                   : (r == 2
                          ? m_map.stride_2()
                          : (r == 3
                                 ? m_map.stride_3()
                                 : (r == 4
                                        ? m_map.stride_4()
                                        : (r == 5
                                               ? m_map.stride_5()
                                               : (r == 6
                                                      ? m_map.stride_6()
                                                      : m_map.stride_7())))))));
  }

  template <typename iType>
  KOKKOS_INLINE_FUNCTION void stride(iType* const s) const {
    m_map.stride(s);
  }

  //----------------------------------------
  // Range span is the span which contains all members.

  typedef typename map_type::reference_type reference_type;
  typedef typename map_type::pointer_type pointer_type;

  enum {
    reference_type_is_lvalue_reference =
        std::is_lvalue_reference<reference_type>::value
  };

  KOKKOS_INLINE_FUNCTION constexpr size_t span() const { return m_map.span(); }
  KOKKOS_INLINE_FUNCTION bool span_is_contiguous() const {
    return m_map.span_is_contiguous();
  }
  KOKKOS_INLINE_FUNCTION constexpr pointer_type data() const {
    return m_map.data();
  }

  //----------------------------------------
  // Allow specializations to query their specialized map

  KOKKOS_INLINE_FUNCTION
  const Kokkos::Impl::ViewMapping<traits, void>& implementation_map() const {
    return m_map;
  }

  //----------------------------------------

 private:
  static constexpr bool is_layout_left =
      std::is_same<typename traits::array_layout, Kokkos::LayoutLeft>::value;

  static constexpr bool is_layout_right =
      std::is_same<typename traits::array_layout, Kokkos::LayoutRight>::value;

  static constexpr bool is_layout_stride =
      std::is_same<typename traits::array_layout, Kokkos::LayoutStride>::value;

  static constexpr bool is_default_map =
      std::is_same<typename traits::specialize, void>::value &&
      (is_layout_left || is_layout_right || is_layout_stride);

  template <class Space, bool = Kokkos::Impl::MemorySpaceAccess<
                             Space, typename traits::memory_space>::accessible>
  struct verify_space {
    KOKKOS_FORCEINLINE_FUNCTION static void check() {}
  };

  template <class Space>
  struct verify_space<Space, false> {
    KOKKOS_FORCEINLINE_FUNCTION static void check() {
      Kokkos::abort(
          "Kokkos::View ERROR: attempt to access inaccessible memory space");
    };
  };

#if defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK)

#define KOKKOS_IMPL_OFFSETVIEW_OPERATOR_VERIFY(ARG)              \
  OffsetView::template verify_space<                             \
      Kokkos::Impl::ActiveExecutionMemorySpace>::check();        \
  Kokkos::Experimental::Impl::offsetview_verify_operator_bounds< \
      typename traits::memory_space>                             \
      ARG;

#else

#define KOKKOS_IMPL_OFFSETVIEW_OPERATOR_VERIFY(ARG) \
  OffsetView::template verify_space<                \
      Kokkos::Impl::ActiveExecutionMemorySpace>::check();

#endif
 public:
  //------------------------------
  // Rank 0 operator()

  KOKKOS_FORCEINLINE_FUNCTION
  reference_type operator()() const { return m_map.reference(); }
  //------------------------------
  // Rank 1 operator()

  template <typename I0>
  KOKKOS_FORCEINLINE_FUNCTION
      typename std::enable_if<(Kokkos::Impl::are_integral<I0>::value &&
                               (1 == Rank) && !is_default_map),
                              reference_type>::type
      operator()(const I0& i0) const {
    KOKKOS_IMPL_OFFSETVIEW_OPERATOR_VERIFY((m_track, m_map, m_begins, i0))
    const size_t j0 = i0 - m_begins[0];
    return m_map.reference(j0);
  }

  template <typename I0>
  KOKKOS_FORCEINLINE_FUNCTION
      typename std::enable_if<(Kokkos::Impl::are_integral<I0>::value &&
                               (1 == Rank) && is_default_map &&
                               !is_layout_stride),
                              reference_type>::type
      operator()(const I0& i0) const {
    KOKKOS_IMPL_OFFSETVIEW_OPERATOR_VERIFY((m_track, m_map, m_begins, i0))
    const size_t j0 = i0 - m_begins[0];
    return m_map.m_impl_handle[j0];
  }

  template <typename I0>
  KOKKOS_FORCEINLINE_FUNCTION
      typename std::enable_if<(Kokkos::Impl::are_integral<I0>::value &&
                               (1 == Rank) && is_default_map &&
                               is_layout_stride),
                              reference_type>::type
      operator()(const I0& i0) const {
    KOKKOS_IMPL_OFFSETVIEW_OPERATOR_VERIFY((m_track, m_map, m_begins, i0))
    const size_t j0 = i0 - m_begins[0];
    return m_map.m_impl_handle[m_map.m_impl_offset.m_stride.S0 * j0];
  }
  //------------------------------
  // Rank 1 operator[]

  template <typename I0>
  KOKKOS_FORCEINLINE_FUNCTION
      typename std::enable_if<(Kokkos::Impl::are_integral<I0>::value &&
                               (1 == Rank) && !is_default_map),
                              reference_type>::type
      operator[](const I0& i0) const {
    KOKKOS_IMPL_OFFSETVIEW_OPERATOR_VERIFY((m_track, m_map, m_begins, i0))
    const size_t j0 = i0 - m_begins[0];
    return m_map.reference(j0);
  }

  template <typename I0>
  KOKKOS_FORCEINLINE_FUNCTION
      typename std::enable_if<(Kokkos::Impl::are_integral<I0>::value &&
                               (1 == Rank) && is_default_map &&
                               !is_layout_stride),
                              reference_type>::type
      operator[](const I0& i0) const {
    KOKKOS_IMPL_OFFSETVIEW_OPERATOR_VERIFY((m_track, m_map, m_begins, i0))
    const size_t j0 = i0 - m_begins[0];
    return m_map.m_impl_handle[j0];
  }

  template <typename I0>
  KOKKOS_FORCEINLINE_FUNCTION
      typename std::enable_if<(Kokkos::Impl::are_integral<I0>::value &&
                               (1 == Rank) && is_default_map &&
                               is_layout_stride),
                              reference_type>::type
      operator[](const I0& i0) const {
    KOKKOS_IMPL_OFFSETVIEW_OPERATOR_VERIFY((m_track, m_map, m_begins, i0))
    const size_t j0 = i0 - m_begins[0];
    return m_map.m_impl_handle[m_map.m_impl_offset.m_stride.S0 * j0];
  }

  //------------------------------
  // Rank 2

  template <typename I0, typename I1>
  KOKKOS_FORCEINLINE_FUNCTION
      typename std::enable_if<(Kokkos::Impl::are_integral<I0, I1>::value &&
                               (2 == Rank) && !is_default_map),
                              reference_type>::type
      operator()(const I0& i0, const I1& i1) const {
    KOKKOS_IMPL_OFFSETVIEW_OPERATOR_VERIFY((m_track, m_map, m_begins, i0, i1))
    const size_t j0 = i0 - m_begins[0];
    const size_t j1 = i1 - m_begins[1];
    return m_map.reference(j0, j1);
  }

  template <typename I0, typename I1>
  KOKKOS_FORCEINLINE_FUNCTION
      typename std::enable_if<(Kokkos::Impl::are_integral<I0, I1>::value &&
                               (2 == Rank) && is_default_map &&
                               is_layout_left && (traits::rank_dynamic == 0)),
                              reference_type>::type
      operator()(const I0& i0, const I1& i1) const {
    KOKKOS_IMPL_OFFSETVIEW_OPERATOR_VERIFY((m_track, m_map, m_begins, i0, i1))
    const size_t j0 = i0 - m_begins[0];
    const size_t j1 = i1 - m_begins[1];
    return m_map.m_impl_handle[j0 + m_map.m_impl_offset.m_dim.N0 * j1];
  }

  template <typename I0, typename I1>
  KOKKOS_FORCEINLINE_FUNCTION
      typename std::enable_if<(Kokkos::Impl::are_integral<I0, I1>::value &&
                               (2 == Rank) && is_default_map &&
                               is_layout_left && (traits::rank_dynamic != 0)),
                              reference_type>::type
      operator()(const I0& i0, const I1& i1) const {
    KOKKOS_IMPL_OFFSETVIEW_OPERATOR_VERIFY((m_track, m_map, m_begins, i0, i1))
    const size_t j0 = i0 - m_begins[0];
    const size_t j1 = i1 - m_begins[1];
    return m_map.m_impl_handle[j0 + m_map.m_impl_offset.m_stride * j1];
  }

  template <typename I0, typename I1>
  KOKKOS_FORCEINLINE_FUNCTION
      typename std::enable_if<(Kokkos::Impl::are_integral<I0, I1>::value &&
                               (2 == Rank) && is_default_map &&
                               is_layout_right && (traits::rank_dynamic == 0)),
                              reference_type>::type
      operator()(const I0& i0, const I1& i1) const {
    KOKKOS_IMPL_OFFSETVIEW_OPERATOR_VERIFY((m_track, m_map, m_begins, i0, i1))
    const size_t j0 = i0 - m_begins[0];
    const size_t j1 = i1 - m_begins[1];
    return m_map.m_impl_handle[j1 + m_map.m_impl_offset.m_dim.N1 * j0];
  }

  template <typename I0, typename I1>
  KOKKOS_FORCEINLINE_FUNCTION
      typename std::enable_if<(Kokkos::Impl::are_integral<I0, I1>::value &&
                               (2 == Rank) && is_default_map &&
                               is_layout_right && (traits::rank_dynamic != 0)),
                              reference_type>::type
      operator()(const I0& i0, const I1& i1) const {
    KOKKOS_IMPL_OFFSETVIEW_OPERATOR_VERIFY((m_track, m_map, m_begins, i0, i1))
    const size_t j0 = i0 - m_begins[0];
    const size_t j1 = i1 - m_begins[1];
    return m_map.m_impl_handle[j1 + m_map.m_impl_offset.m_stride * j0];
  }

  template <typename I0, typename I1>
  KOKKOS_FORCEINLINE_FUNCTION
      typename std::enable_if<(Kokkos::Impl::are_integral<I0, I1>::value &&
                               (2 == Rank) && is_default_map &&
                               is_layout_stride),
                              reference_type>::type
      operator()(const I0& i0, const I1& i1) const {
    KOKKOS_IMPL_OFFSETVIEW_OPERATOR_VERIFY((m_track, m_map, m_begins, i0, i1))
    const size_t j0 = i0 - m_begins[0];
    const size_t j1 = i1 - m_begins[1];
    return m_map.m_impl_handle[j0 * m_map.m_impl_offset.m_stride.S0 +
                               j1 * m_map.m_impl_offset.m_stride.S1];
  }

  //------------------------------
  // Rank 3

  template <typename I0, typename I1, typename I2>
  KOKKOS_FORCEINLINE_FUNCTION
      typename std::enable_if<(Kokkos::Impl::are_integral<I0, I1, I2>::value &&
                               (3 == Rank) && is_default_map),
                              reference_type>::type
      operator()(const I0& i0, const I1& i1, const I2& i2) const {
    KOKKOS_IMPL_OFFSETVIEW_OPERATOR_VERIFY(
        (m_track, m_map, m_begins, i0, i1, i2))
    const size_t j0 = i0 - m_begins[0];
    const size_t j1 = i1 - m_begins[1];
    const size_t j2 = i2 - m_begins[2];
    return m_map.m_impl_handle[m_map.m_impl_offset(j0, j1, j2)];
  }

  template <typename I0, typename I1, typename I2>
  KOKKOS_FORCEINLINE_FUNCTION
      typename std::enable_if<(Kokkos::Impl::are_integral<I0, I1, I2>::value &&
                               (3 == Rank) && !is_default_map),
                              reference_type>::type
      operator()(const I0& i0, const I1& i1, const I2& i2) const {
    KOKKOS_IMPL_OFFSETVIEW_OPERATOR_VERIFY(
        (m_track, m_map, m_begins, i0, i1, i2))
    const size_t j0 = i0 - m_begins[0];
    const size_t j1 = i1 - m_begins[1];
    const size_t j2 = i2 - m_begins[2];
    return m_map.reference(j0, j1, j2);
  }

  //------------------------------
  // Rank 4

  template <typename I0, typename I1, typename I2, typename I3>
  KOKKOS_FORCEINLINE_FUNCTION typename std::enable_if<
      (Kokkos::Impl::are_integral<I0, I1, I2, I3>::value && (4 == Rank) &&
       is_default_map),
      reference_type>::type
  operator()(const I0& i0, const I1& i1, const I2& i2, const I3& i3) const {
    KOKKOS_IMPL_OFFSETVIEW_OPERATOR_VERIFY(
        (m_track, m_map, m_begins, i0, i1, i2, i3))
    const size_t j0 = i0 - m_begins[0];
    const size_t j1 = i1 - m_begins[1];
    const size_t j2 = i2 - m_begins[2];
    const size_t j3 = i3 - m_begins[3];
    return m_map.m_impl_handle[m_map.m_impl_offset(j0, j1, j2, j3)];
  }

  template <typename I0, typename I1, typename I2, typename I3>
  KOKKOS_FORCEINLINE_FUNCTION typename std::enable_if<
      (Kokkos::Impl::are_integral<I0, I1, I2, I3>::value && (4 == Rank) &&
       !is_default_map),
      reference_type>::type
  operator()(const I0& i0, const I1& i1, const I2& i2, const I3& i3) const {
    KOKKOS_IMPL_OFFSETVIEW_OPERATOR_VERIFY(
        (m_track, m_map, m_begins, i0, i1, i2, i3))
    const size_t j0 = i0 - m_begins[0];
    const size_t j1 = i1 - m_begins[1];
    const size_t j2 = i2 - m_begins[2];
    const size_t j3 = i3 - m_begins[3];
    return m_map.reference(j0, j1, j2, j3);
  }

  //------------------------------
  // Rank 5

  template <typename I0, typename I1, typename I2, typename I3, typename I4>
  KOKKOS_FORCEINLINE_FUNCTION typename std::enable_if<
      (Kokkos::Impl::are_integral<I0, I1, I2, I3, I4>::value && (5 == Rank) &&
       is_default_map),
      reference_type>::type
  operator()(const I0& i0, const I1& i1, const I2& i2, const I3& i3,
             const I4& i4) const {
    KOKKOS_IMPL_OFFSETVIEW_OPERATOR_VERIFY(
        (m_track, m_map, m_begins, i0, i1, i2, i3, i4))
    const size_t j0 = i0 - m_begins[0];
    const size_t j1 = i1 - m_begins[1];
    const size_t j2 = i2 - m_begins[2];
    const size_t j3 = i3 - m_begins[3];
    const size_t j4 = i4 - m_begins[4];
    return m_map.m_impl_handle[m_map.m_impl_offset(j0, j1, j2, j3, j4)];
  }

  template <typename I0, typename I1, typename I2, typename I3, typename I4>
  KOKKOS_FORCEINLINE_FUNCTION typename std::enable_if<
      (Kokkos::Impl::are_integral<I0, I1, I2, I3, I4>::value && (5 == Rank) &&
       !is_default_map),
      reference_type>::type
  operator()(const I0& i0, const I1& i1, const I2& i2, const I3& i3,
             const I4& i4) const {
    KOKKOS_IMPL_OFFSETVIEW_OPERATOR_VERIFY(
        (m_track, m_map, m_begins, i0, i1, i2, i3, i4))
    const size_t j0 = i0 - m_begins[0];
    const size_t j1 = i1 - m_begins[1];
    const size_t j2 = i2 - m_begins[2];
    const size_t j3 = i3 - m_begins[3];
    const size_t j4 = i4 - m_begins[4];
    return m_map.reference(j0, j1, j2, j3, j4);
  }

  //------------------------------
  // Rank 6

  template <typename I0, typename I1, typename I2, typename I3, typename I4,
            typename I5>
  KOKKOS_FORCEINLINE_FUNCTION typename std::enable_if<
      (Kokkos::Impl::are_integral<I0, I1, I2, I3, I4, I5>::value &&
       (6 == Rank) && is_default_map),
      reference_type>::type
  operator()(const I0& i0, const I1& i1, const I2& i2, const I3& i3,
             const I4& i4, const I5& i5) const {
    KOKKOS_IMPL_OFFSETVIEW_OPERATOR_VERIFY(
        (m_track, m_map, m_begins, i0, i1, i2, i3, i4, i5))
    const size_t j0 = i0 - m_begins[0];
    const size_t j1 = i1 - m_begins[1];
    const size_t j2 = i2 - m_begins[2];
    const size_t j3 = i3 - m_begins[3];
    const size_t j4 = i4 - m_begins[4];
    const size_t j5 = i5 - m_begins[5];
    return m_map.m_impl_handle[m_map.m_impl_offset(j0, j1, j2, j3, j4, j5)];
  }

  template <typename I0, typename I1, typename I2, typename I3, typename I4,
            typename I5>
  KOKKOS_FORCEINLINE_FUNCTION typename std::enable_if<
      (Kokkos::Impl::are_integral<I0, I1, I2, I3, I4, I5>::value &&
       (6 == Rank) && !is_default_map),
      reference_type>::type
  operator()(const I0& i0, const I1& i1, const I2& i2, const I3& i3,
             const I4& i4, const I5& i5) const {
    KOKKOS_IMPL_OFFSETVIEW_OPERATOR_VERIFY(
        (m_track, m_map, m_begins, i0, i1, i2, i3, i4, i5))
    const size_t j0 = i0 - m_begins[0];
    const size_t j1 = i1 - m_begins[1];
    const size_t j2 = i2 - m_begins[2];
    const size_t j3 = i3 - m_begins[3];
    const size_t j4 = i4 - m_begins[4];
    const size_t j5 = i5 - m_begins[5];
    return m_map.reference(j0, j1, j2, j3, j4, j5);
  }

  //------------------------------
  // Rank 7

  template <typename I0, typename I1, typename I2, typename I3, typename I4,
            typename I5, typename I6>
  KOKKOS_FORCEINLINE_FUNCTION typename std::enable_if<
      (Kokkos::Impl::are_integral<I0, I1, I2, I3, I4, I5, I6>::value &&
       (7 == Rank) && is_default_map),
      reference_type>::type
  operator()(const I0& i0, const I1& i1, const I2& i2, const I3& i3,
             const I4& i4, const I5& i5, const I6& i6) const {
    KOKKOS_IMPL_OFFSETVIEW_OPERATOR_VERIFY(
        (m_track, m_map, m_begins, i0, i1, i2, i3, i4, i5, i6))
    const size_t j0 = i0 - m_begins[0];
    const size_t j1 = i1 - m_begins[1];
    const size_t j2 = i2 - m_begins[2];
    const size_t j3 = i3 - m_begins[3];
    const size_t j4 = i4 - m_begins[4];
    const size_t j5 = i5 - m_begins[5];
    const size_t j6 = i6 - m_begins[6];
    return m_map.m_impl_handle[m_map.m_impl_offset(j0, j1, j2, j3, j4, j5, j6)];
  }

  template <typename I0, typename I1, typename I2, typename I3, typename I4,
            typename I5, typename I6>
  KOKKOS_FORCEINLINE_FUNCTION typename std::enable_if<
      (Kokkos::Impl::are_integral<I0, I1, I2, I3, I4, I5, I6>::value &&
       (7 == Rank) && !is_default_map),
      reference_type>::type
  operator()(const I0& i0, const I1& i1, const I2& i2, const I3& i3,
             const I4& i4, const I5& i5, const I6& i6) const {
    KOKKOS_IMPL_OFFSETVIEW_OPERATOR_VERIFY(
        (m_track, m_map, m_begins, i0, i1, i2, i3, i4, i5, i6))
    const size_t j0 = i0 - m_begins[0];
    const size_t j1 = i1 - m_begins[1];
    const size_t j2 = i2 - m_begins[2];
    const size_t j3 = i3 - m_begins[3];
    const size_t j4 = i4 - m_begins[4];
    const size_t j5 = i5 - m_begins[5];
    const size_t j6 = i6 - m_begins[6];
    return m_map.reference(j0, j1, j2, j3, j4, j5, j6);
  }

  //------------------------------
  // Rank 8

  template <typename I0, typename I1, typename I2, typename I3, typename I4,
            typename I5, typename I6, typename I7>
  KOKKOS_FORCEINLINE_FUNCTION typename std::enable_if<
      (Kokkos::Impl::are_integral<I0, I1, I2, I3, I4, I5, I6, I7>::value &&
       (8 == Rank) && is_default_map),
      reference_type>::type
  operator()(const I0& i0, const I1& i1, const I2& i2, const I3& i3,
             const I4& i4, const I5& i5, const I6& i6, const I7& i7) const {
    KOKKOS_IMPL_OFFSETVIEW_OPERATOR_VERIFY(
        (m_track, m_map, m_begins, i0, i1, i2, i3, i4, i5, i6, i7))
    const size_t j0 = i0 - m_begins[0];
    const size_t j1 = i1 - m_begins[1];
    const size_t j2 = i2 - m_begins[2];
    const size_t j3 = i3 - m_begins[3];
    const size_t j4 = i4 - m_begins[4];
    const size_t j5 = i5 - m_begins[5];
    const size_t j6 = i6 - m_begins[6];
    const size_t j7 = i7 - m_begins[7];
    return m_map
        .m_impl_handle[m_map.m_impl_offset(j0, j1, j2, j3, j4, j5, j6, j7)];
  }

  template <typename I0, typename I1, typename I2, typename I3, typename I4,
            typename I5, typename I6, typename I7>
  KOKKOS_FORCEINLINE_FUNCTION typename std::enable_if<
      (Kokkos::Impl::are_integral<I0, I1, I2, I3, I4, I5, I6, I7>::value &&
       (8 == Rank) && !is_default_map),
      reference_type>::type
  operator()(const I0& i0, const I1& i1, const I2& i2, const I3& i3,
             const I4& i4, const I5& i5, const I6& i6, const I7& i7) const {
    KOKKOS_IMPL_OFFSETVIEW_OPERATOR_VERIFY(
        (m_track, m_map, m_begins, i0, i1, i2, i3, i4, i5, i6, i7))
    const size_t j0 = i0 - m_begins[0];
    const size_t j1 = i1 - m_begins[1];
    const size_t j2 = i2 - m_begins[2];
    const size_t j3 = i3 - m_begins[3];
    const size_t j4 = i4 - m_begins[4];
    const size_t j5 = i5 - m_begins[5];
    const size_t j6 = i6 - m_begins[6];
    const size_t j7 = i7 - m_begins[7];
    return m_map.reference(j0, j1, j2, j3, j4, j5, j6, j7);
  }

#undef KOKKOS_IMPL_OFFSETVIEW_OPERATOR_VERIFY

  //----------------------------------------
  // Standard destructor, constructors, and assignment operators

  KOKKOS_DEFAULTED_FUNCTION
  ~OffsetView() = default;

  KOKKOS_INLINE_FUNCTION
  OffsetView() : m_track(), m_map() {
    for (size_t i = 0; i < Rank; ++i) m_begins[i] = KOKKOS_INVALID_OFFSET;
  }

  KOKKOS_INLINE_FUNCTION
  OffsetView(const OffsetView& rhs)
      : m_track(rhs.m_track, traits::is_managed),
        m_map(rhs.m_map),
        m_begins(rhs.m_begins) {}

  KOKKOS_INLINE_FUNCTION
  OffsetView(OffsetView&& rhs)
      : m_track(std::move(rhs.m_track)),
        m_map(std::move(rhs.m_map)),
        m_begins(std::move(rhs.m_begins)) {}

  KOKKOS_INLINE_FUNCTION
  OffsetView& operator=(const OffsetView& rhs) {
    m_track  = rhs.m_track;
    m_map    = rhs.m_map;
    m_begins = rhs.m_begins;
    return *this;
  }

  KOKKOS_INLINE_FUNCTION
  OffsetView& operator=(OffsetView&& rhs) {
    m_track  = std::move(rhs.m_track);
    m_map    = std::move(rhs.m_map);
    m_begins = std::move(rhs.m_begins);
    return *this;
  }

  // interoperability with View
 private:
  typedef View<typename traits::scalar_array_type,
               typename traits::array_layout, typename traits::device_type,
               typename traits::memory_traits>
      view_type;

 public:
  KOKKOS_INLINE_FUNCTION
  view_type view() const {
    view_type v(m_track, m_map);
    return v;
  }

  template <class RT, class... RP>
  KOKKOS_INLINE_FUNCTION OffsetView(const View<RT, RP...>& aview)
      : m_track(aview.impl_track()), m_map() {
    typedef typename OffsetView<RT, RP...>::traits SrcTraits;
    typedef Kokkos::Impl::ViewMapping<traits, SrcTraits, void> Mapping;
    static_assert(Mapping::is_assignable,
                  "Incompatible OffsetView copy construction");
    Mapping::assign(m_map, aview.impl_map(), m_track);

    for (int i = 0; i < aview.Rank; ++i) {
      m_begins[i] = 0;
    }
  }

  template <class RT, class... RP>
  KOKKOS_INLINE_FUNCTION OffsetView(const View<RT, RP...>& aview,
                                    const index_list_type& minIndices)
      : m_track(aview.impl_track()), m_map() {
    typedef typename OffsetView<RT, RP...>::traits SrcTraits;
    typedef Kokkos::Impl::ViewMapping<traits, SrcTraits, void> Mapping;
    static_assert(Mapping::is_assignable,
                  "Incompatible OffsetView copy construction");
    Mapping::assign(m_map, aview.impl_map(), m_track);

#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
    Kokkos::Experimental::Impl::runtime_check_rank_host(
        traits::rank_dynamic, Rank, minIndices, label());
#else
    Kokkos::Experimental::Impl::runtime_check_rank_device(traits::rank_dynamic,
                                                          Rank, minIndices);

#endif

    for (size_t i = 0; i < minIndices.size(); ++i) {
      m_begins[i] = minIndices.begin()[i];
    }
  }
  template <class RT, class... RP>
  KOKKOS_INLINE_FUNCTION OffsetView(const View<RT, RP...>& aview,
                                    const begins_type& beg)
      : m_track(aview.impl_track()), m_map(), m_begins(beg) {
    typedef typename OffsetView<RT, RP...>::traits SrcTraits;
    typedef Kokkos::Impl::ViewMapping<traits, SrcTraits, void> Mapping;
    static_assert(Mapping::is_assignable,
                  "Incompatible OffsetView copy construction");
    Mapping::assign(m_map, aview.impl_map(), m_track);

    //#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
    //        Kokkos::Experimental::Impl::runtime_check_rank_host(traits::rank_dynamic,
    //        Rank, minIndices, label());
    //#else
    //        Kokkos::Experimental::Impl::runtime_check_rank_device(traits::rank_dynamic,
    //        Rank, minIndices);
    //
    //#endif
  }

  // may assign unmanaged from managed.

  template <class RT, class... RP>
  KOKKOS_INLINE_FUNCTION OffsetView(const OffsetView<RT, RP...>& rhs)
      : m_track(rhs.m_track, traits::is_managed),
        m_map(),
        m_begins(rhs.m_begins) {
    typedef typename OffsetView<RT, RP...>::traits SrcTraits;
    typedef Kokkos::Impl::ViewMapping<traits, SrcTraits, void> Mapping;
    static_assert(Mapping::is_assignable,
                  "Incompatible OffsetView copy construction");
    Mapping::assign(m_map, rhs.m_map, rhs.m_track);  // swb what about assign?
  }

 private:
  enum class subtraction_failure {
    none,
    negative,
    overflow,
  };

  // Subtraction should return a non-negative number and not overflow
  KOKKOS_INLINE_FUNCTION static subtraction_failure check_subtraction(
      int64_t lhs, int64_t rhs) {
    if (lhs < rhs) return subtraction_failure::negative;

    if (static_cast<uint64_t>(-1) / static_cast<uint64_t>(2) <
        static_cast<uint64_t>(lhs) - static_cast<uint64_t>(rhs))
      return subtraction_failure::overflow;

    return subtraction_failure::none;
  }

  // Need a way to get at an element from both begins_type (aka Kokkos::Array
  // which doesn't have iterators) and index_list_type (aka
  // std::initializer_list which doesn't have .data() or operator[]).
  // Returns by value
  KOKKOS_INLINE_FUNCTION
  static int64_t at(const begins_type& a, size_t pos) { return a[pos]; }

  KOKKOS_INLINE_FUNCTION
  static int64_t at(index_list_type a, size_t pos) {
    return *(a.begin() + pos);
  }

#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
  // Check that begins < ends for all elements
  // B, E can be begins_type and/or index_list_type
  template <typename B, typename E>
  KOKKOS_INLINE_FUNCTION static subtraction_failure
  runtime_check_begins_ends_host(const B& begins, const E& ends) {
    std::string message;
    if (begins.size() != Rank)
      message +=
          "begins.size() "
          "(" +
          std::to_string(begins.size()) +
          ")"
          " != Rank "
          "(" +
          std::to_string(Rank) +
          ")"
          "\n";

    if (ends.size() != Rank)
      message +=
          "ends.size() "
          "(" +
          std::to_string(begins.size()) +
          ")"
          " != Rank "
          "(" +
          std::to_string(Rank) +
          ")"
          "\n";

    // If there are no errors so far, then rank == Rank
    // Otherwise, check as much as possible
    size_t rank = begins.size() < ends.size() ? begins.size() : ends.size();
    for (size_t i = 0; i != rank; ++i) {
      subtraction_failure sf = check_subtraction(at(ends, i), at(begins, i));
      if (sf != subtraction_failure::none) {
        message +=
            "("
            "ends[" +
            std::to_string(i) +
            "]"
            " "
            "(" +
            std::to_string(at(ends, i)) +
            ")"
            " - "
            "begins[" +
            std::to_string(i) +
            "]"
            " "
            "(" +
            std::to_string(at(begins, i)) +
            ")"
            ")";
        switch (sf) {
          case subtraction_failure::negative:
            message += " must be non-negative\n";
            break;
          case subtraction_failure::overflow: message += " overflows\n"; break;
          default: break;
        }
      }
    }

    if (!message.empty()) {
      message =
          "Kokkos::Experimental::OffsetView ERROR: for unmanaged OffsetView\n" +
          message;
      Kokkos::Impl::throw_runtime_exception(message);
    }

    return subtraction_failure::none;
  }
#endif  // KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST

  // Check the begins < ends for all elements
  template <typename B, typename E>
  KOKKOS_INLINE_FUNCTION static subtraction_failure
  runtime_check_begins_ends_device(const B& begins, const E& ends) {
    if (begins.size() != Rank)
      Kokkos::abort(
          "Kokkos::Experimental::OffsetView ERROR: for unmanaged "
          "OffsetView: begins has bad Rank");
    if (ends.size() != Rank)
      Kokkos::abort(
          "Kokkos::Experimental::OffsetView ERROR: for unmanaged "
          "OffsetView: ends has bad Rank");

    for (size_t i = 0; i != begins.size(); ++i) {
      switch (check_subtraction(at(ends, i), at(begins, i))) {
        case subtraction_failure::negative:
          Kokkos::abort(
              "Kokkos::Experimental::OffsetView ERROR: for unmanaged "
              "OffsetView: bad range");
          break;
        case subtraction_failure::overflow:
          Kokkos::abort(
              "Kokkos::Experimental::OffsetView ERROR: for unmanaged "
              "OffsetView: range overflows");
          break;
        default: break;
      }
    }

    return subtraction_failure::none;
  }

  // Constructor around unmanaged data after checking begins < ends for all
  // elements
  // Each of B, E can be begins_type and/or index_list_type
  // Precondition: begins.size() == ends.size() == m_begins.size() == Rank
  template <typename B, typename E>
  KOKKOS_INLINE_FUNCTION OffsetView(const pointer_type& p, const B& begins_,
                                    const E& ends_,
                                    subtraction_failure)
      : m_track()  // no tracking
        ,
        m_map(Kokkos::Impl::ViewCtorProp<pointer_type>(p),
              typename traits::array_layout(
                  Rank > 0 ? at(ends_, 0) - at(begins_, 0) : 0,
                  Rank > 1 ? at(ends_, 1) - at(begins_, 1) : 0,
                  Rank > 2 ? at(ends_, 2) - at(begins_, 2) : 0,
                  Rank > 3 ? at(ends_, 3) - at(begins_, 3) : 0,
                  Rank > 4 ? at(ends_, 4) - at(begins_, 4) : 0,
                  Rank > 5 ? at(ends_, 5) - at(begins_, 5) : 0,
                  Rank > 6 ? at(ends_, 6) - at(begins_, 6) : 0,
                  Rank > 7 ? at(ends_, 7) - at(begins_, 7) : 0)) {
    for (size_t i = 0; i != m_begins.size(); ++i) {
      m_begins[i] = at(begins_, i);
    };
  }

 public:
  // Constructor around unmanaged data
  // Four overloads, as both begins and ends can be either
  // begins_type or index_list_type
  KOKKOS_INLINE_FUNCTION
  OffsetView(const pointer_type& p, const begins_type& begins_,
             const begins_type& ends_)
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
      : OffsetView(p, begins_, ends_,
                   runtime_check_begins_ends_host(begins_, ends_))
#else
      : OffsetView(p, begins_, ends_,
                   runtime_check_begins_ends_device(begins_, ends_))
#endif
  {
  }

  KOKKOS_INLINE_FUNCTION
  OffsetView(const pointer_type& p, const begins_type& begins_,
             index_list_type ends_)
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
      : OffsetView(p, begins_, ends_,
                   runtime_check_begins_ends_host(begins_, ends_))
#else
      : OffsetView(p, begins_, ends_,
                   runtime_check_begins_ends_device(begins_, ends_))
#endif
  {
  }

  KOKKOS_INLINE_FUNCTION
  OffsetView(const pointer_type& p, index_list_type begins_,
             const begins_type& ends_)
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
      : OffsetView(p, begins_, ends_,
                   runtime_check_begins_ends_host(begins_, ends_))
#else
      : OffsetView(p, begins_, ends_,
                   runtime_check_begins_ends_device(begins_, ends_))
#endif
  {
  }

  KOKKOS_INLINE_FUNCTION
  OffsetView(const pointer_type& p, index_list_type begins_,
             index_list_type ends_)
#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
      : OffsetView(p, begins_, ends_,
                   runtime_check_begins_ends_host(begins_, ends_))
#else
      : OffsetView(p, begins_, ends_,
                   runtime_check_begins_ends_device(begins_, ends_))
#endif
  {
  }

  //----------------------------------------
  // Allocation tracking properties
  KOKKOS_INLINE_FUNCTION
  int use_count() const { return m_track.use_count(); }

  inline const std::string label() const {
    return m_track.template get_label<typename traits::memory_space>();
  }

  template <typename Label>
  explicit inline OffsetView(
      const Label& arg_label,
      typename std::enable_if<Kokkos::Impl::is_view_label<Label>::value,
                              const index_list_type>::type range0,
      const index_list_type range1 = KOKKOS_INVALID_INDEX_RANGE,
      const index_list_type range2 = KOKKOS_INVALID_INDEX_RANGE,
      const index_list_type range3 = KOKKOS_INVALID_INDEX_RANGE,
      const index_list_type range4 = KOKKOS_INVALID_INDEX_RANGE,
      const index_list_type range5 = KOKKOS_INVALID_INDEX_RANGE,
      const index_list_type range6 = KOKKOS_INVALID_INDEX_RANGE,
      const index_list_type range7 = KOKKOS_INVALID_INDEX_RANGE

      )
      : OffsetView(Kokkos::Impl::ViewCtorProp<std::string>(arg_label),
                   typename traits::array_layout(
                       range0.begin()[1] - range0.begin()[0] + 1,
                       range1.begin()[1] - range1.begin()[0] + 1,
                       range2.begin()[1] - range2.begin()[0] + 1,
                       range3.begin()[1] - range3.begin()[0] + 1,
                       range4.begin()[1] - range4.begin()[0] + 1,
                       range5.begin()[1] - range5.begin()[0] + 1,
                       range6.begin()[1] - range6.begin()[0] + 1,
                       range7.begin()[1] - range7.begin()[0] + 1),
                   {range0.begin()[0], range1.begin()[0], range2.begin()[0],
                    range3.begin()[0], range4.begin()[0], range5.begin()[0],
                    range6.begin()[0], range7.begin()[0]}) {}

  template <class... P>
  explicit KOKKOS_INLINE_FUNCTION OffsetView(
      const Kokkos::Impl::ViewCtorProp<P...>& arg_prop,
      typename std::enable_if<Kokkos::Impl::ViewCtorProp<P...>::has_pointer,
                              typename traits::array_layout>::type const&
          arg_layout,
      const index_list_type minIndices)
      : m_track()  // No memory tracking
        ,
        m_map(arg_prop, arg_layout) {
    for (size_t i = 0; i < minIndices.size(); ++i) {
      m_begins[i] = minIndices.begin()[i];
    }
    static_assert(
        std::is_same<pointer_type, typename Kokkos::Impl::ViewCtorProp<
                                       P...>::pointer_type>::value,
        "When constructing OffsetView to wrap user memory, you must supply "
        "matching pointer type");
  }

  template <class... P>
  explicit inline OffsetView(
      const Kokkos::Impl::ViewCtorProp<P...>& arg_prop,
      typename std::enable_if<!Kokkos::Impl::ViewCtorProp<P...>::has_pointer,
                              typename traits::array_layout>::type const&
          arg_layout,
      const index_list_type minIndices)
      : m_track(),
        m_map()

  {
    for (size_t i = 0; i < Rank; ++i) m_begins[i] = minIndices.begin()[i];

    // Append layout and spaces if not input
    typedef Kokkos::Impl::ViewCtorProp<P...> alloc_prop_input;

    // use 'std::integral_constant<unsigned,I>' for non-types
    // to avoid duplicate class error.
    typedef Kokkos::Impl::ViewCtorProp<
        P...,
        typename std::conditional<alloc_prop_input::has_label,
                                  std::integral_constant<unsigned, 0>,
                                  typename std::string>::type,
        typename std::conditional<
            alloc_prop_input::has_memory_space,
            std::integral_constant<unsigned, 1>,
            typename traits::device_type::memory_space>::type,
        typename std::conditional<
            alloc_prop_input::has_execution_space,
            std::integral_constant<unsigned, 2>,
            typename traits::device_type::execution_space>::type>
        alloc_prop;

    static_assert(traits::is_managed,
                  "OffsetView allocation constructor requires managed memory");

    if (alloc_prop::initialize &&
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
        !alloc_prop::execution_space::is_initialized()
#else
        !alloc_prop::execution_space::impl_is_initialized()
#endif
    ) {
      // If initializing view data then
      // the execution space must be initialized.
      Kokkos::Impl::throw_runtime_exception(
          "Constructing OffsetView and initializing data with uninitialized "
          "execution space");
    }

    // Copy the input allocation properties with possibly defaulted properties
    alloc_prop prop_copy(arg_prop);

    //------------------------------------------------------------
#if defined(KOKKOS_ENABLE_CUDA)
    // If allocating in CudaUVMSpace must fence before and after
    // the allocation to protect against possible concurrent access
    // on the CPU and the GPU.
    // Fence using the trait's executon space (which will be Kokkos::Cuda)
    // to avoid incomplete type errors from usng Kokkos::Cuda directly.
    if (std::is_same<Kokkos::CudaUVMSpace,
                     typename traits::device_type::memory_space>::value) {
      typename traits::device_type::memory_space::execution_space().fence();
    }
#endif
    //------------------------------------------------------------

    Kokkos::Impl::SharedAllocationRecord<>* record =
        m_map.allocate_shared(prop_copy, arg_layout);

    //------------------------------------------------------------
#if defined(KOKKOS_ENABLE_CUDA)
    if (std::is_same<Kokkos::CudaUVMSpace,
                     typename traits::device_type::memory_space>::value) {
      typename traits::device_type::memory_space::execution_space().fence();
    }
#endif
    //------------------------------------------------------------

    // Setup and initialization complete, start tracking
    m_track.assign_allocated_record_to_uninitialized(record);

#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
    Kokkos::Experimental::Impl::runtime_check_rank_host(
        traits::rank_dynamic, Rank, minIndices, label());
#else
    Kokkos::Experimental::Impl::runtime_check_rank_device(traits::rank_dynamic,
                                                          Rank, minIndices);

#endif
  }
};

/** \brief Temporary free function rank()
 *         until rank() is implemented
 *         in the View
 */
template <typename D, class... P>
KOKKOS_INLINE_FUNCTION constexpr unsigned rank(const OffsetView<D, P...>& V) {
  return V.Rank;
}  // Temporary until added to view

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
namespace Impl {

template <class T>
KOKKOS_INLINE_FUNCTION
    typename std::enable_if<std::is_integral<T>::value, T>::type
    shift_input(const T arg, const int64_t offset) {
  return arg - offset;
}

KOKKOS_INLINE_FUNCTION
Kokkos::Impl::ALL_t shift_input(const Kokkos::Impl::ALL_t arg,
                                const int64_t /*offset*/) {
  return arg;
}

template <class T>
KOKKOS_INLINE_FUNCTION typename std::enable_if<std::is_integral<T>::value,
                                               Kokkos::pair<T, T> >::type
shift_input(const Kokkos::pair<T, T> arg, const int64_t offset) {
  return Kokkos::make_pair<T, T>(arg.first - offset, arg.second - offset);
}
template <class T>
inline
    typename std::enable_if<std::is_integral<T>::value, std::pair<T, T> >::type
    shift_input(const std::pair<T, T> arg, const int64_t offset) {
  return std::make_pair<T, T>(arg.first - offset, arg.second - offset);
}

template <size_t N, class Arg, class A>
KOKKOS_INLINE_FUNCTION void map_arg_to_new_begin(
    const size_t i, Kokkos::Array<int64_t, N>& subviewBegins,
    typename std::enable_if<N != 0, const Arg>::type shiftedArg, const Arg arg,
    const A viewBegins, size_t& counter) {
  if (!std::is_integral<Arg>::value) {
    subviewBegins[counter] = shiftedArg == arg ? viewBegins[i] : 0;
    counter++;
  }
}

template <size_t N, class Arg, class A>
KOKKOS_INLINE_FUNCTION void map_arg_to_new_begin(
    const size_t /*i*/, Kokkos::Array<int64_t, N>& /*subviewBegins*/,
    typename std::enable_if<N == 0, const Arg>::type /*shiftedArg*/,
    const Arg /*arg*/, const A /*viewBegins*/, size_t& /*counter*/) {}

template <class D, class... P, class T>
KOKKOS_INLINE_FUNCTION
    typename Kokkos::Experimental::Impl::GetOffsetViewTypeFromViewType<
        typename Kokkos::Impl::ViewMapping<void /* deduce subview type from
                                                   source view traits */
                                           ,
                                           ViewTraits<D, P...>, T>::type>::type
    subview_offset(const OffsetView<D, P...>& src, T arg) {
  auto theView = src.view();
  auto begins  = src.begins();

  T shiftedArg = shift_input(arg, begins[0]);

  constexpr size_t rank =
      Kokkos::Impl::ViewMapping<void /* deduce subview type from source view
                                        traits */
                                ,
                                ViewTraits<D, P...>, T>::type::Rank;

  auto theSubview = Kokkos::subview(theView, shiftedArg);

  Kokkos::Array<int64_t, rank> subviewBegins;
  size_t counter = 0;
  Kokkos::Experimental::Impl::map_arg_to_new_begin(0, subviewBegins, shiftedArg,
                                                   arg, begins, counter);

  typename Kokkos::Experimental::Impl::GetOffsetViewTypeFromViewType<
      typename Kokkos::Impl::ViewMapping<void /* deduce subview type from source
                                                 view traits */
                                         ,
                                         ViewTraits<D, P...>, T>::type>::type
      offsetView(theSubview, subviewBegins);

  return offsetView;
}

template <class D, class... P, class T0, class T1>
KOKKOS_INLINE_FUNCTION
    typename Kokkos::Experimental::Impl::GetOffsetViewTypeFromViewType<
        typename Kokkos::Impl::ViewMapping<
            void /* deduce subview type from source view traits */
            ,
            ViewTraits<D, P...>, T0, T1>::type>::type
    subview_offset(const Kokkos::Experimental::OffsetView<D, P...>& src,
                   T0 arg0, T1 arg1) {
  auto theView = src.view();
  auto begins  = src.begins();

  T0 shiftedArg0 = shift_input(arg0, begins[0]);
  T1 shiftedArg1 = shift_input(arg1, begins[1]);

  auto theSubview = Kokkos::subview(theView, shiftedArg0, shiftedArg1);
  constexpr size_t rank =
      Kokkos::Impl::ViewMapping<void /* deduce subview type from source view
                                        traits */
                                ,
                                ViewTraits<D, P...>, T0, T1>::type::Rank;

  Kokkos::Array<int64_t, rank> subviewBegins;
  size_t counter = 0;
  Kokkos::Experimental::Impl::map_arg_to_new_begin(
      0, subviewBegins, shiftedArg0, arg0, begins, counter);
  Kokkos::Experimental::Impl::map_arg_to_new_begin(
      1, subviewBegins, shiftedArg1, arg1, begins, counter);

  typename Kokkos::Experimental::Impl::GetOffsetViewTypeFromViewType<
      typename Kokkos::Impl::ViewMapping<
          void /* deduce subview type from source view traits */
          ,
          ViewTraits<D, P...>, T0, T1>::type>::type offsetView(theSubview,
                                                               subviewBegins);

  return offsetView;
}

template <class D, class... P, class T0, class T1, class T2>
KOKKOS_INLINE_FUNCTION
    typename Kokkos::Experimental::Impl::GetOffsetViewTypeFromViewType<
        typename Kokkos::Impl::ViewMapping<
            void /* deduce subview type from source view traits */
            ,
            ViewTraits<D, P...>, T0, T1, T2>::type>::type
    subview_offset(const OffsetView<D, P...>& src, T0 arg0, T1 arg1, T2 arg2) {
  auto theView = src.view();
  auto begins  = src.begins();

  T0 shiftedArg0 = shift_input(arg0, begins[0]);
  T1 shiftedArg1 = shift_input(arg1, begins[1]);
  T2 shiftedArg2 = shift_input(arg2, begins[2]);

  auto theSubview =
      Kokkos::subview(theView, shiftedArg0, shiftedArg1, shiftedArg2);

  constexpr size_t rank =
      Kokkos::Impl::ViewMapping<void /* deduce subview type from source view
                                        traits */
                                ,
                                ViewTraits<D, P...>, T0, T1, T2>::type::Rank;

  Kokkos::Array<int64_t, rank> subviewBegins;

  size_t counter = 0;
  Kokkos::Experimental::Impl::map_arg_to_new_begin(
      0, subviewBegins, shiftedArg0, arg0, begins, counter);
  Kokkos::Experimental::Impl::map_arg_to_new_begin(
      1, subviewBegins, shiftedArg1, arg1, begins, counter);
  Kokkos::Experimental::Impl::map_arg_to_new_begin(
      2, subviewBegins, shiftedArg2, arg2, begins, counter);

  typename Kokkos::Experimental::Impl::GetOffsetViewTypeFromViewType<
      typename Kokkos::Impl::ViewMapping<
          void /* deduce subview type from source view traits */
          ,
          ViewTraits<D, P...>, T0, T1, T2>::type>::type
      offsetView(theSubview, subviewBegins);

  return offsetView;
}

template <class D, class... P, class T0, class T1, class T2, class T3>
KOKKOS_INLINE_FUNCTION
    typename Kokkos::Experimental::Impl::GetOffsetViewTypeFromViewType<
        typename Kokkos::Impl::ViewMapping<
            void /* deduce subview type from source view traits */
            ,
            ViewTraits<D, P...>, T0, T1, T2, T3>::type>::type
    subview_offset(const OffsetView<D, P...>& src, T0 arg0, T1 arg1, T2 arg2,
                   T3 arg3) {
  auto theView = src.view();
  auto begins  = src.begins();

  T0 shiftedArg0 = shift_input(arg0, begins[0]);
  T1 shiftedArg1 = shift_input(arg1, begins[1]);
  T2 shiftedArg2 = shift_input(arg2, begins[2]);
  T3 shiftedArg3 = shift_input(arg3, begins[3]);

  auto theSubview = Kokkos::subview(theView, shiftedArg0, shiftedArg1,
                                    shiftedArg2, shiftedArg3);

  constexpr size_t rank = Kokkos::Impl::ViewMapping<
      void /* deduce subview type from source view traits */
      ,
      ViewTraits<D, P...>, T0, T1, T2, T3>::type::Rank;
  Kokkos::Array<int64_t, rank> subviewBegins;

  size_t counter = 0;
  Kokkos::Experimental::Impl::map_arg_to_new_begin(
      0, subviewBegins, shiftedArg0, arg0, begins, counter);
  Kokkos::Experimental::Impl::map_arg_to_new_begin(
      1, subviewBegins, shiftedArg1, arg1, begins, counter);
  Kokkos::Experimental::Impl::map_arg_to_new_begin(
      2, subviewBegins, shiftedArg2, arg2, begins, counter);
  Kokkos::Experimental::Impl::map_arg_to_new_begin(
      3, subviewBegins, shiftedArg3, arg3, begins, counter);

  typename Kokkos::Experimental::Impl::GetOffsetViewTypeFromViewType<
      typename Kokkos::Impl::ViewMapping<
          void /* deduce subview type from source view traits */
          ,
          ViewTraits<D, P...>, T0, T1, T2, T3>::type>::type
      offsetView(theSubview, subviewBegins);

  return offsetView;
}

template <class D, class... P, class T0, class T1, class T2, class T3, class T4>
KOKKOS_INLINE_FUNCTION
    typename Kokkos::Experimental::Impl::GetOffsetViewTypeFromViewType<
        typename Kokkos::Impl::ViewMapping<
            void /* deduce subview type from source view traits */
            ,
            ViewTraits<D, P...>, T0, T1, T2, T3, T4>::type>::type
    subview_offset(const OffsetView<D, P...>& src, T0 arg0, T1 arg1, T2 arg2,
                   T3 arg3, T4 arg4) {
  auto theView = src.view();
  auto begins  = src.begins();

  T0 shiftedArg0 = shift_input(arg0, begins[0]);
  T1 shiftedArg1 = shift_input(arg1, begins[1]);
  T2 shiftedArg2 = shift_input(arg2, begins[2]);
  T3 shiftedArg3 = shift_input(arg3, begins[3]);
  T4 shiftedArg4 = shift_input(arg4, begins[4]);

  auto theSubview = Kokkos::subview(theView, shiftedArg0, shiftedArg1,
                                    shiftedArg2, shiftedArg3, shiftedArg4);

  constexpr size_t rank = Kokkos::Impl::ViewMapping<
      void /* deduce subview type from source view traits */
      ,
      ViewTraits<D, P...>, T0, T1, T2, T3, T4>::type::Rank;
  Kokkos::Array<int64_t, rank> subviewBegins;

  size_t counter = 0;
  Kokkos::Experimental::Impl::map_arg_to_new_begin(
      0, subviewBegins, shiftedArg0, arg0, begins, counter);
  Kokkos::Experimental::Impl::map_arg_to_new_begin(
      1, subviewBegins, shiftedArg1, arg1, begins, counter);
  Kokkos::Experimental::Impl::map_arg_to_new_begin(
      2, subviewBegins, shiftedArg2, arg2, begins, counter);
  Kokkos::Experimental::Impl::map_arg_to_new_begin(
      3, subviewBegins, shiftedArg3, arg3, begins, counter);
  Kokkos::Experimental::Impl::map_arg_to_new_begin(
      4, subviewBegins, shiftedArg4, arg4, begins, counter);

  typename Kokkos::Experimental::Impl::GetOffsetViewTypeFromViewType<
      typename Kokkos::Impl::ViewMapping<
          void /* deduce subview type from source view traits */
          ,
          ViewTraits<D, P...>, T0, T1, T2, T3, T4>::type>::type
      offsetView(theSubview, subviewBegins);

  return offsetView;
}

template <class D, class... P, class T0, class T1, class T2, class T3, class T4,
          class T5>
KOKKOS_INLINE_FUNCTION
    typename Kokkos::Experimental::Impl::GetOffsetViewTypeFromViewType<
        typename Kokkos::Impl::ViewMapping<
            void /* deduce subview type from source view traits */
            ,
            ViewTraits<D, P...>, T0, T1, T2, T3, T4, T5>::type>::type
    subview_offset(const OffsetView<D, P...>& src, T0 arg0, T1 arg1, T2 arg2,
                   T3 arg3, T4 arg4, T5 arg5) {
  auto theView = src.view();
  auto begins  = src.begins();

  T0 shiftedArg0 = shift_input(arg0, begins[0]);
  T1 shiftedArg1 = shift_input(arg1, begins[1]);
  T2 shiftedArg2 = shift_input(arg2, begins[2]);
  T3 shiftedArg3 = shift_input(arg3, begins[3]);
  T4 shiftedArg4 = shift_input(arg4, begins[4]);
  T5 shiftedArg5 = shift_input(arg5, begins[5]);

  auto theSubview =
      Kokkos::subview(theView, shiftedArg0, shiftedArg1, shiftedArg2,
                      shiftedArg3, shiftedArg4, shiftedArg5);

  constexpr size_t rank = Kokkos::Impl::ViewMapping<
      void /* deduce subview type from source view traits */
      ,
      ViewTraits<D, P...>, T0, T1, T2, T3, T4, T5>::type::Rank;

  Kokkos::Array<int64_t, rank> subviewBegins;

  size_t counter = 0;
  Kokkos::Experimental::Impl::map_arg_to_new_begin(
      0, subviewBegins, shiftedArg0, arg0, begins, counter);
  Kokkos::Experimental::Impl::map_arg_to_new_begin(
      1, subviewBegins, shiftedArg1, arg1, begins, counter);
  Kokkos::Experimental::Impl::map_arg_to_new_begin(
      2, subviewBegins, shiftedArg2, arg2, begins, counter);
  Kokkos::Experimental::Impl::map_arg_to_new_begin(
      3, subviewBegins, shiftedArg3, arg3, begins, counter);
  Kokkos::Experimental::Impl::map_arg_to_new_begin(
      4, subviewBegins, shiftedArg4, arg4, begins, counter);
  Kokkos::Experimental::Impl::map_arg_to_new_begin(
      5, subviewBegins, shiftedArg5, arg5, begins, counter);

  typename Kokkos::Experimental::Impl::GetOffsetViewTypeFromViewType<
      typename Kokkos::Impl::ViewMapping<
          void /* deduce subview type from source view traits */
          ,
          ViewTraits<D, P...>, T0, T1, T2, T3, T4, T5>::type>::type
      offsetView(theSubview, subviewBegins);

  return offsetView;
}
template <class D, class... P, class T0, class T1, class T2, class T3, class T4,
          class T5, class T6>
KOKKOS_INLINE_FUNCTION
    typename Kokkos::Experimental::Impl::GetOffsetViewTypeFromViewType<
        typename Kokkos::Impl::ViewMapping<
            void /* deduce subview type from source view traits */
            ,
            ViewTraits<D, P...>, T0, T1, T2, T3, T4, T5, T6>::type>::type
    subview_offset(const OffsetView<D, P...>& src, T0 arg0, T1 arg1, T2 arg2,
                   T3 arg3, T4 arg4, T5 arg5, T6 arg6) {
  auto theView = src.view();
  auto begins  = src.begins();

  T0 shiftedArg0 = shift_input(arg0, begins[0]);
  T1 shiftedArg1 = shift_input(arg1, begins[1]);
  T2 shiftedArg2 = shift_input(arg2, begins[2]);
  T3 shiftedArg3 = shift_input(arg3, begins[3]);
  T4 shiftedArg4 = shift_input(arg4, begins[4]);
  T5 shiftedArg5 = shift_input(arg5, begins[5]);
  T6 shiftedArg6 = shift_input(arg6, begins[6]);

  auto theSubview =
      Kokkos::subview(theView, shiftedArg0, shiftedArg1, shiftedArg2,
                      shiftedArg3, shiftedArg4, shiftedArg5, shiftedArg6);

  constexpr size_t rank = Kokkos::Impl::ViewMapping<
      void /* deduce subview type from source view traits */
      ,
      ViewTraits<D, P...>, T0, T1, T2, T3, T4, T5, T6>::type::Rank;

  Kokkos::Array<int64_t, rank> subviewBegins;

  size_t counter = 0;
  Kokkos::Experimental::Impl::map_arg_to_new_begin(
      0, subviewBegins, shiftedArg0, arg0, begins, counter);
  Kokkos::Experimental::Impl::map_arg_to_new_begin(
      1, subviewBegins, shiftedArg1, arg1, begins, counter);
  Kokkos::Experimental::Impl::map_arg_to_new_begin(
      2, subviewBegins, shiftedArg2, arg2, begins, counter);
  Kokkos::Experimental::Impl::map_arg_to_new_begin(
      3, subviewBegins, shiftedArg3, arg3, begins, counter);
  Kokkos::Experimental::Impl::map_arg_to_new_begin(
      4, subviewBegins, shiftedArg4, arg4, begins, counter);
  Kokkos::Experimental::Impl::map_arg_to_new_begin(
      5, subviewBegins, shiftedArg5, arg5, begins, counter);
  Kokkos::Experimental::Impl::map_arg_to_new_begin(
      6, subviewBegins, shiftedArg6, arg6, begins, counter);

  typename Kokkos::Experimental::Impl::GetOffsetViewTypeFromViewType<
      typename Kokkos::Impl::ViewMapping<
          void /* deduce subview type from source view traits */
          ,
          ViewTraits<D, P...>, T0, T1, T2, T3, T4, T5, T6>::type>::type
      offsetView(theSubview, subviewBegins);

  return offsetView;
}

template <class D, class... P, class T0, class T1, class T2, class T3, class T4,
          class T5, class T6, class T7>
KOKKOS_INLINE_FUNCTION
    typename Kokkos::Experimental::Impl::GetOffsetViewTypeFromViewType<
        typename Kokkos::Impl::ViewMapping<
            void /* deduce subview type from source view traits */
            ,
            ViewTraits<D, P...>, T0, T1, T2, T3, T4, T5, T6, T7>::type>::type
    subview_offset(const OffsetView<D, P...>& src, T0 arg0, T1 arg1, T2 arg2,
                   T3 arg3, T4 arg4, T5 arg5, T6 arg6, T7 arg7) {
  auto theView = src.view();
  auto begins  = src.begins();

  T0 shiftedArg0 = shift_input(arg0, begins[0]);
  T1 shiftedArg1 = shift_input(arg1, begins[1]);
  T2 shiftedArg2 = shift_input(arg2, begins[2]);
  T3 shiftedArg3 = shift_input(arg3, begins[3]);
  T4 shiftedArg4 = shift_input(arg4, begins[4]);
  T5 shiftedArg5 = shift_input(arg5, begins[5]);
  T6 shiftedArg6 = shift_input(arg6, begins[6]);
  T7 shiftedArg7 = shift_input(arg7, begins[7]);

  auto theSubview = Kokkos::subview(theView, shiftedArg0, shiftedArg1,
                                    shiftedArg2, shiftedArg3, shiftedArg4,
                                    shiftedArg5, shiftedArg6, shiftedArg7);

  constexpr size_t rank = Kokkos::Impl::ViewMapping<
      void /* deduce subview type from source view traits */
      ,
      ViewTraits<D, P...>, T0, T1, T2, T3, T4, T5, T6, T7>::type::Rank;

  Kokkos::Array<int64_t, rank> subviewBegins;

  size_t counter = 0;
  Kokkos::Experimental::Impl::map_arg_to_new_begin(
      0, subviewBegins, shiftedArg0, arg0, begins, counter);
  Kokkos::Experimental::Impl::map_arg_to_new_begin(
      1, subviewBegins, shiftedArg1, arg1, begins, counter);
  Kokkos::Experimental::Impl::map_arg_to_new_begin(
      2, subviewBegins, shiftedArg2, arg2, begins, counter);
  Kokkos::Experimental::Impl::map_arg_to_new_begin(
      3, subviewBegins, shiftedArg3, arg3, begins, counter);
  Kokkos::Experimental::Impl::map_arg_to_new_begin(
      4, subviewBegins, shiftedArg4, arg4, begins, counter);
  Kokkos::Experimental::Impl::map_arg_to_new_begin(
      5, subviewBegins, shiftedArg5, arg5, begins, counter);
  Kokkos::Experimental::Impl::map_arg_to_new_begin(
      6, subviewBegins, shiftedArg6, arg6, begins, counter);
  Kokkos::Experimental::Impl::map_arg_to_new_begin(
      7, subviewBegins, shiftedArg7, arg7, begins, counter);

  typename Kokkos::Experimental::Impl::GetOffsetViewTypeFromViewType<
      typename Kokkos::Impl::ViewMapping<
          void /* deduce subview type from source view traits */
          ,
          ViewTraits<D, P...>, T0, T1, T2, T3, T4, T5, T6, T7>::type>::type
      offsetView(theSubview, subviewBegins);

  return offsetView;
}
}  // namespace Impl

template <class D, class... P, class... Args>
KOKKOS_INLINE_FUNCTION
    typename Kokkos::Experimental::Impl::GetOffsetViewTypeFromViewType<
        typename Kokkos::Impl::ViewMapping<
            void /* deduce subview type from source view traits */
            ,
            ViewTraits<D, P...>, Args...>::type>::type
    subview(const OffsetView<D, P...>& src, Args... args) {
  static_assert(
      OffsetView<D, P...>::Rank == sizeof...(Args),
      "subview requires one argument for each source OffsetView rank");

  return Kokkos::Experimental::Impl::subview_offset(src, args...);
}

}  // namespace Experimental
}  // namespace Kokkos
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {
template <class LT, class... LP, class RT, class... RP>
KOKKOS_INLINE_FUNCTION bool operator==(const OffsetView<LT, LP...>& lhs,
                                       const OffsetView<RT, RP...>& rhs) {
  // Same data, layout, dimensions
  typedef ViewTraits<LT, LP...> lhs_traits;
  typedef ViewTraits<RT, RP...> rhs_traits;

  return std::is_same<typename lhs_traits::const_value_type,
                      typename rhs_traits::const_value_type>::value &&
         std::is_same<typename lhs_traits::array_layout,
                      typename rhs_traits::array_layout>::value &&
         std::is_same<typename lhs_traits::memory_space,
                      typename rhs_traits::memory_space>::value &&
         unsigned(lhs_traits::rank) == unsigned(rhs_traits::rank) &&
         lhs.data() == rhs.data() && lhs.span() == rhs.span() &&
         lhs.extent(0) == rhs.extent(0) && lhs.extent(1) == rhs.extent(1) &&
         lhs.extent(2) == rhs.extent(2) && lhs.extent(3) == rhs.extent(3) &&
         lhs.extent(4) == rhs.extent(4) && lhs.extent(5) == rhs.extent(5) &&
         lhs.extent(6) == rhs.extent(6) && lhs.extent(7) == rhs.extent(7) &&
         lhs.begin(0) == rhs.begin(0) && lhs.begin(1) == rhs.begin(1) &&
         lhs.begin(2) == rhs.begin(2) && lhs.begin(3) == rhs.begin(3) &&
         lhs.begin(4) == rhs.begin(4) && lhs.begin(5) == rhs.begin(5) &&
         lhs.begin(6) == rhs.begin(6) && lhs.begin(7) == rhs.begin(7);
}

template <class LT, class... LP, class RT, class... RP>
KOKKOS_INLINE_FUNCTION bool operator!=(const OffsetView<LT, LP...>& lhs,
                                       const OffsetView<RT, RP...>& rhs) {
  return !(operator==(lhs, rhs));
}

template <class LT, class... LP, class RT, class... RP>
KOKKOS_INLINE_FUNCTION bool operator==(const View<LT, LP...>& lhs,
                                       const OffsetView<RT, RP...>& rhs) {
  // Same data, layout, dimensions
  typedef ViewTraits<LT, LP...> lhs_traits;
  typedef ViewTraits<RT, RP...> rhs_traits;

  return std::is_same<typename lhs_traits::const_value_type,
                      typename rhs_traits::const_value_type>::value &&
         std::is_same<typename lhs_traits::array_layout,
                      typename rhs_traits::array_layout>::value &&
         std::is_same<typename lhs_traits::memory_space,
                      typename rhs_traits::memory_space>::value &&
         unsigned(lhs_traits::rank) == unsigned(rhs_traits::rank) &&
         lhs.data() == rhs.data() && lhs.span() == rhs.span() &&
         lhs.extent(0) == rhs.extent(0) && lhs.extent(1) == rhs.extent(1) &&
         lhs.extent(2) == rhs.extent(2) && lhs.extent(3) == rhs.extent(3) &&
         lhs.extent(4) == rhs.extent(4) && lhs.extent(5) == rhs.extent(5) &&
         lhs.extent(6) == rhs.extent(6) && lhs.extent(7) == rhs.extent(7);
}

template <class LT, class... LP, class RT, class... RP>
KOKKOS_INLINE_FUNCTION bool operator==(const OffsetView<LT, LP...>& lhs,
                                       const View<RT, RP...>& rhs) {
  return rhs == lhs;
}

}  // namespace Experimental
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Experimental {
template <class DT, class... DP>
inline void deep_copy(
    const OffsetView<DT, DP...>& dst,
    typename ViewTraits<DT, DP...>::const_value_type& value,
    typename std::enable_if<std::is_same<
        typename ViewTraits<DT, DP...>::specialize, void>::value>::type* =
        nullptr) {
  static_assert(
      std::is_same<typename ViewTraits<DT, DP...>::non_const_value_type,
                   typename ViewTraits<DT, DP...>::value_type>::value,
      "deep_copy requires non-const type");

  auto dstView = dst.view();
  Kokkos::deep_copy(dstView, value);
}

template <class DT, class... DP, class ST, class... SP>
inline void deep_copy(
    const OffsetView<DT, DP...>& dst, const OffsetView<ST, SP...>& value,
    typename std::enable_if<std::is_same<
        typename ViewTraits<DT, DP...>::specialize, void>::value>::type* =
        nullptr) {
  static_assert(
      std::is_same<typename ViewTraits<DT, DP...>::value_type,
                   typename ViewTraits<ST, SP...>::non_const_value_type>::value,
      "deep_copy requires matching non-const destination type");

  auto dstView = dst.view();
  Kokkos::deep_copy(dstView, value.view());
}
template <class DT, class... DP, class ST, class... SP>
inline void deep_copy(
    const OffsetView<DT, DP...>& dst, const View<ST, SP...>& value,
    typename std::enable_if<std::is_same<
        typename ViewTraits<DT, DP...>::specialize, void>::value>::type* =
        nullptr) {
  static_assert(
      std::is_same<typename ViewTraits<DT, DP...>::value_type,
                   typename ViewTraits<ST, SP...>::non_const_value_type>::value,
      "deep_copy requires matching non-const destination type");

  auto dstView = dst.view();
  Kokkos::deep_copy(dstView, value);
}

template <class DT, class... DP, class ST, class... SP>
inline void deep_copy(
    const View<DT, DP...>& dst, const OffsetView<ST, SP...>& value,
    typename std::enable_if<std::is_same<
        typename ViewTraits<DT, DP...>::specialize, void>::value>::type* =
        nullptr) {
  static_assert(
      std::is_same<typename ViewTraits<DT, DP...>::value_type,
                   typename ViewTraits<ST, SP...>::non_const_value_type>::value,
      "deep_copy requires matching non-const destination type");

  Kokkos::deep_copy(dst, value.view());
}
namespace Impl {

// Deduce Mirror Types
template <class Space, class T, class... P>
struct MirrorOffsetViewType {
  // The incoming view_type
  typedef typename Kokkos::Experimental::OffsetView<T, P...> src_view_type;
  // The memory space for the mirror view
  typedef typename Space::memory_space memory_space;
  // Check whether it is the same memory space
  enum {
    is_same_memspace =
        std::is_same<memory_space, typename src_view_type::memory_space>::value
  };
  // The array_layout
  typedef typename src_view_type::array_layout array_layout;
  // The data type (we probably want it non-const since otherwise we can't even
  // deep_copy to it.
  typedef typename src_view_type::non_const_data_type data_type;
  // The destination view type if it is not the same memory space
  typedef Kokkos::Experimental::OffsetView<data_type, array_layout, Space>
      dest_view_type;
  // If it is the same memory_space return the existsing view_type
  // This will also keep the unmanaged trait if necessary
  typedef typename std::conditional<is_same_memspace, src_view_type,
                                    dest_view_type>::type view_type;
};

template <class Space, class T, class... P>
struct MirrorOffsetType {
  // The incoming view_type
  typedef typename Kokkos::Experimental::OffsetView<T, P...> src_view_type;
  // The memory space for the mirror view
  typedef typename Space::memory_space memory_space;
  // Check whether it is the same memory space
  enum {
    is_same_memspace =
        std::is_same<memory_space, typename src_view_type::memory_space>::value
  };
  // The array_layout
  typedef typename src_view_type::array_layout array_layout;
  // The data type (we probably want it non-const since otherwise we can't even
  // deep_copy to it.
  typedef typename src_view_type::non_const_data_type data_type;
  // The destination view type if it is not the same memory space
  typedef Kokkos::Experimental::OffsetView<data_type, array_layout, Space>
      view_type;
};

}  // namespace Impl

template <class T, class... P>
inline typename Kokkos::Experimental::OffsetView<T, P...>::HostMirror
create_mirror(
    const Kokkos::Experimental::OffsetView<T, P...>& src,
    typename std::enable_if<
        !std::is_same<typename Kokkos::ViewTraits<T, P...>::array_layout,
                      Kokkos::LayoutStride>::value>::type* = 0) {
  typedef OffsetView<T, P...> src_type;
  typedef typename src_type::HostMirror dst_type;

  return dst_type(
      Kokkos::Impl::ViewCtorProp<std::string>(
          std::string(src.label()).append("_mirror")),
      typename Kokkos::ViewTraits<T, P...>::array_layout(
          src.extent(0), src.extent(1), src.extent(2), src.extent(3),
          src.extent(4), src.extent(5), src.extent(6), src.extent(7)),
      {src.begin(0), src.begin(1), src.begin(2), src.begin(3), src.begin(4),
       src.begin(5), src.begin(6), src.begin(7)});
}

template <class T, class... P>
inline typename Kokkos::Experimental::OffsetView<T, P...>::HostMirror
create_mirror(
    const Kokkos::Experimental::OffsetView<T, P...>& src,
    typename std::enable_if<
        std::is_same<typename Kokkos::ViewTraits<T, P...>::array_layout,
                     Kokkos::LayoutStride>::value>::type* = 0) {
  typedef OffsetView<T, P...> src_type;
  typedef typename src_type::HostMirror dst_type;

  Kokkos::LayoutStride layout;

  layout.dimension[0] = src.extent(0);
  layout.dimension[1] = src.extent(1);
  layout.dimension[2] = src.extent(2);
  layout.dimension[3] = src.extent(3);
  layout.dimension[4] = src.extent(4);
  layout.dimension[5] = src.extent(5);
  layout.dimension[6] = src.extent(6);
  layout.dimension[7] = src.extent(7);

  layout.stride[0] = src.stride_0();
  layout.stride[1] = src.stride_1();
  layout.stride[2] = src.stride_2();
  layout.stride[3] = src.stride_3();
  layout.stride[4] = src.stride_4();
  layout.stride[5] = src.stride_5();
  layout.stride[6] = src.stride_6();
  layout.stride[7] = src.stride_7();

  return dst_type(std::string(src.label()).append("_mirror"), layout,
                  {src.begin(0), src.begin(1), src.begin(2), src.begin(3),
                   src.begin(4), src.begin(5), src.begin(6), src.begin(7)});
}

// Create a mirror in a new space (specialization for different space)
template <class Space, class T, class... P>
typename Kokkos::Experimental::Impl::MirrorOffsetType<Space, T, P...>::view_type
create_mirror(const Space&,
              const Kokkos::Experimental::OffsetView<T, P...>& src) {
  return typename Kokkos::Experimental::Impl::MirrorOffsetType<
      Space, T, P...>::view_type(src.label(), src.layout(),
                                 {src.begin(0), src.begin(1), src.begin(2),
                                  src.begin(3), src.begin(4), src.begin(5),
                                  src.begin(6), src.begin(7)});
}

template <class T, class... P>
inline typename Kokkos::Experimental::OffsetView<T, P...>::HostMirror
create_mirror_view(
    const typename Kokkos::Experimental::OffsetView<T, P...>& src,
    typename std::enable_if<
        (std::is_same<
             typename Kokkos::Experimental::OffsetView<T, P...>::memory_space,
             typename Kokkos::Experimental::OffsetView<
                 T, P...>::HostMirror::memory_space>::value &&
         std::is_same<
             typename Kokkos::Experimental::OffsetView<T, P...>::data_type,
             typename Kokkos::Experimental::OffsetView<
                 T, P...>::HostMirror::data_type>::value)>::type* = nullptr) {
  return src;
}

template <class T, class... P>
inline typename Kokkos::Experimental::OffsetView<T, P...>::HostMirror
create_mirror_view(
    const Kokkos::Experimental::OffsetView<T, P...>& src,
    typename std::enable_if<
        !(std::is_same<
              typename Kokkos::Experimental::OffsetView<T, P...>::memory_space,
              typename Kokkos::Experimental::OffsetView<
                  T, P...>::HostMirror::memory_space>::value &&
          std::is_same<
              typename Kokkos::Experimental::OffsetView<T, P...>::data_type,
              typename Kokkos::Experimental::OffsetView<
                  T, P...>::HostMirror::data_type>::value)>::type* = 0) {
  return Kokkos::Experimental::create_mirror(src);
}

// Create a mirror view in a new space (specialization for same space)
template <class Space, class T, class... P>
typename Kokkos::Experimental::Impl::MirrorOffsetViewType<Space, T,
                                                          P...>::view_type
create_mirror_view(const Space&,
                   const Kokkos::Experimental::OffsetView<T, P...>& src,
                   typename std::enable_if<Impl::MirrorOffsetViewType<
                       Space, T, P...>::is_same_memspace>::type* = 0) {
  return src;
}

// Create a mirror view in a new space (specialization for different space)
template <class Space, class T, class... P>
typename Kokkos::Experimental::Impl::MirrorOffsetViewType<Space, T,
                                                          P...>::view_type
create_mirror_view(const Space&,
                   const Kokkos::Experimental::OffsetView<T, P...>& src,
                   typename std::enable_if<!Impl::MirrorOffsetViewType<
                       Space, T, P...>::is_same_memspace>::type* = 0) {
  return typename Kokkos::Experimental::Impl::MirrorOffsetViewType<
      Space, T, P...>::view_type(src.label(), src.layout(),
                                 {src.begin(0), src.begin(1), src.begin(2),
                                  src.begin(3), src.begin(4), src.begin(5),
                                  src.begin(6), src.begin(7)});
}
//
//  // Create a mirror view and deep_copy in a new space (specialization for
//  same space) template<class Space, class T, class ... P> typename
//  Kokkos::Experimental::Impl::MirrorViewType<Space,T,P ...>::view_type
//  create_mirror_view_and_copy(const Space& , const
//  Kokkos::Experimental::OffsetView<T,P...> & src
//                              , std::string const& name = ""
//                                  , typename
//                                  std::enable_if<Impl::MirrorViewType<Space,T,P
//                                  ...>::is_same_memspace>::type* = 0 ) {
//    (void)name;
//    return src;
//  }
//
//  // Create a mirror view and deep_copy in a new space (specialization for
//  different space) template<class Space, class T, class ... P> typename
//  Kokkos::Experimental::Impl::MirrorViewType<Space,T,P ...>::view_type
//  create_mirror_view_and_copy(const Space& , const
//  Kokkos::Experimental::OffsetView<T,P...> & src
//                              , std::string const& name = ""
//                                  , typename
//                                  std::enable_if<!Impl::MirrorViewType<Space,T,P
//                                  ...>::is_same_memspace>::type* = 0 ) {
//    using Mirror = typename
//    Kokkos::Experimental::Impl::MirrorViewType<Space,T,P ...>::view_type;
//    std::string label = name.empty() ? src.label() : name;
//    auto mirror = Mirror(ViewAllocateWithoutInitializing(label), src.layout(),
//                         { src.begin(0), src.begin(1), src.begin(2),
//                         src.begin(3), src.begin(4),
//                             src.begin(5), src.begin(6), src.begin(7) });
//    deep_copy(mirror, src);
//    return mirror;
//  }

}  // namespace Experimental
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* KOKKOS_OFFSETVIEW_HPP_ */
