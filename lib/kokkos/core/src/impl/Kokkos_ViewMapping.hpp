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

#ifndef KOKKOS_EXPERIMENTAL_VIEW_MAPPING_HPP
#define KOKKOS_EXPERIMENTAL_VIEW_MAPPING_HPP

#include <cstring>
#include <type_traits>
#include <initializer_list>

#include <Kokkos_Core_fwd.hpp>
#include <Kokkos_DetectionIdiom.hpp>
#include <Kokkos_Pair.hpp>
#include <Kokkos_Layout.hpp>
#include <Kokkos_Extents.hpp>
#include <impl/Kokkos_Error.hpp>
#include <impl/Kokkos_Traits.hpp>
#include <impl/Kokkos_ViewTracker.hpp>
#include <impl/Kokkos_ViewCtor.hpp>
#include <impl/Kokkos_Atomic_View.hpp>
#include <impl/Kokkos_Tools.hpp>
#include <impl/Kokkos_StringManipulation.hpp>
#include <impl/Kokkos_ZeroMemset_fwd.hpp>
#include <impl/Kokkos_ViewDataAnalysis.hpp>
#include <View/Kokkos_ViewAlloc.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {

struct ALL_t {
  KOKKOS_INLINE_FUNCTION
  constexpr const ALL_t& operator()() const { return *this; }

  KOKKOS_INLINE_FUNCTION
  constexpr bool operator==(const ALL_t&) const { return true; }
};

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
namespace Impl {
// TODO This alias declaration forces us to fully qualify ALL_t inside the
// Kokkos::Impl namespace to avoid deprecation warnings. Replace the
// fully-qualified name when we remove Kokkos::Impl::ALL_t.
using ALL_t KOKKOS_DEPRECATED_WITH_COMMENT("Use Kokkos::ALL_t instead!") =
    Kokkos::ALL_t;
}  // namespace Impl
#endif
}  // namespace Kokkos

namespace Kokkos {
namespace Impl {

template <class T>
struct is_integral_extent_type {
  enum : bool { value = std::is_same<T, Kokkos::ALL_t>::value ? 1 : 0 };
};

template <class iType>
struct is_integral_extent_type<std::pair<iType, iType>> {
  enum : bool { value = std::is_integral<iType>::value ? 1 : 0 };
};

template <class iType>
struct is_integral_extent_type<Kokkos::pair<iType, iType>> {
  enum : bool { value = std::is_integral<iType>::value ? 1 : 0 };
};

// Assuming '2 == initializer_list<iType>::size()'
template <class iType>
struct is_integral_extent_type<std::initializer_list<iType>> {
  enum : bool { value = std::is_integral<iType>::value ? 1 : 0 };
};

template <unsigned I, class... Args>
struct is_integral_extent {
  // get_type is void when sizeof...(Args) <= I
  using type = std::remove_cv_t<std::remove_reference_t<
      typename Kokkos::Impl::get_type<I, Args...>::type>>;

  enum : bool { value = is_integral_extent_type<type>::value };

  static_assert(value || std::is_integral<type>::value ||
                    std::is_void<type>::value,
                "subview argument must be either integral or integral extent");
};

// Rules for subview arguments and layouts matching

template <class LayoutDest, class LayoutSrc, int RankDest, int RankSrc,
          int CurrentArg, class... SubViewArgs>
struct SubviewLegalArgsCompileTime;

// Rules which allow LayoutLeft to LayoutLeft assignment

template <int RankDest, int RankSrc, int CurrentArg, class Arg,
          class... SubViewArgs>
struct SubviewLegalArgsCompileTime<Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                   RankDest, RankSrc, CurrentArg, Arg,
                                   SubViewArgs...> {
  enum {
    value = (((CurrentArg == RankDest - 1) &&
              (Kokkos::Impl::is_integral_extent_type<Arg>::value)) ||
             ((CurrentArg >= RankDest) && (std::is_integral<Arg>::value)) ||
             ((CurrentArg < RankDest) &&
              (std::is_same<Arg, Kokkos::ALL_t>::value)) ||
             ((CurrentArg == 0) &&
              (Kokkos::Impl::is_integral_extent_type<Arg>::value))) &&
            (SubviewLegalArgsCompileTime<Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                         RankDest, RankSrc, CurrentArg + 1,
                                         SubViewArgs...>::value)
  };
};

template <int RankDest, int RankSrc, int CurrentArg, class Arg>
struct SubviewLegalArgsCompileTime<Kokkos::LayoutLeft, Kokkos::LayoutLeft,
                                   RankDest, RankSrc, CurrentArg, Arg> {
  enum {
    value = ((CurrentArg == RankDest - 1) || (std::is_integral<Arg>::value)) &&
            (CurrentArg == RankSrc - 1)
  };
};

// Rules which allow LayoutRight to LayoutRight assignment

template <int RankDest, int RankSrc, int CurrentArg, class Arg,
          class... SubViewArgs>
struct SubviewLegalArgsCompileTime<Kokkos::LayoutRight, Kokkos::LayoutRight,
                                   RankDest, RankSrc, CurrentArg, Arg,
                                   SubViewArgs...> {
  enum {
    value = (((CurrentArg == RankSrc - RankDest) &&
              (Kokkos::Impl::is_integral_extent_type<Arg>::value)) ||
             ((CurrentArg < RankSrc - RankDest) &&
              (std::is_integral<Arg>::value)) ||
             ((CurrentArg >= RankSrc - RankDest) &&
              (std::is_same<Arg, Kokkos::ALL_t>::value))) &&
            (SubviewLegalArgsCompileTime<Kokkos::LayoutRight,
                                         Kokkos::LayoutRight, RankDest, RankSrc,
                                         CurrentArg + 1, SubViewArgs...>::value)
  };
};

template <int RankDest, int RankSrc, int CurrentArg, class Arg>
struct SubviewLegalArgsCompileTime<Kokkos::LayoutRight, Kokkos::LayoutRight,
                                   RankDest, RankSrc, CurrentArg, Arg> {
  enum {
    value = ((CurrentArg == RankSrc - 1) &&
             (std::is_same<Arg, Kokkos::ALL_t>::value))
  };
};

// Rules which allow assignment to LayoutStride

template <int RankDest, int RankSrc, int CurrentArg, class... SubViewArgs>
struct SubviewLegalArgsCompileTime<Kokkos::LayoutStride, Kokkos::LayoutLeft,
                                   RankDest, RankSrc, CurrentArg,
                                   SubViewArgs...> {
  enum : bool { value = true };
};

template <int RankDest, int RankSrc, int CurrentArg, class... SubViewArgs>
struct SubviewLegalArgsCompileTime<Kokkos::LayoutStride, Kokkos::LayoutRight,
                                   RankDest, RankSrc, CurrentArg,
                                   SubViewArgs...> {
  enum : bool { value = true };
};

template <int RankDest, int RankSrc, int CurrentArg, class... SubViewArgs>
struct SubviewLegalArgsCompileTime<Kokkos::LayoutStride, Kokkos::LayoutStride,
                                   RankDest, RankSrc, CurrentArg,
                                   SubViewArgs...> {
  enum : bool { value = true };
};

template <unsigned DomainRank, unsigned RangeRank>
struct SubviewExtents {
 private:
  // Cannot declare zero-length arrays
  // '+' is used to silence GCC 7.2.0 -Wduplicated-branches warning when
  // RangeRank=1
  enum { InternalRangeRank = RangeRank ? RangeRank : +1u };

  size_t m_begin[DomainRank];
  size_t m_length[InternalRangeRank];
  unsigned m_index[InternalRangeRank];

  template <size_t... DimArgs>
  KOKKOS_FORCEINLINE_FUNCTION bool set(unsigned, unsigned,
                                       const ViewDimension<DimArgs...>&) {
    return true;
  }

  template <class T, size_t... DimArgs, class... Args>
  KOKKOS_FORCEINLINE_FUNCTION bool set(unsigned domain_rank,
                                       unsigned range_rank,
                                       const ViewDimension<DimArgs...>& dim,
                                       const T& val, Args... args) {
    const size_t v = static_cast<size_t>(val);

    m_begin[domain_rank] = v;

    return set(domain_rank + 1, range_rank, dim, args...)
#if defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK)
           && (v < dim.extent(domain_rank))
#endif
        ;
  }

  // ALL_t
  template <size_t... DimArgs, class... Args>
  KOKKOS_FORCEINLINE_FUNCTION bool set(unsigned domain_rank,
                                       unsigned range_rank,
                                       const ViewDimension<DimArgs...>& dim,
                                       Kokkos::ALL_t, Args... args) {
    m_begin[domain_rank] = 0;
    m_length[range_rank] = dim.extent(domain_rank);
    m_index[range_rank]  = domain_rank;

    return set(domain_rank + 1, range_rank + 1, dim, args...);
  }

  // std::pair range
  template <class T, size_t... DimArgs, class... Args>
  KOKKOS_FORCEINLINE_FUNCTION bool set(unsigned domain_rank,
                                       unsigned range_rank,
                                       const ViewDimension<DimArgs...>& dim,
                                       const std::pair<T, T>& val,
                                       Args... args) {
    const size_t b = static_cast<size_t>(val.first);
    const size_t e = static_cast<size_t>(val.second);

    m_begin[domain_rank] = b;
    m_length[range_rank] = e - b;
    m_index[range_rank]  = domain_rank;

    return set(domain_rank + 1, range_rank + 1, dim, args...)
#if defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK)
           && (e <= b + dim.extent(domain_rank))
#endif
        ;
  }

  // Kokkos::pair range
  template <class T, size_t... DimArgs, class... Args>
  KOKKOS_FORCEINLINE_FUNCTION bool set(unsigned domain_rank,
                                       unsigned range_rank,
                                       const ViewDimension<DimArgs...>& dim,
                                       const Kokkos::pair<T, T>& val,
                                       Args... args) {
    const size_t b = static_cast<size_t>(val.first);
    const size_t e = static_cast<size_t>(val.second);

    m_begin[domain_rank] = b;
    m_length[range_rank] = e - b;
    m_index[range_rank]  = domain_rank;

    return set(domain_rank + 1, range_rank + 1, dim, args...)
#if defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK)
           && (e <= b + dim.extent(domain_rank))
#endif
        ;
  }

  // { begin , end } range
  template <class T, size_t... DimArgs, class... Args>
  KOKKOS_FORCEINLINE_FUNCTION bool set(unsigned domain_rank,
                                       unsigned range_rank,
                                       const ViewDimension<DimArgs...>& dim,
                                       const std::initializer_list<T>& val,
                                       Args... args) {
    const size_t b = static_cast<size_t>(val.begin()[0]);
    const size_t e = static_cast<size_t>(val.begin()[1]);

    m_begin[domain_rank] = b;
    m_length[range_rank] = e - b;
    m_index[range_rank]  = domain_rank;

    return set(domain_rank + 1, range_rank + 1, dim, args...)
#if defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK)
           && (val.size() == 2) && (e <= b + dim.extent(domain_rank))
#endif
        ;
  }

  //------------------------------

#if defined(KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK)

  template <size_t... DimArgs>
  void error(char*, int, unsigned, unsigned,
             const ViewDimension<DimArgs...>&) const {}

  template <class T, size_t... DimArgs, class... Args>
  void error(char* buf, int buf_len, unsigned domain_rank, unsigned range_rank,
             const ViewDimension<DimArgs...>& dim, const T& val,
             Args... args) const {
    const int n = std::min(
        buf_len,
        snprintf(buf, buf_len, " %lu < %lu %c", static_cast<unsigned long>(val),
                 static_cast<unsigned long>(dim.extent(domain_rank)),
                 int(sizeof...(Args) ? ',' : ')')));

    error(buf + n, buf_len - n, domain_rank + 1, range_rank, dim, args...);
  }

  // std::pair range
  template <size_t... DimArgs, class... Args>
  void error(char* buf, int buf_len, unsigned domain_rank, unsigned range_rank,
             const ViewDimension<DimArgs...>& dim, Kokkos::ALL_t,
             Args... args) const {
    const int n = std::min(buf_len, snprintf(buf, buf_len, " Kokkos::ALL %c",
                                             int(sizeof...(Args) ? ',' : ')')));

    error(buf + n, buf_len - n, domain_rank + 1, range_rank + 1, dim, args...);
  }

  // std::pair range
  template <class T, size_t... DimArgs, class... Args>
  void error(char* buf, int buf_len, unsigned domain_rank, unsigned range_rank,
             const ViewDimension<DimArgs...>& dim, const std::pair<T, T>& val,
             Args... args) const {
    // d <= e - b
    const int n = std::min(
        buf_len, snprintf(buf, buf_len, " %lu <= %lu - %lu %c",
                          static_cast<unsigned long>(dim.extent(domain_rank)),
                          static_cast<unsigned long>(val.second),
                          static_cast<unsigned long>(val.first),
                          int(sizeof...(Args) ? ',' : ')')));

    error(buf + n, buf_len - n, domain_rank + 1, range_rank + 1, dim, args...);
  }

  // Kokkos::pair range
  template <class T, size_t... DimArgs, class... Args>
  void error(char* buf, int buf_len, unsigned domain_rank, unsigned range_rank,
             const ViewDimension<DimArgs...>& dim,
             const Kokkos::pair<T, T>& val, Args... args) const {
    // d <= e - b
    const int n = std::min(
        buf_len, snprintf(buf, buf_len, " %lu <= %lu - %lu %c",
                          static_cast<unsigned long>(dim.extent(domain_rank)),
                          static_cast<unsigned long>(val.second),
                          static_cast<unsigned long>(val.first),
                          int(sizeof...(Args) ? ',' : ')')));

    error(buf + n, buf_len - n, domain_rank + 1, range_rank + 1, dim, args...);
  }

  // { begin , end } range
  template <class T, size_t... DimArgs, class... Args>
  void error(char* buf, int buf_len, unsigned domain_rank, unsigned range_rank,
             const ViewDimension<DimArgs...>& dim,
             const std::initializer_list<T>& val, Args... args) const {
    // d <= e - b
    int n = 0;
    if (val.size() == 2) {
      n = std::min(buf_len,
                   snprintf(buf, buf_len, " %lu <= %lu - %lu %c",
                            static_cast<unsigned long>(dim.extent(domain_rank)),
                            static_cast<unsigned long>(val.begin()[0]),
                            static_cast<unsigned long>(val.begin()[1]),
                            int(sizeof...(Args) ? ',' : ')')));
    } else {
      n = std::min(buf_len, snprintf(buf, buf_len, " { ... }.size() == %u %c",
                                     unsigned(val.size()),
                                     int(sizeof...(Args) ? ',' : ')')));
    }

    error(buf + n, buf_len - n, domain_rank + 1, range_rank + 1, dim, args...);
  }

  template <size_t... DimArgs, class... Args>
  KOKKOS_FORCEINLINE_FUNCTION void error(const ViewDimension<DimArgs...>& dim,
                                         Args... args) const {
    KOKKOS_IF_ON_HOST(
        (enum {LEN = 1024}; char buffer[LEN];

         const int n = snprintf(buffer, LEN, "Kokkos::subview bounds error (");
         error(buffer + n, LEN - n, 0, 0, dim, args...);

         Kokkos::Impl::throw_runtime_exception(std::string(buffer));))

    KOKKOS_IF_ON_DEVICE(((void)dim;
                         Kokkos::abort("Kokkos::subview bounds error");
                         [](Args...) {}(args...);))
  }

#else

  template <size_t... DimArgs, class... Args>
  KOKKOS_FORCEINLINE_FUNCTION void error(const ViewDimension<DimArgs...>&,
                                         Args...) const {}

#endif

 public:
  template <size_t... DimArgs, class... Args>
  KOKKOS_INLINE_FUNCTION SubviewExtents(const ViewDimension<DimArgs...>& dim,
                                        Args... args) {
    static_assert(DomainRank == sizeof...(DimArgs));
    static_assert(DomainRank == sizeof...(Args));

    // Verifies that all arguments, up to 8, are integral types,
    // integral extents, or don't exist.
    static_assert(RangeRank ==
                  unsigned(is_integral_extent<0, Args...>::value) +
                      unsigned(is_integral_extent<1, Args...>::value) +
                      unsigned(is_integral_extent<2, Args...>::value) +
                      unsigned(is_integral_extent<3, Args...>::value) +
                      unsigned(is_integral_extent<4, Args...>::value) +
                      unsigned(is_integral_extent<5, Args...>::value) +
                      unsigned(is_integral_extent<6, Args...>::value) +
                      unsigned(is_integral_extent<7, Args...>::value));

    if (RangeRank == 0) {
      m_length[0] = 0;
      m_index[0]  = ~0u;
    }

    if (!set(0, 0, dim, args...)) error(dim, args...);
  }

  template <typename iType>
  KOKKOS_FORCEINLINE_FUNCTION constexpr size_t domain_offset(
      const iType i) const {
    return unsigned(i) < DomainRank ? m_begin[i] : 0;
  }

  template <typename iType>
  KOKKOS_FORCEINLINE_FUNCTION constexpr size_t range_extent(
      const iType i) const {
    return unsigned(i) < InternalRangeRank ? m_length[i] : 0;
  }

  template <typename iType>
  KOKKOS_FORCEINLINE_FUNCTION constexpr unsigned range_index(
      const iType i) const {
    return unsigned(i) < InternalRangeRank ? m_index[i] : ~0u;
  }
};

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {
//----------------------------------------------------------------------------
// LayoutLeft AND ( 1 >= rank OR 0 == rank_dynamic ) : no padding / striding
template <class Dimension>
struct ViewOffset<
    Dimension, Kokkos::LayoutLeft,
    std::enable_if_t<(1 >= Dimension::rank || 0 == Dimension::rank_dynamic)>> {
  using is_mapping_plugin = std::true_type;
  using is_regular        = std::true_type;

  using size_type      = size_t;
  using dimension_type = Dimension;
  using array_layout   = Kokkos::LayoutLeft;

  dimension_type m_dim;

  //----------------------------------------

  // rank 1
  template <typename I0>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(I0 const& i0) const {
    return i0;
  }

  // rank 2
  template <typename I0, typename I1>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(I0 const& i0,
                                                        I1 const& i1) const {
    return i0 + m_dim.N0 * i1;
  }

  // rank 3
  template <typename I0, typename I1, typename I2>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(I0 const& i0,
                                                        I1 const& i1,
                                                        I2 const& i2) const {
    return i0 + m_dim.N0 * (i1 + m_dim.N1 * i2);
  }

  // rank 4
  template <typename I0, typename I1, typename I2, typename I3>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(I0 const& i0,
                                                        I1 const& i1,
                                                        I2 const& i2,
                                                        I3 const& i3) const {
    return i0 + m_dim.N0 * (i1 + m_dim.N1 * (i2 + m_dim.N2 * i3));
  }

  // rank 5
  template <typename I0, typename I1, typename I2, typename I3, typename I4>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(I0 const& i0,
                                                        I1 const& i1,
                                                        I2 const& i2,
                                                        I3 const& i3,
                                                        I4 const& i4) const {
    return i0 +
           m_dim.N0 * (i1 + m_dim.N1 * (i2 + m_dim.N2 * (i3 + m_dim.N3 * i4)));
  }

  // rank 6
  template <typename I0, typename I1, typename I2, typename I3, typename I4,
            typename I5>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(
      I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4,
      I5 const& i5) const {
    return i0 +
           m_dim.N0 *
               (i1 +
                m_dim.N1 *
                    (i2 + m_dim.N2 * (i3 + m_dim.N3 * (i4 + m_dim.N4 * i5))));
  }

  // rank 7
  template <typename I0, typename I1, typename I2, typename I3, typename I4,
            typename I5, typename I6>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(
      I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4,
      I5 const& i5, I6 const& i6) const {
    return i0 +
           m_dim.N0 *
               (i1 + m_dim.N1 *
                         (i2 + m_dim.N2 *
                                   (i3 + m_dim.N3 *
                                             (i4 + m_dim.N4 *
                                                       (i5 + m_dim.N5 * i6)))));
  }

  // rank 8
  template <typename I0, typename I1, typename I2, typename I3, typename I4,
            typename I5, typename I6, typename I7>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(
      I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4,
      I5 const& i5, I6 const& i6, I7 const& i7) const {
    return i0 +
           m_dim.N0 *
               (i1 +
                m_dim.N1 *
                    (i2 + m_dim.N2 *
                              (i3 + m_dim.N3 *
                                        (i4 + m_dim.N4 *
                                                  (i5 + m_dim.N5 *
                                                            (i6 + m_dim.N6 *
                                                                      i7))))));
  }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  constexpr array_layout layout() const {
    constexpr auto r = dimension_type::rank;
    return array_layout((r > 0 ? m_dim.N0 : KOKKOS_INVALID_INDEX),
                        (r > 1 ? m_dim.N1 : KOKKOS_INVALID_INDEX),
                        (r > 2 ? m_dim.N2 : KOKKOS_INVALID_INDEX),
                        (r > 3 ? m_dim.N3 : KOKKOS_INVALID_INDEX),
                        (r > 4 ? m_dim.N4 : KOKKOS_INVALID_INDEX),
                        (r > 5 ? m_dim.N5 : KOKKOS_INVALID_INDEX),
                        (r > 6 ? m_dim.N6 : KOKKOS_INVALID_INDEX),
                        (r > 7 ? m_dim.N7 : KOKKOS_INVALID_INDEX));
  }

  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_0() const {
    return m_dim.N0;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_1() const {
    return m_dim.N1;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_2() const {
    return m_dim.N2;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_3() const {
    return m_dim.N3;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_4() const {
    return m_dim.N4;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_5() const {
    return m_dim.N5;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_6() const {
    return m_dim.N6;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_7() const {
    return m_dim.N7;
  }

  /* Cardinality of the domain index space */
  KOKKOS_INLINE_FUNCTION
  constexpr size_type size() const {
    return size_type(m_dim.N0) * m_dim.N1 * m_dim.N2 * m_dim.N3 * m_dim.N4 *
           m_dim.N5 * m_dim.N6 * m_dim.N7;
  }

  /* Span of the range space */
  KOKKOS_INLINE_FUNCTION
  constexpr size_type span() const {
    return size_type(m_dim.N0) * m_dim.N1 * m_dim.N2 * m_dim.N3 * m_dim.N4 *
           m_dim.N5 * m_dim.N6 * m_dim.N7;
  }

  KOKKOS_INLINE_FUNCTION constexpr bool span_is_contiguous() const {
    return true;
  }

  /* Strides of dimensions */
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_0() const { return 1; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_1() const {
    return m_dim.N0;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_2() const {
    return size_type(m_dim.N0) * m_dim.N1;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_3() const {
    return size_type(m_dim.N0) * m_dim.N1 * m_dim.N2;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_4() const {
    return size_type(m_dim.N0) * m_dim.N1 * m_dim.N2 * m_dim.N3;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_5() const {
    return size_type(m_dim.N0) * m_dim.N1 * m_dim.N2 * m_dim.N3 * m_dim.N4;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_6() const {
    return size_type(m_dim.N0) * m_dim.N1 * m_dim.N2 * m_dim.N3 * m_dim.N4 *
           m_dim.N5;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_7() const {
    return size_type(m_dim.N0) * m_dim.N1 * m_dim.N2 * m_dim.N3 * m_dim.N4 *
           m_dim.N5 * m_dim.N6;
  }

  // Fill the target unbounded array s with the stride.
  // This method differs from stride() in that it does not write the total
  // length to the last index of the array. Preconditions: s must be an array of
  // dimension_type::rank elements
  // FIXME: The version of clang-format in CI fails from maybe_unused
  // clang-format off
  template <typename iType>
  KOKKOS_INLINE_FUNCTION iType
  stride_fill([[maybe_unused]] iType* const s) const {
    iType n = 1;
    if constexpr (0 < dimension_type::rank) {
      s[0] = n;
      n *= m_dim.N0;
    }
    if constexpr (1 < dimension_type::rank) {
      s[1] = n;
      n *= m_dim.N1;
    }
    if constexpr (2 < dimension_type::rank) {
      s[2] = n;
      n *= m_dim.N2;
    }
    if constexpr (3 < dimension_type::rank) {
      s[3] = n;
      n *= m_dim.N3;
    }
    if constexpr (4 < dimension_type::rank) {
      s[4] = n;
      n *= m_dim.N4;
    }
    if constexpr (5 < dimension_type::rank) {
      s[5] = n;
      n *= m_dim.N5;
    }
    if constexpr (6 < dimension_type::rank) {
      s[6] = n;
      n *= m_dim.N6;
    }
    if constexpr (7 < dimension_type::rank) {
      s[7] = n;
      n *= m_dim.N7;
    }
    return n;
  }
  // clang-format on

  // Fill the target unbounded array s with the stride and the total spanned
  // size. This method differs from stride_fill() in that it writes the total
  // spanned size to the last index of the array. Preconditions: s must be an
  // array of dimension_type::rank + 1 elements Stride with [ rank ] value is
  // the total length
  template <typename iType>
  KOKKOS_INLINE_FUNCTION void stride(iType* const s) const {
    s[dimension_type::rank] = stride_fill(s);
  }

  //----------------------------------------

  // MSVC (16.5.5) + CUDA (10.2) did not generate the defaulted functions
  // correct and errors out during compilation. Same for the other places where
  // I changed this.
#ifdef KOKKOS_IMPL_WINDOWS_CUDA
  KOKKOS_FUNCTION ViewOffset() : m_dim(dimension_type()) {}
  KOKKOS_FUNCTION ViewOffset(const ViewOffset& src) { m_dim = src.m_dim; }
  KOKKOS_FUNCTION ViewOffset& operator=(const ViewOffset& src) {
    m_dim = src.m_dim;
    return *this;
  }
#else
  ViewOffset()                  = default;
  ViewOffset(const ViewOffset&) = default;
  ViewOffset& operator=(const ViewOffset&) = default;
#endif

  template <unsigned TrivialScalarSize>
  KOKKOS_INLINE_FUNCTION constexpr ViewOffset(
      std::integral_constant<unsigned, TrivialScalarSize> const&,
      Kokkos::LayoutLeft const& arg_layout)
      : m_dim(arg_layout.dimension[0], 0, 0, 0, 0, 0, 0, 0) {}

  template <class DimRHS>
  KOKKOS_INLINE_FUNCTION constexpr ViewOffset(
      const ViewOffset<DimRHS, Kokkos::LayoutLeft, void>& rhs)
      : m_dim(rhs.m_dim.N0, rhs.m_dim.N1, rhs.m_dim.N2, rhs.m_dim.N3,
              rhs.m_dim.N4, rhs.m_dim.N5, rhs.m_dim.N6, rhs.m_dim.N7) {
    static_assert(int(DimRHS::rank) == int(dimension_type::rank),
                  "ViewOffset assignment requires equal rank");
    // Also requires equal static dimensions ...
  }

  template <class DimRHS>
  KOKKOS_INLINE_FUNCTION constexpr ViewOffset(
      const ViewOffset<DimRHS, Kokkos::LayoutRight, void>& rhs)
      : m_dim(rhs.m_dim.N0, 0, 0, 0, 0, 0, 0, 0) {
    static_assert(((DimRHS::rank == 0 && dimension_type::rank == 0) ||
                   (DimRHS::rank == 1 && dimension_type::rank == 1)),
                  "ViewOffset LayoutLeft and LayoutRight are only compatible "
                  "when rank <= 1");
  }

  template <class DimRHS>
  KOKKOS_INLINE_FUNCTION ViewOffset(
      const ViewOffset<DimRHS, Kokkos::LayoutStride, void>& rhs)
      : m_dim(rhs.m_dim.N0, 0, 0, 0, 0, 0, 0, 0) {
    if (rhs.m_stride.S0 != 1) {
      Kokkos::abort(
          "Kokkos::Impl::ViewOffset assignment of LayoutLeft from LayoutStride "
          " requires stride == 1");
    }
  }

  //----------------------------------------
  // Subview construction

  template <class DimRHS>
  KOKKOS_INLINE_FUNCTION constexpr ViewOffset(
      const ViewOffset<DimRHS, Kokkos::LayoutLeft, void>&,
      const SubviewExtents<DimRHS::rank, dimension_type::rank>& sub)
      : m_dim(sub.range_extent(0), 0, 0, 0, 0, 0, 0, 0) {
    static_assert((0 == dimension_type::rank_dynamic) ||
                      (1 == dimension_type::rank &&
                       1 == dimension_type::rank_dynamic && 1 <= DimRHS::rank),
                  "ViewOffset subview construction requires compatible rank");
  }
};

//----------------------------------------------------------------------------
// LayoutLeft AND ( 1 < rank AND 0 < rank_dynamic ) : has padding / striding
template <class Dimension>
struct ViewOffset<
    Dimension, Kokkos::LayoutLeft,
    std::enable_if_t<(1 < Dimension::rank && 0 < Dimension::rank_dynamic)>> {
  using is_mapping_plugin = std::true_type;
  using is_regular        = std::true_type;

  using size_type      = size_t;
  using dimension_type = Dimension;
  using array_layout   = Kokkos::LayoutLeft;

  dimension_type m_dim;
  size_type m_stride;

  //----------------------------------------

  // rank 1
  template <typename I0>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(I0 const& i0) const {
    return i0;
  }

  // rank 2
  template <typename I0, typename I1>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(I0 const& i0,
                                                        I1 const& i1) const {
    return i0 + m_stride * i1;
  }

  // rank 3
  template <typename I0, typename I1, typename I2>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(I0 const& i0,
                                                        I1 const& i1,
                                                        I2 const& i2) const {
    return i0 + m_stride * (i1 + m_dim.N1 * i2);
  }

  // rank 4
  template <typename I0, typename I1, typename I2, typename I3>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(I0 const& i0,
                                                        I1 const& i1,
                                                        I2 const& i2,
                                                        I3 const& i3) const {
    return i0 + m_stride * (i1 + m_dim.N1 * (i2 + m_dim.N2 * i3));
  }

  // rank 5
  template <typename I0, typename I1, typename I2, typename I3, typename I4>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(I0 const& i0,
                                                        I1 const& i1,
                                                        I2 const& i2,
                                                        I3 const& i3,
                                                        I4 const& i4) const {
    return i0 +
           m_stride * (i1 + m_dim.N1 * (i2 + m_dim.N2 * (i3 + m_dim.N3 * i4)));
  }

  // rank 6
  template <typename I0, typename I1, typename I2, typename I3, typename I4,
            typename I5>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(
      I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4,
      I5 const& i5) const {
    return i0 +
           m_stride *
               (i1 +
                m_dim.N1 *
                    (i2 + m_dim.N2 * (i3 + m_dim.N3 * (i4 + m_dim.N4 * i5))));
  }

  // rank 7
  template <typename I0, typename I1, typename I2, typename I3, typename I4,
            typename I5, typename I6>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(
      I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4,
      I5 const& i5, I6 const& i6) const {
    return i0 +
           m_stride *
               (i1 + m_dim.N1 *
                         (i2 + m_dim.N2 *
                                   (i3 + m_dim.N3 *
                                             (i4 + m_dim.N4 *
                                                       (i5 + m_dim.N5 * i6)))));
  }

  // rank 8
  template <typename I0, typename I1, typename I2, typename I3, typename I4,
            typename I5, typename I6, typename I7>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(
      I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4,
      I5 const& i5, I6 const& i6, I7 const& i7) const {
    return i0 +
           m_stride *
               (i1 +
                m_dim.N1 *
                    (i2 + m_dim.N2 *
                              (i3 + m_dim.N3 *
                                        (i4 + m_dim.N4 *
                                                  (i5 + m_dim.N5 *
                                                            (i6 + m_dim.N6 *
                                                                      i7))))));
  }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  constexpr array_layout layout() const {
    constexpr auto r = dimension_type::rank;
    return array_layout((r > 0 ? m_dim.N0 : KOKKOS_INVALID_INDEX),
                        (r > 1 ? m_dim.N1 : KOKKOS_INVALID_INDEX),
                        (r > 2 ? m_dim.N2 : KOKKOS_INVALID_INDEX),
                        (r > 3 ? m_dim.N3 : KOKKOS_INVALID_INDEX),
                        (r > 4 ? m_dim.N4 : KOKKOS_INVALID_INDEX),
                        (r > 5 ? m_dim.N5 : KOKKOS_INVALID_INDEX),
                        (r > 6 ? m_dim.N6 : KOKKOS_INVALID_INDEX),
                        (r > 7 ? m_dim.N7 : KOKKOS_INVALID_INDEX));
  }

  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_0() const {
    return m_dim.N0;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_1() const {
    return m_dim.N1;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_2() const {
    return m_dim.N2;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_3() const {
    return m_dim.N3;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_4() const {
    return m_dim.N4;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_5() const {
    return m_dim.N5;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_6() const {
    return m_dim.N6;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_7() const {
    return m_dim.N7;
  }

  /* Cardinality of the domain index space */
  KOKKOS_INLINE_FUNCTION
  constexpr size_type size() const {
    return size_type(m_dim.N0) * m_dim.N1 * m_dim.N2 * m_dim.N3 * m_dim.N4 *
           m_dim.N5 * m_dim.N6 * m_dim.N7;
  }

  /* Span of the range space */
  KOKKOS_INLINE_FUNCTION
  constexpr size_type span() const {
    return (m_dim.N0 > size_type(0) ? m_stride : size_type(0)) * m_dim.N1 *
           m_dim.N2 * m_dim.N3 * m_dim.N4 * m_dim.N5 * m_dim.N6 * m_dim.N7;
  }

  KOKKOS_INLINE_FUNCTION constexpr bool span_is_contiguous() const {
    return m_stride == m_dim.N0;
  }

  /* Strides of dimensions */
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_0() const { return 1; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_1() const {
    return m_stride;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_2() const {
    return m_stride * m_dim.N1;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_3() const {
    return m_stride * m_dim.N1 * m_dim.N2;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_4() const {
    return m_stride * m_dim.N1 * m_dim.N2 * m_dim.N3;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_5() const {
    return m_stride * m_dim.N1 * m_dim.N2 * m_dim.N3 * m_dim.N4;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_6() const {
    return m_stride * m_dim.N1 * m_dim.N2 * m_dim.N3 * m_dim.N4 * m_dim.N5;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_7() const {
    return m_stride * m_dim.N1 * m_dim.N2 * m_dim.N3 * m_dim.N4 * m_dim.N5 *
           m_dim.N6;
  }

  // Fill the target unbounded array s with the stride.
  // This method differs from stride() in that it does not write the total
  // length to the last index of the array. Preconditions: s must be an array of
  // dimension_type::rank elements
  // The version of clang-format in CI fails from maybe_unused
  // clang-format off
  template <typename iType>
  KOKKOS_INLINE_FUNCTION iType
  stride_fill([[maybe_unused]] iType* const s) const {
    iType n = 1;
    if constexpr (0 < dimension_type::rank) {
      s[0] = n;
      n *= m_stride;
    }
    if constexpr (1 < dimension_type::rank) {
      s[1] = n;
      n *= m_dim.N1;
    }
    if constexpr (2 < dimension_type::rank) {
      s[2] = n;
      n *= m_dim.N2;
    }
    if constexpr (3 < dimension_type::rank) {
      s[3] = n;
      n *= m_dim.N3;
    }
    if constexpr (4 < dimension_type::rank) {
      s[4] = n;
      n *= m_dim.N4;
    }
    if constexpr (5 < dimension_type::rank) {
      s[5] = n;
      n *= m_dim.N5;
    }
    if constexpr (6 < dimension_type::rank) {
      s[6] = n;
      n *= m_dim.N6;
    }
    if constexpr (7 < dimension_type::rank) {
      s[7] = n;
      n *= m_dim.N7;
    }
    return n;
  }
  // clang-format on

  // Fill the target unbounded array s with the stride and the total spanned
  // size. This method differs from stride_fill() in that it writes the total
  // spanned size to the last index of the array. Preconditions: s must be an
  // array of dimension_type::rank + 1 elements
  template <typename iType>
  KOKKOS_INLINE_FUNCTION void stride(iType* const s) const {
    s[dimension_type::rank] = stride_fill(s);
  }

  //----------------------------------------

 private:
  template <unsigned TrivialScalarSize>
  struct Padding {
    enum {
      div = TrivialScalarSize == 0
                ? 0
                : Kokkos::Impl::MEMORY_ALIGNMENT /
                      (TrivialScalarSize ? TrivialScalarSize : 1)
    };
    enum {
      mod = TrivialScalarSize == 0
                ? 0
                : Kokkos::Impl::MEMORY_ALIGNMENT %
                      (TrivialScalarSize ? TrivialScalarSize : 1)
    };

    // If memory alignment is a multiple of the trivial scalar size then attempt
    // to align.
    enum { align = 0 != TrivialScalarSize && 0 == mod ? div : 0 };
    enum {
      div_ok = (div != 0) ? div : 1
    };  // To valid modulo zero in constexpr

    KOKKOS_INLINE_FUNCTION
    static constexpr size_t stride(size_t const N) {
      return ((align != 0) &&
              ((static_cast<int>(Kokkos::Impl::MEMORY_ALIGNMENT_THRESHOLD) *
                static_cast<int>(align)) < N) &&
              ((N % div_ok) != 0))
                 ? N + align - (N % div_ok)
                 : N;
    }
  };

 public:
  // MSVC (16.5.5) + CUDA (10.2) did not generate the defaulted functions
  // correct and errors out during compilation. Same for the other places where
  // I changed this.
#ifdef KOKKOS_IMPL_WINDOWS_CUDA
  KOKKOS_FUNCTION ViewOffset() : m_dim(dimension_type()), m_stride(0) {}
  KOKKOS_FUNCTION ViewOffset(const ViewOffset& src) {
    m_dim    = src.m_dim;
    m_stride = src.m_stride;
  }
  KOKKOS_FUNCTION ViewOffset& operator=(const ViewOffset& src) {
    m_dim    = src.m_dim;
    m_stride = src.m_stride;
    return *this;
  }
#else

  ViewOffset()                  = default;
  ViewOffset(const ViewOffset&) = default;
  ViewOffset& operator=(const ViewOffset&) = default;
#endif

  /* Enable padding for trivial scalar types with non-zero trivial scalar size
   */
  template <unsigned TrivialScalarSize>
  KOKKOS_INLINE_FUNCTION constexpr ViewOffset(
      std::integral_constant<unsigned, TrivialScalarSize> const&,
      Kokkos::LayoutLeft const& arg_layout)
      : m_dim(arg_layout.dimension[0], arg_layout.dimension[1],
              arg_layout.dimension[2], arg_layout.dimension[3],
              arg_layout.dimension[4], arg_layout.dimension[5],
              arg_layout.dimension[6], arg_layout.dimension[7]),
        m_stride(Padding<TrivialScalarSize>::stride(arg_layout.dimension[0])) {}

  template <class DimRHS>
  KOKKOS_INLINE_FUNCTION constexpr ViewOffset(
      const ViewOffset<DimRHS, Kokkos::LayoutLeft, void>& rhs)
      : m_dim(rhs.m_dim.N0, rhs.m_dim.N1, rhs.m_dim.N2, rhs.m_dim.N3,
              rhs.m_dim.N4, rhs.m_dim.N5, rhs.m_dim.N6, rhs.m_dim.N7),
        m_stride(rhs.stride_1()) {
    static_assert(int(DimRHS::rank) == int(dimension_type::rank),
                  "ViewOffset assignment requires equal rank");
    // Also requires equal static dimensions ...
  }

  template <class DimRHS>
  KOKKOS_INLINE_FUNCTION ViewOffset(
      const ViewOffset<DimRHS, Kokkos::LayoutStride, void>& rhs)
      : m_dim(rhs.m_dim.N0, rhs.m_dim.N1, rhs.m_dim.N2, rhs.m_dim.N3,
              rhs.m_dim.N4, rhs.m_dim.N5, rhs.m_dim.N6, rhs.m_dim.N7),
        m_stride(rhs.stride_1()) {
    if (rhs.m_stride.S0 != 1) {
      Kokkos::abort(
          "Kokkos::Impl::ViewOffset assignment of LayoutLeft from LayoutStride "
          "requires stride == 1");
    }
  }

  //----------------------------------------
  // Subview construction
  // This subview must be 2 == rank and 2 == rank_dynamic
  // due to only having stride #0.
  // The source dimension #0 must be non-zero for stride-one leading dimension.
  // At most subsequent dimension can be non-zero.

  template <class DimRHS>
  KOKKOS_INLINE_FUNCTION constexpr ViewOffset(
      const ViewOffset<DimRHS, Kokkos::LayoutLeft, void>& rhs,
      const SubviewExtents<DimRHS::rank, dimension_type::rank>& sub)
      : m_dim(sub.range_extent(0), sub.range_extent(1), sub.range_extent(2),
              sub.range_extent(3), sub.range_extent(4), sub.range_extent(5),
              sub.range_extent(6), sub.range_extent(7)),
        m_stride(
            (1 == sub.range_index(1)
                 ? rhs.stride_1()
                 : (2 == sub.range_index(1)
                        ? rhs.stride_2()
                        : (3 == sub.range_index(1)
                               ? rhs.stride_3()
                               : (4 == sub.range_index(1)
                                      ? rhs.stride_4()
                                      : (5 == sub.range_index(1)
                                             ? rhs.stride_5()
                                             : (6 == sub.range_index(1)
                                                    ? rhs.stride_6()
                                                    : (7 == sub.range_index(1)
                                                           ? rhs.stride_7()
                                                           : 0)))))))) {
    // static_assert( ( 2 == dimension_type::rank ) &&
    //               ( 2 == dimension_type::rank_dynamic ) &&
    //               ( 2 <= DimRHS::rank )
    //             , "ViewOffset subview construction requires compatible rank"
    //             );
  }
};

//----------------------------------------------------------------------------
// LayoutRight AND ( 1 >= rank OR 0 == rank_dynamic ) : no padding / striding
template <class Dimension>
struct ViewOffset<
    Dimension, Kokkos::LayoutRight,
    std::enable_if_t<(1 >= Dimension::rank || 0 == Dimension::rank_dynamic)>> {
  using is_mapping_plugin = std::true_type;
  using is_regular        = std::true_type;

  using size_type      = size_t;
  using dimension_type = Dimension;
  using array_layout   = Kokkos::LayoutRight;

  dimension_type m_dim;

  //----------------------------------------

  // rank 1
  template <typename I0>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(I0 const& i0) const {
    return i0;
  }

  // rank 2
  template <typename I0, typename I1>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(I0 const& i0,
                                                        I1 const& i1) const {
    return i1 + m_dim.N1 * i0;
  }

  // rank 3
  template <typename I0, typename I1, typename I2>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(I0 const& i0,
                                                        I1 const& i1,
                                                        I2 const& i2) const {
    return i2 + m_dim.N2 * (i1 + m_dim.N1 * (i0));
  }

  // rank 4
  template <typename I0, typename I1, typename I2, typename I3>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(I0 const& i0,
                                                        I1 const& i1,
                                                        I2 const& i2,
                                                        I3 const& i3) const {
    return i3 + m_dim.N3 * (i2 + m_dim.N2 * (i1 + m_dim.N1 * (i0)));
  }

  // rank 5
  template <typename I0, typename I1, typename I2, typename I3, typename I4>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(I0 const& i0,
                                                        I1 const& i1,
                                                        I2 const& i2,
                                                        I3 const& i3,
                                                        I4 const& i4) const {
    return i4 + m_dim.N4 *
                    (i3 + m_dim.N3 * (i2 + m_dim.N2 * (i1 + m_dim.N1 * (i0))));
  }

  // rank 6
  template <typename I0, typename I1, typename I2, typename I3, typename I4,
            typename I5>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(
      I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4,
      I5 const& i5) const {
    return i5 +
           m_dim.N5 *
               (i4 +
                m_dim.N4 *
                    (i3 + m_dim.N3 * (i2 + m_dim.N2 * (i1 + m_dim.N1 * (i0)))));
  }

  // rank 7
  template <typename I0, typename I1, typename I2, typename I3, typename I4,
            typename I5, typename I6>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(
      I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4,
      I5 const& i5, I6 const& i6) const {
    return i6 +
           m_dim.N6 *
               (i5 +
                m_dim.N5 *
                    (i4 +
                     m_dim.N4 *
                         (i3 + m_dim.N3 *
                                   (i2 + m_dim.N2 * (i1 + m_dim.N1 * (i0))))));
  }

  // rank 8
  template <typename I0, typename I1, typename I2, typename I3, typename I4,
            typename I5, typename I6, typename I7>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(
      I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4,
      I5 const& i5, I6 const& i6, I7 const& i7) const {
    return i7 +
           m_dim.N7 *
               (i6 +
                m_dim.N6 *
                    (i5 +
                     m_dim.N5 *
                         (i4 +
                          m_dim.N4 *
                              (i3 +
                               m_dim.N3 *
                                   (i2 + m_dim.N2 * (i1 + m_dim.N1 * (i0)))))));
  }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  constexpr array_layout layout() const {
    constexpr auto r = dimension_type::rank;
    return array_layout((r > 0 ? m_dim.N0 : KOKKOS_INVALID_INDEX),
                        (r > 1 ? m_dim.N1 : KOKKOS_INVALID_INDEX),
                        (r > 2 ? m_dim.N2 : KOKKOS_INVALID_INDEX),
                        (r > 3 ? m_dim.N3 : KOKKOS_INVALID_INDEX),
                        (r > 4 ? m_dim.N4 : KOKKOS_INVALID_INDEX),
                        (r > 5 ? m_dim.N5 : KOKKOS_INVALID_INDEX),
                        (r > 6 ? m_dim.N6 : KOKKOS_INVALID_INDEX),
                        (r > 7 ? m_dim.N7 : KOKKOS_INVALID_INDEX));
  }

  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_0() const {
    return m_dim.N0;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_1() const {
    return m_dim.N1;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_2() const {
    return m_dim.N2;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_3() const {
    return m_dim.N3;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_4() const {
    return m_dim.N4;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_5() const {
    return m_dim.N5;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_6() const {
    return m_dim.N6;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_7() const {
    return m_dim.N7;
  }

  /* Cardinality of the domain index space */
  KOKKOS_INLINE_FUNCTION
  constexpr size_type size() const {
    return size_type(m_dim.N0) * m_dim.N1 * m_dim.N2 * m_dim.N3 * m_dim.N4 *
           m_dim.N5 * m_dim.N6 * m_dim.N7;
  }

  /* Span of the range space */
  KOKKOS_INLINE_FUNCTION
  constexpr size_type span() const {
    return size_type(m_dim.N0) * m_dim.N1 * m_dim.N2 * m_dim.N3 * m_dim.N4 *
           m_dim.N5 * m_dim.N6 * m_dim.N7;
  }

  KOKKOS_INLINE_FUNCTION constexpr bool span_is_contiguous() const {
    return true;
  }

  /* Strides of dimensions */
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_7() const { return 1; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_6() const {
    return m_dim.N7;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_5() const {
    return m_dim.N7 * m_dim.N6;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_4() const {
    return m_dim.N7 * m_dim.N6 * m_dim.N5;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_3() const {
    return m_dim.N7 * m_dim.N6 * m_dim.N5 * m_dim.N4;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_2() const {
    return m_dim.N7 * m_dim.N6 * m_dim.N5 * m_dim.N4 * m_dim.N3;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_1() const {
    return m_dim.N7 * m_dim.N6 * m_dim.N5 * m_dim.N4 * m_dim.N3 * m_dim.N2;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_0() const {
    return m_dim.N7 * m_dim.N6 * m_dim.N5 * m_dim.N4 * m_dim.N3 * m_dim.N2 *
           m_dim.N1;
  }

  // Fill the target unbounded array s with the stride.
  // This method differs from stride() in that it does not write the total
  // length to the last index of the array. Preconditions: s must be an array of
  // dimension_type::rank elements
  // The version of clang-format in CI fails from maybe_unused
  // clang-format off
  template <typename iType>
  KOKKOS_INLINE_FUNCTION iType
  stride_fill([[maybe_unused]] iType* const s) const {
    size_type n = 1;
    if constexpr (7 < dimension_type::rank) {
      s[7] = n;
      n *= m_dim.N7;
    }
    if constexpr (6 < dimension_type::rank) {
      s[6] = n;
      n *= m_dim.N6;
    }
    if constexpr (5 < dimension_type::rank) {
      s[5] = n;
      n *= m_dim.N5;
    }
    if constexpr (4 < dimension_type::rank) {
      s[4] = n;
      n *= m_dim.N4;
    }
    if constexpr (3 < dimension_type::rank) {
      s[3] = n;
      n *= m_dim.N3;
    }
    if constexpr (2 < dimension_type::rank) {
      s[2] = n;
      n *= m_dim.N2;
    }
    if constexpr (1 < dimension_type::rank) {
      s[1] = n;
      n *= m_dim.N1;
    }
    if constexpr (0 < dimension_type::rank) {
      s[0] = n;
    }
    return n * m_dim.N0;
  }
  // clang-format on

  // Fill the target unbounded array s with the stride and the total spanned
  // size. This method differs from stride_fill() in that it writes the total
  // spanned size to the last index of the array. Preconditions: s must be an
  // array of dimension_type::rank + 1 elements
  template <typename iType>
  KOKKOS_INLINE_FUNCTION void stride(iType* const s) const {
    s[dimension_type::rank] = stride_fill(s);
  }

  //----------------------------------------
  // MSVC (16.5.5) + CUDA (10.2) did not generate the defaulted functions
  // correct and errors out during compilation. Same for the other places where
  // I changed this.

#ifdef KOKKOS_IMPL_WINDOWS_CUDA
  KOKKOS_FUNCTION ViewOffset() : m_dim(dimension_type()) {}
  KOKKOS_FUNCTION ViewOffset(const ViewOffset& src) { m_dim = src.m_dim; }
  KOKKOS_FUNCTION ViewOffset& operator=(const ViewOffset& src) {
    m_dim = src.m_dim;
    return *this;
  }
#else

  ViewOffset()                  = default;
  ViewOffset(const ViewOffset&) = default;
  ViewOffset& operator=(const ViewOffset&) = default;
#endif

  template <unsigned TrivialScalarSize>
  KOKKOS_INLINE_FUNCTION constexpr ViewOffset(
      std::integral_constant<unsigned, TrivialScalarSize> const&,
      Kokkos::LayoutRight const& arg_layout)
      : m_dim(arg_layout.dimension[0], 0, 0, 0, 0, 0, 0, 0) {}

  template <class DimRHS>
  KOKKOS_INLINE_FUNCTION constexpr ViewOffset(
      const ViewOffset<DimRHS, Kokkos::LayoutRight, void>& rhs)
      : m_dim(rhs.m_dim.N0, rhs.m_dim.N1, rhs.m_dim.N2, rhs.m_dim.N3,
              rhs.m_dim.N4, rhs.m_dim.N5, rhs.m_dim.N6, rhs.m_dim.N7) {
    static_assert(int(DimRHS::rank) == int(dimension_type::rank),
                  "ViewOffset assignment requires equal rank");
    // Also requires equal static dimensions ...
  }

  template <class DimRHS>
  KOKKOS_INLINE_FUNCTION constexpr ViewOffset(
      const ViewOffset<DimRHS, Kokkos::LayoutLeft, void>& rhs)
      : m_dim(rhs.m_dim.N0, 0, 0, 0, 0, 0, 0, 0) {
    static_assert((DimRHS::rank == 0 && dimension_type::rank == 0) ||
                      (DimRHS::rank == 1 && dimension_type::rank == 1),
                  "ViewOffset LayoutRight and LayoutLeft are only compatible "
                  "when rank <= 1");
  }

  template <class DimRHS>
  KOKKOS_INLINE_FUNCTION ViewOffset(
      const ViewOffset<DimRHS, Kokkos::LayoutStride, void>& rhs)
      : m_dim(rhs.m_dim.N0, 0, 0, 0, 0, 0, 0, 0) {}

  //----------------------------------------
  // Subview construction

  template <class DimRHS>
  KOKKOS_INLINE_FUNCTION constexpr ViewOffset(
      const ViewOffset<DimRHS, Kokkos::LayoutRight, void>&,
      const SubviewExtents<DimRHS::rank, dimension_type::rank>& sub)
      : m_dim(sub.range_extent(0), 0, 0, 0, 0, 0, 0, 0) {
    static_assert((0 == dimension_type::rank_dynamic) ||
                      (1 == dimension_type::rank &&
                       1 == dimension_type::rank_dynamic && 1 <= DimRHS::rank),
                  "ViewOffset subview construction requires compatible rank");
  }
};

//----------------------------------------------------------------------------
// LayoutRight AND ( 1 < rank AND 0 < rank_dynamic ) : has padding / striding
template <class Dimension>
struct ViewOffset<
    Dimension, Kokkos::LayoutRight,
    std::enable_if_t<(1 < Dimension::rank && 0 < Dimension::rank_dynamic)>> {
  using is_mapping_plugin = std::true_type;
  using is_regular        = std::true_type;

  using size_type      = size_t;
  using dimension_type = Dimension;
  using array_layout   = Kokkos::LayoutRight;

  dimension_type m_dim;
  size_type m_stride;

  //----------------------------------------

  // rank 1
  template <typename I0>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(I0 const& i0) const {
    return i0;
  }

  // rank 2
  template <typename I0, typename I1>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(I0 const& i0,
                                                        I1 const& i1) const {
    return i1 + i0 * m_stride;
  }

  // rank 3
  template <typename I0, typename I1, typename I2>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(I0 const& i0,
                                                        I1 const& i1,
                                                        I2 const& i2) const {
    return i2 + m_dim.N2 * (i1) + i0 * m_stride;
  }

  // rank 4
  template <typename I0, typename I1, typename I2, typename I3>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(I0 const& i0,
                                                        I1 const& i1,
                                                        I2 const& i2,
                                                        I3 const& i3) const {
    return i3 + m_dim.N3 * (i2 + m_dim.N2 * (i1)) + i0 * m_stride;
  }

  // rank 5
  template <typename I0, typename I1, typename I2, typename I3, typename I4>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(I0 const& i0,
                                                        I1 const& i1,
                                                        I2 const& i2,
                                                        I3 const& i3,
                                                        I4 const& i4) const {
    return i4 + m_dim.N4 * (i3 + m_dim.N3 * (i2 + m_dim.N2 * (i1))) +
           i0 * m_stride;
  }

  // rank 6
  template <typename I0, typename I1, typename I2, typename I3, typename I4,
            typename I5>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(
      I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4,
      I5 const& i5) const {
    return i5 +
           m_dim.N5 *
               (i4 + m_dim.N4 * (i3 + m_dim.N3 * (i2 + m_dim.N2 * (i1)))) +
           i0 * m_stride;
  }

  // rank 7
  template <typename I0, typename I1, typename I2, typename I3, typename I4,
            typename I5, typename I6>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(
      I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4,
      I5 const& i5, I6 const& i6) const {
    return i6 +
           m_dim.N6 *
               (i5 + m_dim.N5 *
                         (i4 + m_dim.N4 *
                                   (i3 + m_dim.N3 * (i2 + m_dim.N2 * (i1))))) +
           i0 * m_stride;
  }

  // rank 8
  template <typename I0, typename I1, typename I2, typename I3, typename I4,
            typename I5, typename I6, typename I7>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(
      I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4,
      I5 const& i5, I6 const& i6, I7 const& i7) const {
    return i7 +
           m_dim.N7 *
               (i6 +
                m_dim.N6 *
                    (i5 +
                     m_dim.N5 *
                         (i4 + m_dim.N4 *
                                   (i3 + m_dim.N3 * (i2 + m_dim.N2 * (i1)))))) +
           i0 * m_stride;
  }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  constexpr array_layout layout() const {
    constexpr auto r = dimension_type::rank;
    return array_layout((r > 0 ? m_dim.N0 : KOKKOS_INVALID_INDEX),
                        (r > 1 ? m_dim.N1 : KOKKOS_INVALID_INDEX),
                        (r > 2 ? m_dim.N2 : KOKKOS_INVALID_INDEX),
                        (r > 3 ? m_dim.N3 : KOKKOS_INVALID_INDEX),
                        (r > 4 ? m_dim.N4 : KOKKOS_INVALID_INDEX),
                        (r > 5 ? m_dim.N5 : KOKKOS_INVALID_INDEX),
                        (r > 6 ? m_dim.N6 : KOKKOS_INVALID_INDEX),
                        (r > 7 ? m_dim.N7 : KOKKOS_INVALID_INDEX));
  }

  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_0() const {
    return m_dim.N0;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_1() const {
    return m_dim.N1;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_2() const {
    return m_dim.N2;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_3() const {
    return m_dim.N3;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_4() const {
    return m_dim.N4;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_5() const {
    return m_dim.N5;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_6() const {
    return m_dim.N6;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_7() const {
    return m_dim.N7;
  }

  /* Cardinality of the domain index space */
  KOKKOS_INLINE_FUNCTION
  constexpr size_type size() const {
    return size_type(m_dim.N0) * m_dim.N1 * m_dim.N2 * m_dim.N3 * m_dim.N4 *
           m_dim.N5 * m_dim.N6 * m_dim.N7;
  }

  /* Span of the range space */
  KOKKOS_INLINE_FUNCTION
  constexpr size_type span() const {
    return size() > 0 ? size_type(m_dim.N0) * m_stride : 0;
  }

  KOKKOS_INLINE_FUNCTION constexpr bool span_is_contiguous() const {
    return m_stride == m_dim.N7 * m_dim.N6 * m_dim.N5 * m_dim.N4 * m_dim.N3 *
                           m_dim.N2 * m_dim.N1;
  }

  /* Strides of dimensions */
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_7() const { return 1; }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_6() const {
    return m_dim.N7;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_5() const {
    return m_dim.N7 * m_dim.N6;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_4() const {
    return m_dim.N7 * m_dim.N6 * m_dim.N5;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_3() const {
    return m_dim.N7 * m_dim.N6 * m_dim.N5 * m_dim.N4;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_2() const {
    return m_dim.N7 * m_dim.N6 * m_dim.N5 * m_dim.N4 * m_dim.N3;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_1() const {
    return m_dim.N7 * m_dim.N6 * m_dim.N5 * m_dim.N4 * m_dim.N3 * m_dim.N2;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_0() const {
    return m_stride;
  }

  // Fill the target unbounded array s with the stride.
  // This method differs from stride() in that it does not write the total
  // length to the last index of the array. Preconditions: s must be an array of
  // dimension_type::rank elements
  // The version of clang-format in CI fails from maybe_unused
  // clang-format off
  template <typename iType>
  KOKKOS_INLINE_FUNCTION iType
  stride_fill([[maybe_unused]] iType* const s) const {
    size_type n = 1;
    if constexpr (7 < dimension_type::rank) {
      s[7] = n;
      n *= m_dim.N7;
    }
    if constexpr (6 < dimension_type::rank) {
      s[6] = n;
      n *= m_dim.N6;
    }
    if constexpr (5 < dimension_type::rank) {
      s[5] = n;
      n *= m_dim.N5;
    }
    if constexpr (4 < dimension_type::rank) {
      s[4] = n;
      n *= m_dim.N4;
    }
    if constexpr (3 < dimension_type::rank) {
      s[3] = n;
      n *= m_dim.N3;
    }
    if constexpr (2 < dimension_type::rank) {
      s[2] = n;
      n *= m_dim.N2;
    }
    if constexpr (1 < dimension_type::rank) {
      s[1] = n;
    }
    if constexpr (0 < dimension_type::rank) {
      s[0] = m_stride;
    }
    return m_stride * m_dim.N0;
  }
  // clang-format on

  // Fill the target unbounded array s with the stride and the total spanned
  // size. This method differs from stride_fill() in that it writes the total
  // spanned size to the last index of the array. Preconditions: s must be an
  // array of dimension_type::rank + 1 elements
  template <typename iType>
  KOKKOS_INLINE_FUNCTION void stride(iType* const s) const {
    s[dimension_type::rank] = stride_fill(s);
  }

  //----------------------------------------

 private:
  template <unsigned TrivialScalarSize>
  struct Padding {
    enum {
      div = TrivialScalarSize == 0
                ? 0
                : Kokkos::Impl::MEMORY_ALIGNMENT /
                      (TrivialScalarSize ? TrivialScalarSize : 1)
    };
    enum {
      mod = TrivialScalarSize == 0
                ? 0
                : Kokkos::Impl::MEMORY_ALIGNMENT %
                      (TrivialScalarSize ? TrivialScalarSize : 1)
    };

    // If memory alignment is a multiple of the trivial scalar size then attempt
    // to align.
    enum { align = 0 != TrivialScalarSize && 0 == mod ? div : 0 };
    enum {
      div_ok = (div != 0) ? div : 1
    };  // To valid modulo zero in constexpr

    KOKKOS_INLINE_FUNCTION
    static constexpr size_t stride(size_t const N) {
      return ((align != 0) &&
              ((static_cast<int>(Kokkos::Impl::MEMORY_ALIGNMENT_THRESHOLD) *
                static_cast<int>(align)) < N) &&
              ((N % div_ok) != 0))
                 ? N + align - (N % div_ok)
                 : N;
    }
  };

 public:
  // MSVC (16.5.5) + CUDA (10.2) did not generate the defaulted functions
  // correct and errors out during compilation. Same for the other places where
  // I changed this.

#ifdef KOKKOS_IMPL_WINDOWS_CUDA
  KOKKOS_FUNCTION ViewOffset() : m_dim(dimension_type()), m_stride(0) {}
  KOKKOS_FUNCTION ViewOffset(const ViewOffset& src) {
    m_dim    = src.m_dim;
    m_stride = src.m_stride;
  }
  KOKKOS_FUNCTION ViewOffset& operator=(const ViewOffset& src) {
    m_dim    = src.m_dim;
    m_stride = src.m_stride;
    return *this;
  }
#else

  ViewOffset()                  = default;
  ViewOffset(const ViewOffset&) = default;
  ViewOffset& operator=(const ViewOffset&) = default;
#endif

  /* Enable padding for trivial scalar types with non-zero trivial scalar size.
   */
  template <unsigned TrivialScalarSize>
  KOKKOS_INLINE_FUNCTION constexpr ViewOffset(
      std::integral_constant<unsigned, TrivialScalarSize> const&,
      Kokkos::LayoutRight const& arg_layout)
      : m_dim(arg_layout.dimension[0], arg_layout.dimension[1],
              arg_layout.dimension[2], arg_layout.dimension[3],
              arg_layout.dimension[4], arg_layout.dimension[5],
              arg_layout.dimension[6], arg_layout.dimension[7]),
        m_stride(
            Padding<TrivialScalarSize>::
                stride(/* 2 <= rank */
                       m_dim.N1 *
                       (dimension_type::rank == 2
                            ? size_t(1)
                            : m_dim.N2 *
                                  (dimension_type::rank == 3
                                       ? size_t(1)
                                       : m_dim.N3 *
                                             (dimension_type::rank == 4
                                                  ? size_t(1)
                                                  : m_dim.N4 *
                                                        (dimension_type::rank ==
                                                                 5
                                                             ? size_t(1)
                                                             : m_dim.N5 *
                                                                   (dimension_type::
                                                                                rank ==
                                                                            6
                                                                        ? size_t(
                                                                              1)
                                                                        : m_dim.N6 *
                                                                              (dimension_type::
                                                                                           rank ==
                                                                                       7
                                                                                   ? size_t(
                                                                                         1)
                                                                                   : m_dim
                                                                                         .N7)))))))) {
  }

  template <class DimRHS>
  KOKKOS_INLINE_FUNCTION constexpr ViewOffset(
      const ViewOffset<DimRHS, Kokkos::LayoutRight, void>& rhs)
      : m_dim(rhs.m_dim.N0, rhs.m_dim.N1, rhs.m_dim.N2, rhs.m_dim.N3,
              rhs.m_dim.N4, rhs.m_dim.N5, rhs.m_dim.N6, rhs.m_dim.N7),
        m_stride(rhs.stride_0()) {
    static_assert(int(DimRHS::rank) == int(dimension_type::rank),
                  "ViewOffset assignment requires equal rank");
    // Also requires equal static dimensions ...
  }

  template <class DimRHS>
  KOKKOS_INLINE_FUNCTION ViewOffset(
      const ViewOffset<DimRHS, Kokkos::LayoutStride, void>& rhs)
      : m_dim(rhs.m_dim.N0, rhs.m_dim.N1, rhs.m_dim.N2, rhs.m_dim.N3,
              rhs.m_dim.N4, rhs.m_dim.N5, rhs.m_dim.N6, rhs.m_dim.N7),
        m_stride(rhs.stride_0()) {
    if (((dimension_type::rank == 2)
             ? rhs.m_stride.S1
             : ((dimension_type::rank == 3)
                    ? rhs.m_stride.S2
                    : ((dimension_type::rank == 4)
                           ? rhs.m_stride.S3
                           : ((dimension_type::rank == 5)
                                  ? rhs.m_stride.S4
                                  : ((dimension_type::rank == 6)
                                         ? rhs.m_stride.S5
                                         : ((dimension_type::rank == 7)
                                                ? rhs.m_stride.S6
                                                : rhs.m_stride.S7)))))) != 1) {
      Kokkos::abort(
          "Kokkos::Impl::ViewOffset assignment of LayoutRight from "
          "LayoutStride requires right-most stride == 1");
    }
  }

  //----------------------------------------
  // Subview construction
  // Last dimension must be non-zero

  template <class DimRHS>
  KOKKOS_INLINE_FUNCTION constexpr ViewOffset(
      const ViewOffset<DimRHS, Kokkos::LayoutRight, void>& rhs,
      const SubviewExtents<DimRHS::rank, dimension_type::rank>& sub)
      : m_dim(sub.range_extent(0), sub.range_extent(1), sub.range_extent(2),
              sub.range_extent(3), sub.range_extent(4), sub.range_extent(5),
              sub.range_extent(6), sub.range_extent(7)),
        m_stride(
            0 == sub.range_index(0)
                ? rhs.stride_0()
                : (1 == sub.range_index(0)
                       ? rhs.stride_1()
                       : (2 == sub.range_index(0)
                              ? rhs.stride_2()
                              : (3 == sub.range_index(0)
                                     ? rhs.stride_3()
                                     : (4 == sub.range_index(0)
                                            ? rhs.stride_4()
                                            : (5 == sub.range_index(0)
                                                   ? rhs.stride_5()
                                                   : (6 == sub.range_index(0)
                                                          ? rhs.stride_6()
                                                          : 0))))))) {
    /*      // This subview must be 2 == rank and 2 == rank_dynamic
          // due to only having stride #0.
          // The source dimension #0 must be non-zero for stride-one leading
       dimension.
          // At most subsequent dimension can be non-zero.

          static_assert( (( 2 == dimension_type::rank ) &&
                          ( 2 <= DimRHS::rank )) ||
                         ()
                       , "ViewOffset subview construction requires compatible
       rank" );
    */
  }
};

//----------------------------------------------------------------------------
/* Strided array layout only makes sense for 0 < rank */
/* rank = 0 included for DynRankView case */

template <unsigned Rank>
struct ViewStride;

template <>
struct ViewStride<0> {
  static constexpr size_t S0 = 0, S1 = 0, S2 = 0, S3 = 0, S4 = 0, S5 = 0,
                          S6 = 0, S7 = 0;

  ViewStride()                  = default;
  ViewStride(const ViewStride&) = default;
  ViewStride& operator=(const ViewStride&) = default;

  KOKKOS_INLINE_FUNCTION
  constexpr ViewStride(size_t, size_t, size_t, size_t, size_t, size_t, size_t,
                       size_t) {}
};

template <>
struct ViewStride<1> {
  size_t S0;
  static constexpr size_t S1 = 0, S2 = 0, S3 = 0, S4 = 0, S5 = 0, S6 = 0,
                          S7 = 0;

  ViewStride()                  = default;
  ViewStride(const ViewStride&) = default;
  ViewStride& operator=(const ViewStride&) = default;

  KOKKOS_INLINE_FUNCTION
  constexpr ViewStride(size_t aS0, size_t, size_t, size_t, size_t, size_t,
                       size_t, size_t)
      : S0(aS0) {}
};

template <>
struct ViewStride<2> {
  size_t S0, S1;
  static constexpr size_t S2 = 0, S3 = 0, S4 = 0, S5 = 0, S6 = 0, S7 = 0;

  ViewStride()                  = default;
  ViewStride(const ViewStride&) = default;
  ViewStride& operator=(const ViewStride&) = default;

  KOKKOS_INLINE_FUNCTION
  constexpr ViewStride(size_t aS0, size_t aS1, size_t, size_t, size_t, size_t,
                       size_t, size_t)
      : S0(aS0), S1(aS1) {}
};

template <>
struct ViewStride<3> {
  size_t S0, S1, S2;
  static constexpr size_t S3 = 0, S4 = 0, S5 = 0, S6 = 0, S7 = 0;

  ViewStride()                  = default;
  ViewStride(const ViewStride&) = default;
  ViewStride& operator=(const ViewStride&) = default;

  KOKKOS_INLINE_FUNCTION
  constexpr ViewStride(size_t aS0, size_t aS1, size_t aS2, size_t, size_t,
                       size_t, size_t, size_t)
      : S0(aS0), S1(aS1), S2(aS2) {}
};

template <>
struct ViewStride<4> {
  size_t S0, S1, S2, S3;
  static constexpr size_t S4 = 0, S5 = 0, S6 = 0, S7 = 0;

  ViewStride()                  = default;
  ViewStride(const ViewStride&) = default;
  ViewStride& operator=(const ViewStride&) = default;

  KOKKOS_INLINE_FUNCTION
  constexpr ViewStride(size_t aS0, size_t aS1, size_t aS2, size_t aS3, size_t,
                       size_t, size_t, size_t)
      : S0(aS0), S1(aS1), S2(aS2), S3(aS3) {}
};

template <>
struct ViewStride<5> {
  size_t S0, S1, S2, S3, S4;
  static constexpr size_t S5 = 0, S6 = 0, S7 = 0;

  ViewStride()                  = default;
  ViewStride(const ViewStride&) = default;
  ViewStride& operator=(const ViewStride&) = default;

  KOKKOS_INLINE_FUNCTION
  constexpr ViewStride(size_t aS0, size_t aS1, size_t aS2, size_t aS3,
                       size_t aS4, size_t, size_t, size_t)
      : S0(aS0), S1(aS1), S2(aS2), S3(aS3), S4(aS4) {}
};

template <>
struct ViewStride<6> {
  size_t S0, S1, S2, S3, S4, S5;
  static constexpr size_t S6 = 0, S7 = 0;

  ViewStride()                  = default;
  ViewStride(const ViewStride&) = default;
  ViewStride& operator=(const ViewStride&) = default;

  KOKKOS_INLINE_FUNCTION
  constexpr ViewStride(size_t aS0, size_t aS1, size_t aS2, size_t aS3,
                       size_t aS4, size_t aS5, size_t, size_t)
      : S0(aS0), S1(aS1), S2(aS2), S3(aS3), S4(aS4), S5(aS5) {}
};

template <>
struct ViewStride<7> {
  size_t S0, S1, S2, S3, S4, S5, S6;
  static constexpr size_t S7 = 0;

  ViewStride()                  = default;
  ViewStride(const ViewStride&) = default;
  ViewStride& operator=(const ViewStride&) = default;

  KOKKOS_INLINE_FUNCTION
  constexpr ViewStride(size_t aS0, size_t aS1, size_t aS2, size_t aS3,
                       size_t aS4, size_t aS5, size_t aS6, size_t)
      : S0(aS0), S1(aS1), S2(aS2), S3(aS3), S4(aS4), S5(aS5), S6(aS6) {}
};

template <>
struct ViewStride<8> {
  size_t S0, S1, S2, S3, S4, S5, S6, S7;

  ViewStride()                  = default;
  ViewStride(const ViewStride&) = default;
  ViewStride& operator=(const ViewStride&) = default;

  KOKKOS_INLINE_FUNCTION
  constexpr ViewStride(size_t aS0, size_t aS1, size_t aS2, size_t aS3,
                       size_t aS4, size_t aS5, size_t aS6, size_t aS7)
      : S0(aS0),
        S1(aS1),
        S2(aS2),
        S3(aS3),
        S4(aS4),
        S5(aS5),
        S6(aS6),
        S7(aS7) {}
};

template <class Dimension>
struct ViewOffset<Dimension, Kokkos::LayoutStride, void> {
 private:
  using stride_type = ViewStride<Dimension::rank>;

 public:
  using is_mapping_plugin = std::true_type;
  using is_regular        = std::true_type;

  using size_type      = size_t;
  using dimension_type = Dimension;
  using array_layout   = Kokkos::LayoutStride;

  dimension_type m_dim;
  stride_type m_stride;

  //----------------------------------------

  // rank 1
  template <typename I0>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(I0 const& i0) const {
    return i0 * m_stride.S0;
  }

  // rank 2
  template <typename I0, typename I1>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(I0 const& i0,
                                                        I1 const& i1) const {
    return i0 * m_stride.S0 + i1 * m_stride.S1;
  }

  // rank 3
  template <typename I0, typename I1, typename I2>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(I0 const& i0,
                                                        I1 const& i1,
                                                        I2 const& i2) const {
    return i0 * m_stride.S0 + i1 * m_stride.S1 + i2 * m_stride.S2;
  }

  // rank 4
  template <typename I0, typename I1, typename I2, typename I3>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(I0 const& i0,
                                                        I1 const& i1,
                                                        I2 const& i2,
                                                        I3 const& i3) const {
    return i0 * m_stride.S0 + i1 * m_stride.S1 + i2 * m_stride.S2 +
           i3 * m_stride.S3;
  }

  // rank 5
  template <typename I0, typename I1, typename I2, typename I3, typename I4>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(I0 const& i0,
                                                        I1 const& i1,
                                                        I2 const& i2,
                                                        I3 const& i3,
                                                        I4 const& i4) const {
    return i0 * m_stride.S0 + i1 * m_stride.S1 + i2 * m_stride.S2 +
           i3 * m_stride.S3 + i4 * m_stride.S4;
  }

  // rank 6
  template <typename I0, typename I1, typename I2, typename I3, typename I4,
            typename I5>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(
      I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4,
      I5 const& i5) const {
    return i0 * m_stride.S0 + i1 * m_stride.S1 + i2 * m_stride.S2 +
           i3 * m_stride.S3 + i4 * m_stride.S4 + i5 * m_stride.S5;
  }

  // rank 7
  template <typename I0, typename I1, typename I2, typename I3, typename I4,
            typename I5, typename I6>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(
      I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4,
      I5 const& i5, I6 const& i6) const {
    return i0 * m_stride.S0 + i1 * m_stride.S1 + i2 * m_stride.S2 +
           i3 * m_stride.S3 + i4 * m_stride.S4 + i5 * m_stride.S5 +
           i6 * m_stride.S6;
  }

  // rank 8
  template <typename I0, typename I1, typename I2, typename I3, typename I4,
            typename I5, typename I6, typename I7>
  KOKKOS_INLINE_FUNCTION constexpr size_type operator()(
      I0 const& i0, I1 const& i1, I2 const& i2, I3 const& i3, I4 const& i4,
      I5 const& i5, I6 const& i6, I7 const& i7) const {
    return i0 * m_stride.S0 + i1 * m_stride.S1 + i2 * m_stride.S2 +
           i3 * m_stride.S3 + i4 * m_stride.S4 + i5 * m_stride.S5 +
           i6 * m_stride.S6 + i7 * m_stride.S7;
  }

  //----------------------------------------

  KOKKOS_INLINE_FUNCTION
  constexpr array_layout layout() const {
    constexpr auto r = dimension_type::rank;
    return array_layout((r > 0 ? m_dim.N0 : KOKKOS_INVALID_INDEX), m_stride.S0,
                        (r > 1 ? m_dim.N1 : KOKKOS_INVALID_INDEX), m_stride.S1,
                        (r > 2 ? m_dim.N2 : KOKKOS_INVALID_INDEX), m_stride.S2,
                        (r > 3 ? m_dim.N3 : KOKKOS_INVALID_INDEX), m_stride.S3,
                        (r > 4 ? m_dim.N4 : KOKKOS_INVALID_INDEX), m_stride.S4,
                        (r > 5 ? m_dim.N5 : KOKKOS_INVALID_INDEX), m_stride.S5,
                        (r > 6 ? m_dim.N6 : KOKKOS_INVALID_INDEX), m_stride.S6,
                        (r > 7 ? m_dim.N7 : KOKKOS_INVALID_INDEX), m_stride.S7);
  }

  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_0() const {
    return m_dim.N0;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_1() const {
    return m_dim.N1;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_2() const {
    return m_dim.N2;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_3() const {
    return m_dim.N3;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_4() const {
    return m_dim.N4;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_5() const {
    return m_dim.N5;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_6() const {
    return m_dim.N6;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type dimension_7() const {
    return m_dim.N7;
  }

  /* Cardinality of the domain index space */
  KOKKOS_INLINE_FUNCTION
  constexpr size_type size() const {
    return dimension_type::rank == 0
               ? 1
               : size_type(m_dim.N0) * m_dim.N1 * m_dim.N2 * m_dim.N3 *
                     m_dim.N4 * m_dim.N5 * m_dim.N6 * m_dim.N7;
  }

 private:
  KOKKOS_INLINE_FUNCTION
  static constexpr size_type Max(size_type lhs, size_type rhs) {
    return lhs < rhs ? rhs : lhs;
  }

 public:
  /* Span of the range space, largest stride * dimension */
  KOKKOS_INLINE_FUNCTION
  constexpr size_type span() const {
    return dimension_type::rank == 0
               ? 1
               : (size() == size_type(0)
                      ? size_type(0)
                      : Max(m_dim.N0 * m_stride.S0,
                            Max(m_dim.N1 * m_stride.S1,
                                Max(m_dim.N2 * m_stride.S2,
                                    Max(m_dim.N3 * m_stride.S3,
                                        Max(m_dim.N4 * m_stride.S4,
                                            Max(m_dim.N5 * m_stride.S5,
                                                Max(m_dim.N6 * m_stride.S6,
                                                    m_dim.N7 *
                                                        m_stride.S7))))))));
  }

  KOKKOS_INLINE_FUNCTION constexpr bool span_is_contiguous() const {
    return span() == size();
  }

  /* Strides of dimensions */
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_0() const {
    return m_stride.S0;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_1() const {
    return m_stride.S1;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_2() const {
    return m_stride.S2;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_3() const {
    return m_stride.S3;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_4() const {
    return m_stride.S4;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_5() const {
    return m_stride.S5;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_6() const {
    return m_stride.S6;
  }
  KOKKOS_INLINE_FUNCTION constexpr size_type stride_7() const {
    return m_stride.S7;
  }

  // Fill the target unbounded array s with the stride.
  // This method differs from stride() in that it does not write the total
  // length to the last index of the array. Preconditions: s must be an array of
  // dimension_type::rank elements
  // The version of clang-format in CI fails from maybe_unused
  // clang-format off
  template <typename iType>
  KOKKOS_INLINE_FUNCTION iType
  stride_fill([[maybe_unused]] iType* const s) const {
    if constexpr (0 < dimension_type::rank) {
      s[0] = m_stride.S0;
    }
    if constexpr (1 < dimension_type::rank) {
      s[1] = m_stride.S1;
    }
    if constexpr (2 < dimension_type::rank) {
      s[2] = m_stride.S2;
    }
    if constexpr (3 < dimension_type::rank) {
      s[3] = m_stride.S3;
    }
    if constexpr (4 < dimension_type::rank) {
      s[4] = m_stride.S4;
    }
    if constexpr (5 < dimension_type::rank) {
      s[5] = m_stride.S5;
    }
    if constexpr (6 < dimension_type::rank) {
      s[6] = m_stride.S6;
    }
    if constexpr (7 < dimension_type::rank) {
      s[7] = m_stride.S7;
    }
    return span();
  }
  // clang-format on

  // Fill the target unbounded array s with the stride and the total spanned
  // size. This method differs from stride_fill() in that it writes the total
  // spanned size to the last index of the array. Preconditions: s must be an
  // array of dimension_type::rank + 1 elements
  template <typename iType>
  KOKKOS_INLINE_FUNCTION void stride(iType* const s) const {
    s[dimension_type::rank] = stride_fill(s);
  }

  //----------------------------------------
  // MSVC (16.5.5) + CUDA (10.2) did not generate the defaulted functions
  // correct and errors out during compilation. Same for the other places where
  // I changed this.

#ifdef KOKKOS_IMPL_WINDOWS_CUDA
  KOKKOS_FUNCTION ViewOffset()
      : m_dim(dimension_type()), m_stride(stride_type()) {}
  KOKKOS_FUNCTION ViewOffset(const ViewOffset& src) {
    m_dim    = src.m_dim;
    m_stride = src.m_stride;
  }
  KOKKOS_FUNCTION ViewOffset& operator=(const ViewOffset& src) {
    m_dim    = src.m_dim;
    m_stride = src.m_stride;
    return *this;
  }
#else

  ViewOffset()                  = default;
  ViewOffset(const ViewOffset&) = default;
  ViewOffset& operator=(const ViewOffset&) = default;
#endif

  KOKKOS_INLINE_FUNCTION
  constexpr ViewOffset(std::integral_constant<unsigned, 0> const&,
                       Kokkos::LayoutStride const& rhs)
      : m_dim(rhs.dimension[0], rhs.dimension[1], rhs.dimension[2],
              rhs.dimension[3], rhs.dimension[4], rhs.dimension[5],
              rhs.dimension[6], rhs.dimension[7]),
        m_stride(rhs.stride[0], rhs.stride[1], rhs.stride[2], rhs.stride[3],
                 rhs.stride[4], rhs.stride[5], rhs.stride[6], rhs.stride[7]) {}

  template <class DimRHS, class LayoutRHS>
  KOKKOS_INLINE_FUNCTION constexpr ViewOffset(
      const ViewOffset<DimRHS, LayoutRHS, void>& rhs)
      : m_dim(rhs.m_dim.N0, rhs.m_dim.N1, rhs.m_dim.N2, rhs.m_dim.N3,
              rhs.m_dim.N4, rhs.m_dim.N5, rhs.m_dim.N6, rhs.m_dim.N7),
        m_stride(rhs.stride_0(), rhs.stride_1(), rhs.stride_2(), rhs.stride_3(),
                 rhs.stride_4(), rhs.stride_5(), rhs.stride_6(),
                 rhs.stride_7()) {
    static_assert(int(DimRHS::rank) == int(dimension_type::rank),
                  "ViewOffset assignment requires equal rank");
    // Also requires equal static dimensions ...
  }

  //----------------------------------------
  // Subview construction

 private:
  template <class DimRHS, class LayoutRHS>
  KOKKOS_INLINE_FUNCTION static constexpr size_t stride(
      unsigned r, const ViewOffset<DimRHS, LayoutRHS, void>& rhs) {
    return r > 7
               ? 0
               : (r == 0
                      ? rhs.stride_0()
                      : (r == 1
                             ? rhs.stride_1()
                             : (r == 2
                                    ? rhs.stride_2()
                                    : (r == 3
                                           ? rhs.stride_3()
                                           : (r == 4
                                                  ? rhs.stride_4()
                                                  : (r == 5
                                                         ? rhs.stride_5()
                                                         : (r == 6
                                                                ? rhs.stride_6()
                                                                : rhs.stride_7())))))));
  }

 public:
  template <class DimRHS, class LayoutRHS>
  KOKKOS_INLINE_FUNCTION constexpr ViewOffset(
      const ViewOffset<DimRHS, LayoutRHS, void>& rhs,
      const SubviewExtents<DimRHS::rank, dimension_type::rank>& sub)
      // range_extent(r) returns 0 when dimension_type::rank <= r
      : m_dim(sub.range_extent(0), sub.range_extent(1), sub.range_extent(2),
              sub.range_extent(3), sub.range_extent(4), sub.range_extent(5),
              sub.range_extent(6), sub.range_extent(7))
        // range_index(r) returns ~0u when dimension_type::rank <= r
        ,
        m_stride(
            stride(sub.range_index(0), rhs), stride(sub.range_index(1), rhs),
            stride(sub.range_index(2), rhs), stride(sub.range_index(3), rhs),
            stride(sub.range_index(4), rhs), stride(sub.range_index(5), rhs),
            stride(sub.range_index(6), rhs), stride(sub.range_index(7), rhs)) {}
};

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

/** \brief  ViewDataHandle provides the type of the 'data handle' which the view
 *          uses to access data with the [] operator. It also provides
 *          an allocate function and a function to extract a raw ptr from the
 *          data handle. ViewDataHandle also defines an enum ReferenceAble which
 *          specifies whether references/pointers to elements can be taken and a
 *          'return_type' which is what the view operators will give back.
 *          Specialisation of this object allows three things depending
 *          on ViewTraits and compiler options:
 *          (i)   Use special allocator (e.g. huge pages/small pages and pinned
 * memory) (ii)  Use special data handle type (e.g. add Cuda Texture Object)
 *          (iii) Use special access intrinsics (e.g. texture fetch and
 * non-caching loads)
 */
template <class Traits, class Enable = void>
struct ViewDataHandle {
  using value_type  = typename Traits::value_type;
  using handle_type = typename Traits::value_type*;
  using return_type = typename Traits::value_type&;
  using track_type  = Kokkos::Impl::SharedAllocationTracker;

  KOKKOS_INLINE_FUNCTION
  static handle_type assign(value_type* arg_data_ptr,
                            track_type const& /*arg_tracker*/) {
    return handle_type(arg_data_ptr);
  }

  KOKKOS_INLINE_FUNCTION
  static handle_type assign(handle_type const arg_data_ptr, size_t offset) {
    return handle_type(arg_data_ptr + offset);
  }
};

template <class Traits>
struct ViewDataHandle<
    Traits,
    std::enable_if_t<(std::is_same<typename Traits::non_const_value_type,
                                   typename Traits::value_type>::value &&
                      std::is_void<typename Traits::specialize>::value &&
                      Traits::memory_traits::is_atomic)>> {
  using value_type  = typename Traits::value_type;
  using handle_type = typename Kokkos::Impl::AtomicViewDataHandle<Traits>;
  using return_type = typename Kokkos::Impl::AtomicDataElement<Traits>;
  using track_type  = Kokkos::Impl::SharedAllocationTracker;

  KOKKOS_INLINE_FUNCTION
  static handle_type assign(value_type* arg_data_ptr,
                            track_type const& /*arg_tracker*/) {
    return handle_type(arg_data_ptr);
  }

  template <class SrcHandleType>
  KOKKOS_INLINE_FUNCTION static handle_type assign(
      const SrcHandleType& arg_handle, size_t offset) {
    return handle_type(arg_handle + offset);
  }
};

template <class Traits>
struct ViewDataHandle<
    Traits,
    std::enable_if_t<(std::is_void<typename Traits::specialize>::value &&
                      (!Traits::memory_traits::is_aligned) &&
                      Traits::memory_traits::is_restrict &&
                      (!Traits::memory_traits::is_atomic))>> {
  using value_type  = typename Traits::value_type;
  using handle_type = typename Traits::value_type* KOKKOS_RESTRICT;
  using return_type = typename Traits::value_type& KOKKOS_RESTRICT;
  using track_type  = Kokkos::Impl::SharedAllocationTracker;

  KOKKOS_INLINE_FUNCTION
  static value_type* assign(value_type* arg_data_ptr,
                            track_type const& /*arg_tracker*/) {
    return (value_type*)(arg_data_ptr);
  }

  KOKKOS_INLINE_FUNCTION
  static value_type* assign(handle_type const arg_data_ptr, size_t offset) {
    return (value_type*)(arg_data_ptr + offset);
  }
};

template <class Traits>
struct ViewDataHandle<
    Traits,
    std::enable_if_t<(std::is_void<typename Traits::specialize>::value &&
                      Traits::memory_traits::is_aligned &&
                      (!Traits::memory_traits::is_restrict) &&
                      (!Traits::memory_traits::is_atomic))>> {
  using value_type = typename Traits::value_type;
  // typedef work-around for intel compilers error #3186: expected typedef
  // declaration
  // NOLINTNEXTLINE(modernize-use-using)
  typedef value_type* KOKKOS_IMPL_ALIGN_PTR(KOKKOS_MEMORY_ALIGNMENT)
      handle_type;
  using return_type = typename Traits::value_type&;
  using track_type  = Kokkos::Impl::SharedAllocationTracker;

  KOKKOS_INLINE_FUNCTION
  static handle_type assign(value_type* arg_data_ptr,
                            track_type const& /*arg_tracker*/) {
    if (reinterpret_cast<uintptr_t>(arg_data_ptr) % Impl::MEMORY_ALIGNMENT) {
      Kokkos::abort(
          "Assigning NonAligned View or Pointer to Kokkos::View with Aligned "
          "attribute");
    }
    return handle_type(arg_data_ptr);
  }

  KOKKOS_INLINE_FUNCTION
  static handle_type assign(handle_type const arg_data_ptr, size_t offset) {
    if (reinterpret_cast<uintptr_t>(arg_data_ptr + offset) %
        Impl::MEMORY_ALIGNMENT) {
      Kokkos::abort(
          "Assigning NonAligned View or Pointer to Kokkos::View with Aligned "
          "attribute");
    }
    return handle_type(arg_data_ptr + offset);
  }
};

template <class Traits>
struct ViewDataHandle<
    Traits,
    std::enable_if_t<(std::is_void<typename Traits::specialize>::value &&
                      Traits::memory_traits::is_aligned &&
                      Traits::memory_traits::is_restrict &&
                      (!Traits::memory_traits::is_atomic))>> {
  using value_type = typename Traits::value_type;
  // typedef work-around for intel compilers error #3186: expected typedef
  // declaration
  // NOLINTNEXTLINE(modernize-use-using)
  typedef value_type* KOKKOS_IMPL_ALIGN_PTR(KOKKOS_MEMORY_ALIGNMENT)
      handle_type;
  using return_type = typename Traits::value_type& KOKKOS_RESTRICT;
  using track_type  = Kokkos::Impl::SharedAllocationTracker;

  KOKKOS_INLINE_FUNCTION
  static value_type* assign(value_type* arg_data_ptr,
                            track_type const& /*arg_tracker*/) {
    if (reinterpret_cast<uintptr_t>(arg_data_ptr) % Impl::MEMORY_ALIGNMENT) {
      Kokkos::abort(
          "Assigning NonAligned View or Pointer to Kokkos::View with Aligned "
          "attribute");
    }
    return (value_type*)(arg_data_ptr);
  }

  KOKKOS_INLINE_FUNCTION
  static value_type* assign(handle_type const arg_data_ptr, size_t offset) {
    if (reinterpret_cast<uintptr_t>(arg_data_ptr + offset) %
        Impl::MEMORY_ALIGNMENT) {
      Kokkos::abort(
          "Assigning NonAligned View or Pointer to Kokkos::View with Aligned "
          "attribute");
    }
    return (value_type*)(arg_data_ptr + offset);
  }
};
}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {
//----------------------------------------------------------------------------
/** \brief  View mapping for non-specialized data type and standard layout */
template <class Traits>
class ViewMapping<
    Traits,
    std::enable_if_t<(
        std::is_void<typename Traits::specialize>::value &&
        ViewOffset<typename Traits::dimension, typename Traits::array_layout,
                   void>::is_mapping_plugin::value)>> {
 public:
  using offset_type = ViewOffset<typename Traits::dimension,
                                 typename Traits::array_layout, void>;

  using handle_type = typename ViewDataHandle<Traits>::handle_type;

  handle_type m_impl_handle;
  offset_type m_impl_offset;

 private:
  template <class, class...>
  friend class ViewMapping;

  KOKKOS_INLINE_FUNCTION
  ViewMapping(const handle_type& arg_handle, const offset_type& arg_offset)
      : m_impl_handle(arg_handle), m_impl_offset(arg_offset) {}

 public:
  using printable_label_typedef = void;
  enum { is_managed = Traits::is_managed };

  //----------------------------------------
  // Domain dimensions

  static constexpr unsigned Rank = Traits::dimension::rank;

  template <typename iType>
  KOKKOS_INLINE_FUNCTION constexpr size_t extent(const iType& r) const {
    return m_impl_offset.m_dim.extent(r);
  }

  static KOKKOS_INLINE_FUNCTION constexpr size_t static_extent(
      const unsigned r) noexcept {
    using dim_type = typename offset_type::dimension_type;
    return dim_type::static_extent(r);
  }

  KOKKOS_INLINE_FUNCTION constexpr typename Traits::array_layout layout()
      const {
    return m_impl_offset.layout();
  }

  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_0() const {
    return m_impl_offset.dimension_0();
  }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_1() const {
    return m_impl_offset.dimension_1();
  }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_2() const {
    return m_impl_offset.dimension_2();
  }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_3() const {
    return m_impl_offset.dimension_3();
  }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_4() const {
    return m_impl_offset.dimension_4();
  }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_5() const {
    return m_impl_offset.dimension_5();
  }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_6() const {
    return m_impl_offset.dimension_6();
  }
  KOKKOS_INLINE_FUNCTION constexpr size_t dimension_7() const {
    return m_impl_offset.dimension_7();
  }

  // Is a regular layout with uniform striding for each index.
  using is_regular = typename offset_type::is_regular;

  KOKKOS_INLINE_FUNCTION constexpr size_t stride_0() const {
    return m_impl_offset.stride_0();
  }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_1() const {
    return m_impl_offset.stride_1();
  }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_2() const {
    return m_impl_offset.stride_2();
  }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_3() const {
    return m_impl_offset.stride_3();
  }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_4() const {
    return m_impl_offset.stride_4();
  }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_5() const {
    return m_impl_offset.stride_5();
  }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_6() const {
    return m_impl_offset.stride_6();
  }
  KOKKOS_INLINE_FUNCTION constexpr size_t stride_7() const {
    return m_impl_offset.stride_7();
  }

  // Fill the target unbounded array s with the stride and the total spanned
  // size. This method differs from stride_fill() in that it writes the total
  // spanned size to the last index of the array. Preconditions: s must be an
  // array of dimension_type::rank + 1 elements
  template <typename iType>
  KOKKOS_INLINE_FUNCTION void stride(iType* const s) const {
    m_impl_offset.stride(s);
  }

  // Fill the target unbounded array s with the stride.
  // This method differs from stride() in that it does not write the total
  // length to the last index of the array. Preconditions: s must be an array of
  // dimension_type::rank elements
  template <typename iType>
  KOKKOS_INLINE_FUNCTION iType stride_fill(iType* const s) const {
    return m_impl_offset.stride_fill(s);
  }

  //----------------------------------------
  // Range span

  /** \brief  Span of the mapped range */
  KOKKOS_INLINE_FUNCTION constexpr size_t span() const {
    return m_impl_offset.span();
  }

  /** \brief  Is the mapped range span contiguous */
  KOKKOS_INLINE_FUNCTION constexpr bool span_is_contiguous() const {
    return m_impl_offset.span_is_contiguous();
  }

  using reference_type = typename ViewDataHandle<Traits>::return_type;
  using pointer_type   = typename Traits::value_type*;

  /** \brief  Query raw pointer to memory */
  KOKKOS_INLINE_FUNCTION constexpr pointer_type data() const {
    return m_impl_handle;
  }

  //----------------------------------------
  // The View class performs all rank and bounds checking before
  // calling these element reference methods.

  KOKKOS_FORCEINLINE_FUNCTION
  reference_type reference() const { return m_impl_handle[0]; }

  template <typename I0>
  KOKKOS_FORCEINLINE_FUNCTION
      std::enable_if_t<(std::is_integral<I0>::value &&
                        // if layout is neither stride nor irregular,
                        // then just use the handle directly
                        !(std::is_same<typename Traits::array_layout,
                                       Kokkos::LayoutStride>::value ||
                          !is_regular::value)),
                       reference_type>
      reference(const I0& i0) const {
    return m_impl_handle[i0];
  }

  template <typename I0>
  KOKKOS_FORCEINLINE_FUNCTION
      std::enable_if_t<(std::is_integral<I0>::value &&
                        // if the layout is strided or irregular, then
                        // we have to use the offset
                        (std::is_same<typename Traits::array_layout,
                                      Kokkos::LayoutStride>::value ||
                         !is_regular::value)),
                       reference_type>
      reference(const I0& i0) const {
    return m_impl_handle[m_impl_offset(i0)];
  }

  template <typename I0, typename I1>
  KOKKOS_FORCEINLINE_FUNCTION reference_type reference(const I0& i0,
                                                       const I1& i1) const {
    return m_impl_handle[m_impl_offset(i0, i1)];
  }

  template <typename I0, typename I1, typename I2>
  KOKKOS_FORCEINLINE_FUNCTION reference_type reference(const I0& i0,
                                                       const I1& i1,
                                                       const I2& i2) const {
    return m_impl_handle[m_impl_offset(i0, i1, i2)];
  }

  template <typename I0, typename I1, typename I2, typename I3>
  KOKKOS_FORCEINLINE_FUNCTION reference_type
  reference(const I0& i0, const I1& i1, const I2& i2, const I3& i3) const {
    return m_impl_handle[m_impl_offset(i0, i1, i2, i3)];
  }

  template <typename I0, typename I1, typename I2, typename I3, typename I4>
  KOKKOS_FORCEINLINE_FUNCTION reference_type reference(const I0& i0,
                                                       const I1& i1,
                                                       const I2& i2,
                                                       const I3& i3,
                                                       const I4& i4) const {
    return m_impl_handle[m_impl_offset(i0, i1, i2, i3, i4)];
  }

  template <typename I0, typename I1, typename I2, typename I3, typename I4,
            typename I5>
  KOKKOS_FORCEINLINE_FUNCTION reference_type
  reference(const I0& i0, const I1& i1, const I2& i2, const I3& i3,
            const I4& i4, const I5& i5) const {
    return m_impl_handle[m_impl_offset(i0, i1, i2, i3, i4, i5)];
  }

  template <typename I0, typename I1, typename I2, typename I3, typename I4,
            typename I5, typename I6>
  KOKKOS_FORCEINLINE_FUNCTION reference_type
  reference(const I0& i0, const I1& i1, const I2& i2, const I3& i3,
            const I4& i4, const I5& i5, const I6& i6) const {
    return m_impl_handle[m_impl_offset(i0, i1, i2, i3, i4, i5, i6)];
  }

  template <typename I0, typename I1, typename I2, typename I3, typename I4,
            typename I5, typename I6, typename I7>
  KOKKOS_FORCEINLINE_FUNCTION reference_type
  reference(const I0& i0, const I1& i1, const I2& i2, const I3& i3,
            const I4& i4, const I5& i5, const I6& i6, const I7& i7) const {
    return m_impl_handle[m_impl_offset(i0, i1, i2, i3, i4, i5, i6, i7)];
  }

  //----------------------------------------

 private:
  enum { MemorySpanMask = 8 - 1 /* Force alignment on 8 byte boundary */ };
  enum { MemorySpanSize = sizeof(typename Traits::value_type) };

 public:
  /** \brief  Span, in bytes, of the referenced memory */
  KOKKOS_INLINE_FUNCTION constexpr size_t memory_span() const {
    return (m_impl_offset.span() * sizeof(typename Traits::value_type) +
            MemorySpanMask) &
           ~size_t(MemorySpanMask);
  }

  //----------------------------------------

  KOKKOS_DEFAULTED_FUNCTION ~ViewMapping() = default;
  KOKKOS_INLINE_FUNCTION ViewMapping() : m_impl_handle(), m_impl_offset() {}

  KOKKOS_DEFAULTED_FUNCTION ViewMapping(const ViewMapping&) = default;
  KOKKOS_DEFAULTED_FUNCTION ViewMapping& operator=(const ViewMapping&) =
      default;

  KOKKOS_DEFAULTED_FUNCTION ViewMapping(ViewMapping&&) = default;
  KOKKOS_DEFAULTED_FUNCTION ViewMapping& operator=(ViewMapping&&) = default;

  //----------------------------------------

  /**\brief  Span, in bytes, of the required memory */
  KOKKOS_INLINE_FUNCTION
  static constexpr size_t memory_span(
      typename Traits::array_layout const& arg_layout) {
    using padding = std::integral_constant<unsigned int, 0>;
    return (offset_type(padding(), arg_layout).span() * MemorySpanSize +
            MemorySpanMask) &
           ~size_t(MemorySpanMask);
  }

  /**\brief  Wrap a span of memory */
  template <class... P>
  KOKKOS_INLINE_FUNCTION ViewMapping(
      Kokkos::Impl::ViewCtorProp<P...> const& arg_prop,
      typename Traits::array_layout const& arg_layout)
      : m_impl_handle(Impl::get_property<Impl::PointerTag>(arg_prop)),
        m_impl_offset(std::integral_constant<unsigned, 0>(), arg_layout) {}

  /**\brief  Assign data */
  KOKKOS_INLINE_FUNCTION
  void assign_data(pointer_type arg_ptr) {
    m_impl_handle = handle_type(arg_ptr);
  }

  //----------------------------------------
  /*  Allocate and construct mapped array.
   *  Allocate via shared allocation record and
   *  return that record for allocation tracking.
   */
  template <class... P>
  Kokkos::Impl::SharedAllocationRecord<>* allocate_shared(
      Kokkos::Impl::ViewCtorProp<P...> const& arg_prop,
      typename Traits::array_layout const& arg_layout,
      bool execution_space_specified) {
    using alloc_prop = Kokkos::Impl::ViewCtorProp<P...>;

    using execution_space = typename alloc_prop::execution_space;
    using memory_space    = typename Traits::memory_space;
    static_assert(
        SpaceAccessibility<execution_space, memory_space>::accessible);
    using device_type  = Kokkos::Device<execution_space, memory_space>;
    using value_type   = typename Traits::value_type;
    using functor_type = std::conditional_t<
        alloc_prop::sequential_host_init,
        ViewValueFunctorSequentialHostInit<device_type, value_type>,
        ViewValueFunctor<device_type, value_type>>;
    using record_type =
        Kokkos::Impl::SharedAllocationRecord<memory_space, functor_type>;

    // Query the mapping for byte-size of allocation.
    // If padding is allowed then pass in sizeof value type
    // for padding computation.
    using padding = std::integral_constant<
        unsigned int, alloc_prop::allow_padding ? sizeof(value_type) : 0>;

    m_impl_offset = offset_type(padding(), arg_layout);

    const size_t alloc_size =
        (m_impl_offset.span() * MemorySpanSize + MemorySpanMask) &
        ~size_t(MemorySpanMask);
    const std::string& alloc_name =
        Impl::get_property<Impl::LabelTag>(arg_prop);
    const execution_space& exec_space =
        Impl::get_property<Impl::ExecutionSpaceTag>(arg_prop);
    const memory_space& mem_space =
        Impl::get_property<Impl::MemorySpaceTag>(arg_prop);

    // Create shared memory tracking record with allocate memory from the memory
    // space
    record_type* const record =
        execution_space_specified
            ? record_type::allocate(exec_space, mem_space, alloc_name,
                                    alloc_size)
            : record_type::allocate(mem_space, alloc_name, alloc_size);

    m_impl_handle = handle_type(reinterpret_cast<pointer_type>(record->data()));

    functor_type functor =
        execution_space_specified
            ? functor_type(exec_space, (value_type*)m_impl_handle,
                           m_impl_offset.span(), alloc_name)
            : functor_type((value_type*)m_impl_handle, m_impl_offset.span(),
                           alloc_name);

    //  Only initialize if the allocation is non-zero.
    //  May be zero if one of the dimensions is zero.
    if constexpr (alloc_prop::initialize)
      if (alloc_size) {
        // Assume destruction is only required when construction is requested.
        // The ViewValueFunctor has both value construction and destruction
        // operators.
        record->m_destroy = std::move(functor);

        // Construct values
        record->m_destroy.construct_shared_allocation();
      }

    return record;
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
/** \brief  Assign compatible default mappings */

template <class DstTraits, class SrcTraits>
class ViewMapping<
    DstTraits, SrcTraits,
    std::enable_if_t<(
        !(std::is_same<typename SrcTraits::array_layout, LayoutStride>::
              value) &&  // Added to have a new specialization for SrcType of
                         // LayoutStride
        // default mappings
        std::is_void<typename DstTraits::specialize>::value &&
        std::is_void<typename SrcTraits::specialize>::value &&
        (
            // same layout
            std::is_same<typename DstTraits::array_layout,
                         typename SrcTraits::array_layout>::value ||
            // known layout
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
                           Kokkos::LayoutStride>::value))))>> {
 private:
  enum {
    is_assignable_space = Kokkos::Impl::MemorySpaceAccess<
        typename DstTraits::memory_space,
        typename SrcTraits::memory_space>::assignable
  };

  enum {
    is_assignable_value_type =
        std::is_same<typename DstTraits::value_type,
                     typename SrcTraits::value_type>::value ||
        std::is_same<typename DstTraits::value_type,
                     typename SrcTraits::const_value_type>::value
  };

  enum {
    is_assignable_dimension =
        ViewDimensionAssignable<typename DstTraits::dimension,
                                typename SrcTraits::dimension>::value
  };

  enum {
    is_assignable_layout =
        std::is_same<typename DstTraits::array_layout,
                     typename SrcTraits::array_layout>::value ||
        std::is_same<typename DstTraits::array_layout,
                     Kokkos::LayoutStride>::value ||
        (DstTraits::dimension::rank == 0) || (DstTraits::dimension::rank == 1)
  };

 public:
  enum {
    is_assignable_data_type =
        is_assignable_value_type && is_assignable_dimension
  };
  enum {
    is_assignable = is_assignable_space && is_assignable_value_type &&
                    is_assignable_dimension && is_assignable_layout
  };

  using TrackType = Kokkos::Impl::SharedAllocationTracker;
  using DstType   = ViewMapping<DstTraits, void>;
  using SrcType   = ViewMapping<SrcTraits, void>;

  KOKKOS_INLINE_FUNCTION
  static void assign(DstType& dst, const SrcType& src,
                     const TrackType& src_track) {
    static_assert(is_assignable_space,
                  "View assignment must have compatible spaces");

    static_assert(
        is_assignable_value_type,
        "View assignment must have same value type or const = non-const");

    static_assert(is_assignable_dimension,
                  "View assignment must have compatible dimensions");

    static_assert(
        is_assignable_layout,
        "View assignment must have compatible layout or have rank <= 1");

    using dst_offset_type = typename DstType::offset_type;

    if (size_t(DstTraits::dimension::rank_dynamic) <
        size_t(SrcTraits::dimension::rank_dynamic)) {
      using dst_dim   = typename DstTraits::dimension;
      bool assignable = ((1 > DstTraits::dimension::rank_dynamic &&
                          1 <= SrcTraits::dimension::rank_dynamic)
                             ? dst_dim::ArgN0 == src.dimension_0()
                             : true) &&
                        ((2 > DstTraits::dimension::rank_dynamic &&
                          2 <= SrcTraits::dimension::rank_dynamic)
                             ? dst_dim::ArgN1 == src.dimension_1()
                             : true) &&
                        ((3 > DstTraits::dimension::rank_dynamic &&
                          3 <= SrcTraits::dimension::rank_dynamic)
                             ? dst_dim::ArgN2 == src.dimension_2()
                             : true) &&
                        ((4 > DstTraits::dimension::rank_dynamic &&
                          4 <= SrcTraits::dimension::rank_dynamic)
                             ? dst_dim::ArgN3 == src.dimension_3()
                             : true) &&
                        ((5 > DstTraits::dimension::rank_dynamic &&
                          5 <= SrcTraits::dimension::rank_dynamic)
                             ? dst_dim::ArgN4 == src.dimension_4()
                             : true) &&
                        ((6 > DstTraits::dimension::rank_dynamic &&
                          6 <= SrcTraits::dimension::rank_dynamic)
                             ? dst_dim::ArgN5 == src.dimension_5()
                             : true) &&
                        ((7 > DstTraits::dimension::rank_dynamic &&
                          7 <= SrcTraits::dimension::rank_dynamic)
                             ? dst_dim::ArgN6 == src.dimension_6()
                             : true) &&
                        ((8 > DstTraits::dimension::rank_dynamic &&
                          8 <= SrcTraits::dimension::rank_dynamic)
                             ? dst_dim::ArgN7 == src.dimension_7()
                             : true);
      if (!assignable)
        Kokkos::abort(
            "View Assignment: trying to assign runtime dimension to non "
            "matching compile time dimension.");
    }
    dst.m_impl_offset = dst_offset_type(src.m_impl_offset);
    dst.m_impl_handle = Kokkos::Impl::ViewDataHandle<DstTraits>::assign(
        src.m_impl_handle, src_track);
  }
};

//----------------------------------------------------------------------------
// Create new specialization for SrcType of LayoutStride. Runtime check for
// compatible layout
template <class DstTraits, class SrcTraits>
class ViewMapping<
    DstTraits, SrcTraits,
    std::enable_if_t<(
        std::is_same<typename SrcTraits::array_layout,
                     Kokkos::LayoutStride>::value &&
        std::is_void<typename DstTraits::specialize>::value &&
        std::is_void<typename SrcTraits::specialize>::value &&
        (
            // same layout
            std::is_same<typename DstTraits::array_layout,
                         typename SrcTraits::array_layout>::value ||
            // known layout
            (std::is_same<typename DstTraits::array_layout,
                          Kokkos::LayoutLeft>::value ||
             std::is_same<typename DstTraits::array_layout,
                          Kokkos::LayoutRight>::value ||
             std::is_same<typename DstTraits::array_layout,
                          Kokkos::LayoutStride>::value)))>> {
 private:
  enum {
    is_assignable_space = Kokkos::Impl::MemorySpaceAccess<
        typename DstTraits::memory_space,
        typename SrcTraits::memory_space>::assignable
  };

  enum {
    is_assignable_value_type =
        std::is_same<typename DstTraits::value_type,
                     typename SrcTraits::value_type>::value ||
        std::is_same<typename DstTraits::value_type,
                     typename SrcTraits::const_value_type>::value
  };

  enum {
    is_assignable_dimension =
        ViewDimensionAssignable<typename DstTraits::dimension,
                                typename SrcTraits::dimension>::value
  };

 public:
  enum {
    is_assignable_data_type =
        is_assignable_value_type && is_assignable_dimension
  };
  enum {
    is_assignable = is_assignable_space && is_assignable_value_type &&
                    is_assignable_dimension
  };

  using TrackType = Kokkos::Impl::SharedAllocationTracker;
  using DstType   = ViewMapping<DstTraits, void>;
  using SrcType   = ViewMapping<SrcTraits, void>;

  KOKKOS_INLINE_FUNCTION
  static bool assignable_layout_check(DstType&,
                                      const SrcType& src)  // Runtime check
  {
    size_t strides[9];
    bool assignable = true;
    src.stride(strides);
    size_t exp_stride = 1;
    if (std::is_same<typename DstTraits::array_layout,
                     Kokkos::LayoutLeft>::value) {
      for (int i = 0; i < (int)src.Rank; i++) {
        if (i > 0) exp_stride *= src.extent(i - 1);
        if (strides[i] != exp_stride) {
          assignable = false;
          break;
        }
      }
    } else if (std::is_same<typename DstTraits::array_layout,
                            Kokkos::LayoutRight>::value) {
      for (int i = 0; i < (int)src.Rank; i++) {
        if (i > 0) exp_stride *= src.extent(src.Rank - i);
        if (strides[src.Rank - 1 - i] != exp_stride) {
          assignable = false;
          break;
        }
      }
    }
    return assignable;
  }

  KOKKOS_INLINE_FUNCTION
  static void assign(DstType& dst, const SrcType& src,
                     const TrackType& src_track) {
    static_assert(is_assignable_space,
                  "View assignment must have compatible spaces");

    static_assert(
        is_assignable_value_type,
        "View assignment must have same value type or const = non-const");

    static_assert(is_assignable_dimension,
                  "View assignment must have compatible dimensions");

    bool assignable_layout = assignable_layout_check(dst, src);  // Runtime
                                                                 // check
    if (!assignable_layout)
      Kokkos::abort("View assignment must have compatible layouts\n");

    using dst_offset_type = typename DstType::offset_type;

    if (size_t(DstTraits::dimension::rank_dynamic) <
        size_t(SrcTraits::dimension::rank_dynamic)) {
      using dst_dim   = typename DstTraits::dimension;
      bool assignable = ((1 > DstTraits::dimension::rank_dynamic &&
                          1 <= SrcTraits::dimension::rank_dynamic)
                             ? dst_dim::ArgN0 == src.dimension_0()
                             : true) &&
                        ((2 > DstTraits::dimension::rank_dynamic &&
                          2 <= SrcTraits::dimension::rank_dynamic)
                             ? dst_dim::ArgN1 == src.dimension_1()
                             : true) &&
                        ((3 > DstTraits::dimension::rank_dynamic &&
                          3 <= SrcTraits::dimension::rank_dynamic)
                             ? dst_dim::ArgN2 == src.dimension_2()
                             : true) &&
                        ((4 > DstTraits::dimension::rank_dynamic &&
                          4 <= SrcTraits::dimension::rank_dynamic)
                             ? dst_dim::ArgN3 == src.dimension_3()
                             : true) &&
                        ((5 > DstTraits::dimension::rank_dynamic &&
                          5 <= SrcTraits::dimension::rank_dynamic)
                             ? dst_dim::ArgN4 == src.dimension_4()
                             : true) &&
                        ((6 > DstTraits::dimension::rank_dynamic &&
                          6 <= SrcTraits::dimension::rank_dynamic)
                             ? dst_dim::ArgN5 == src.dimension_5()
                             : true) &&
                        ((7 > DstTraits::dimension::rank_dynamic &&
                          7 <= SrcTraits::dimension::rank_dynamic)
                             ? dst_dim::ArgN6 == src.dimension_6()
                             : true) &&
                        ((8 > DstTraits::dimension::rank_dynamic &&
                          8 <= SrcTraits::dimension::rank_dynamic)
                             ? dst_dim::ArgN7 == src.dimension_7()
                             : true);
      if (!assignable)
        Kokkos::abort(
            "View Assignment: trying to assign runtime dimension to non "
            "matching compile time dimension.");
    }
    dst.m_impl_offset = dst_offset_type(src.m_impl_offset);
    dst.m_impl_handle = Kokkos::Impl::ViewDataHandle<DstTraits>::assign(
        src.m_impl_handle, src_track);
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// Subview mapping.
// Deduce destination view type from source view traits and subview arguments

template <class, class ValueType, class Exts, class... Args>
struct SubViewDataTypeImpl;

/* base case */
template <class ValueType>
struct SubViewDataTypeImpl<void, ValueType, Kokkos::Experimental::Extents<>> {
  using type = ValueType;
};

/* for integral args, subview doesn't have that dimension */
template <class ValueType, size_t Ext, size_t... Exts, class Integral,
          class... Args>
struct SubViewDataTypeImpl<
    std::enable_if_t<std::is_integral<std::decay_t<Integral>>::value>,
    ValueType, Kokkos::Experimental::Extents<Ext, Exts...>, Integral, Args...>
    : SubViewDataTypeImpl<void, ValueType,
                          Kokkos::Experimental::Extents<Exts...>, Args...> {};

/* for ALL slice, subview has the same dimension */
template <class ValueType, size_t Ext, size_t... Exts, class... Args>
struct SubViewDataTypeImpl<void, ValueType,
                           Kokkos::Experimental::Extents<Ext, Exts...>,
                           Kokkos::ALL_t, Args...>
    : SubViewDataTypeImpl<void, typename ApplyExtent<ValueType, Ext>::type,
                          Kokkos::Experimental::Extents<Exts...>, Args...> {};

/* for pair-style slice, subview has dynamic dimension, since pair doesn't give
 * static sizes */
/* Since we don't allow interleaving of dynamic and static extents, make all of
 * the dimensions to the left dynamic  */
template <class ValueType, size_t Ext, size_t... Exts, class PairLike,
          class... Args>
struct SubViewDataTypeImpl<
    std::enable_if_t<is_pair_like<PairLike>::value>, ValueType,
    Kokkos::Experimental::Extents<Ext, Exts...>, PairLike, Args...>
    : SubViewDataTypeImpl<
          void, typename make_all_extents_into_pointers<ValueType>::type*,
          Kokkos::Experimental::Extents<Exts...>, Args...> {};

template <class ValueType, class Exts, class... Args>
struct SubViewDataType : SubViewDataTypeImpl<void, ValueType, Exts, Args...> {};

//----------------------------------------------------------------------------

template <class SrcTraits, class... Args>
class ViewMapping<
    std::enable_if_t<(std::is_void<typename SrcTraits::specialize>::value &&
                      (std::is_same<typename SrcTraits::array_layout,
                                    Kokkos::LayoutLeft>::value ||
                       std::is_same<typename SrcTraits::array_layout,
                                    Kokkos::LayoutRight>::value ||
                       std::is_same<typename SrcTraits::array_layout,
                                    Kokkos::LayoutStride>::value))>,
    SrcTraits, Args...> {
 private:
  static_assert(SrcTraits::rank == sizeof...(Args),
                "Subview mapping requires one argument for each dimension of "
                "source View");

  enum {
    RZ = false,
    R0 = bool(is_integral_extent<0, Args...>::value),
    R1 = bool(is_integral_extent<1, Args...>::value),
    R2 = bool(is_integral_extent<2, Args...>::value),
    R3 = bool(is_integral_extent<3, Args...>::value),
    R4 = bool(is_integral_extent<4, Args...>::value),
    R5 = bool(is_integral_extent<5, Args...>::value),
    R6 = bool(is_integral_extent<6, Args...>::value),
    R7 = bool(is_integral_extent<7, Args...>::value)
  };

  enum {
    rank = unsigned(R0) + unsigned(R1) + unsigned(R2) + unsigned(R3) +
           unsigned(R4) + unsigned(R5) + unsigned(R6) + unsigned(R7)
  };

  // Whether right-most rank is a range.
  enum {
    R0_rev =
        (0 == SrcTraits::rank
             ? RZ
             : (1 == SrcTraits::rank
                    ? R0
                    : (2 == SrcTraits::rank
                           ? R1
                           : (3 == SrcTraits::rank
                                  ? R2
                                  : (4 == SrcTraits::rank
                                         ? R3
                                         : (5 == SrcTraits::rank
                                                ? R4
                                                : (6 == SrcTraits::rank
                                                       ? R5
                                                       : (7 == SrcTraits::rank
                                                              ? R6
                                                              : R7))))))))
  };

  // Subview's layout
  using array_layout = std::conditional_t<
      (            /* Same array layout IF */
       (rank == 0) /* output rank zero */
       || SubviewLegalArgsCompileTime<typename SrcTraits::array_layout,
                                      typename SrcTraits::array_layout, rank,
                                      SrcTraits::rank, 0, Args...>::value ||
       // OutputRank 1 or 2, InputLayout Left, Interval 0
       // because single stride one or second index has a stride.
       (rank <= 2 && R0 &&
        std::is_same<typename SrcTraits::array_layout,
                     Kokkos::LayoutLeft>::value)  // replace with input rank
       ||
       // OutputRank 1 or 2, InputLayout Right, Interval [InputRank-1]
       // because single stride one or second index has a stride.
       (rank <= 2 && R0_rev &&
        std::is_same<typename SrcTraits::array_layout,
                     Kokkos::LayoutRight>::value)  // replace input rank
       ),
      typename SrcTraits::array_layout, Kokkos::LayoutStride>;

  using value_type = typename SrcTraits::value_type;

  using data_type =
      typename SubViewDataType<value_type,
                               typename Kokkos::Impl::ParseViewExtents<
                                   typename SrcTraits::data_type>::type,
                               Args...>::type;

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

  // The presumed type is 'ViewMapping< traits_type , void >'
  // However, a compatible ViewMapping is acceptable.
  template <class DstTraits>
  KOKKOS_INLINE_FUNCTION static void assign(
      ViewMapping<DstTraits, void>& dst,
      ViewMapping<SrcTraits, void> const& src, Args... args) {
    static_assert(ViewMapping<DstTraits, traits_type, void>::is_assignable,
                  "Subview destination type must be compatible with subview "
                  "derived type");

    using DstType = ViewMapping<DstTraits, void>;

    using dst_offset_type = typename DstType::offset_type;

    const SubviewExtents<SrcTraits::rank, rank> extents(src.m_impl_offset.m_dim,
                                                        args...);

    dst.m_impl_offset = dst_offset_type(src.m_impl_offset, extents);

    dst.m_impl_handle = ViewDataHandle<DstTraits>::assign(
        src.m_impl_handle,
        src.m_impl_offset(extents.domain_offset(0), extents.domain_offset(1),
                          extents.domain_offset(2), extents.domain_offset(3),
                          extents.domain_offset(4), extents.domain_offset(5),
                          extents.domain_offset(6), extents.domain_offset(7)));
  }
};

//----------------------------------------------------------------------------

}  // namespace Impl
}  // namespace Kokkos

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Impl {

template <class Map, class... Indices, std::size_t... Enumerate>
KOKKOS_FUNCTION bool within_range(Map const& map,
                                  std::index_sequence<Enumerate...>,
                                  Indices... indices) {
  return (((std::size_t)indices < map.extent(Enumerate)) && ...);
}

template <class... Indices>
KOKKOS_FUNCTION constexpr char* append_formatted_multidimensional_index(
    char* dest, Indices... indices) {
  char* d = dest;
  strcat(d, "[");
  (
      [&] {
        d += strlen(d);
        to_chars_i(d,
                   d + 20,  // 20 digits ought to be enough
                   indices);
        strcat(d, ",");
      }(),
      ...);
  d[strlen(d) - 1] = ']';  // overwrite trailing comma
  return dest;
}

template <class Map, class... Indices, std::size_t... Enumerate>
KOKKOS_FUNCTION void print_extents(char* dest, Map const& map,
                                   std::index_sequence<Enumerate...>) {
  append_formatted_multidimensional_index(dest, map.extent(Enumerate)...);
}

template <class T>
using printable_label_typedef_t = typename T::printable_label_typedef;

template <class MemorySpace, class ViewType, class MapType, class... Args>
KOKKOS_INLINE_FUNCTION void view_verify_operator_bounds(
    Kokkos::Impl::ViewTracker<ViewType> const& tracker, const MapType& map,
    Args... args) {
  if (!within_range(map, std::make_index_sequence<sizeof...(Args)>(),
                    args...)) {
    char err[256] = "";
    strcat(err, "Kokkos::View ERROR: out of bounds access");
    strcat(err, " label=(\"");
    KOKKOS_IF_ON_HOST(
        if (tracker.m_tracker.has_record()) {
          strncat(err, tracker.m_tracker.template get_label<void>().c_str(),
                  128);
        } else { strcat(err, "**UNMANAGED**"); })
    KOKKOS_IF_ON_DEVICE([&] {
      // Check #1: is there a SharedAllocationRecord?  (we won't use it, but
      // if its not there then there isn't a corresponding
      // SharedAllocationHeader containing a label).  This check should cover
      // the case of Views that don't have the Unmanaged trait but were
      // initialized by pointer.
      if (!tracker.m_tracker.has_record()) {
        strcat(err, "**UNMANAGED**");
        return;
      }
      // Check #2: does the ViewMapping have the printable_label_typedef
      // defined? See above that only the non-specialized standard-layout
      // ViewMapping has this defined by default. The existence of this
      // alias indicates the existence of MapType::is_managed
      if constexpr (is_detected_v<printable_label_typedef_t, MapType>) {
        // Check #3: is the View managed as determined by the MemoryTraits?
        if constexpr (MapType::is_managed != 0) {
          SharedAllocationHeader const* const header =
              SharedAllocationHeader::get_header(
                  static_cast<void const*>(map.data()));
          char const* const label = header->label();
          strcat(err, label);
          return;
        }
        strcat(err, "**UNAVAILABLE**");
      }
    }();)
    strcat(err, "\") with indices ");
    append_formatted_multidimensional_index(err, args...);
    strcat(err, " but extents ");
    print_extents(err, map, std::make_index_sequence<sizeof...(Args)>());
    Kokkos::abort(err);
  }
}

// primary template: memory space is accessible, do nothing.
template <class MemorySpace, class AccessSpace,
          bool = SpaceAccessibility<AccessSpace, MemorySpace>::accessible>
struct RuntimeCheckViewMemoryAccessViolation {
  template <class Track, class Map>
  KOKKOS_FUNCTION RuntimeCheckViewMemoryAccessViolation(char const* const,
                                                        Track const&,
                                                        Map const&) {}
};

// explicit specialization: memory access violation will occur, call abort with
// the specified error message.
template <class MemorySpace, class AccessSpace>
struct RuntimeCheckViewMemoryAccessViolation<MemorySpace, AccessSpace, false> {
  template <class Track, class Map>
  KOKKOS_FUNCTION RuntimeCheckViewMemoryAccessViolation(char const* const msg,
                                                        Track const& track,
                                                        Map const&) {
    char err[256] = "";
    strncat(err, msg, 64);
    strcat(err, " (label=\"");

    KOKKOS_IF_ON_HOST(({
      auto const tracker = track.m_tracker;

      if (tracker.has_record()) {
        strncat(err, tracker.template get_label<void>().c_str(), 128);
      } else {
        strcat(err, "**UNMANAGED**");
      }
    }))

    KOKKOS_IF_ON_DEVICE(({
      strcat(err, "**UNAVAILABLE**");
      (void)track;
    }))

    strcat(err, "\")");

    Kokkos::abort(err);
  }
};

template <class MemorySpace, class Track, class Map, class... Ignore>
KOKKOS_FUNCTION void runtime_check_memory_access_violation(
    char const* const msg, Track const& track, Map const& map, Ignore...) {
  KOKKOS_IF_ON_HOST(
      ((void)RuntimeCheckViewMemoryAccessViolation<MemorySpace,
                                                   DefaultHostExecutionSpace>(
           msg, track, map);))
  KOKKOS_IF_ON_DEVICE(
      ((void)RuntimeCheckViewMemoryAccessViolation<MemorySpace,
                                                   DefaultExecutionSpace>(
           msg, track, map);))
}

} /* namespace Impl */
} /* namespace Kokkos */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_EXPERIMENTAL_VIEW_MAPPING_HPP */
