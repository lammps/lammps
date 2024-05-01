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

#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
static_assert(false,
              "Including non-public Kokkos header files is not allowed.");
#endif

#ifndef KOKKOS_VIEW_DATA_ANALYSIS_HPP
#define KOKKOS_VIEW_DATA_ANALYSIS_HPP

#include <Kokkos_Macros.hpp>

namespace Kokkos::Impl {

template <unsigned I, size_t... Args>
struct variadic_size_t {
  enum : size_t { value = KOKKOS_INVALID_INDEX };
};

template <size_t Val, size_t... Args>
struct variadic_size_t<0, Val, Args...> {
  enum : size_t { value = Val };
};

template <unsigned I, size_t Val, size_t... Args>
struct variadic_size_t<I, Val, Args...> {
  enum : size_t { value = variadic_size_t<I - 1, Args...>::value };
};

template <size_t... Args>
struct rank_dynamic;

template <>
struct rank_dynamic<> {
  enum : unsigned { value = 0 };
};

template <size_t Val, size_t... Args>
struct rank_dynamic<Val, Args...> {
  enum : unsigned { value = (Val == 0 ? 1 : 0) + rank_dynamic<Args...>::value };
};

#define KOKKOS_IMPL_VIEW_DIMENSION(R)                                       \
  template <size_t V, unsigned>                                             \
  struct ViewDimension##R {                                                 \
    static constexpr size_t ArgN##R = (V != KOKKOS_INVALID_INDEX ? V : 1);  \
    static constexpr size_t N##R    = (V != KOKKOS_INVALID_INDEX ? V : 1);  \
    KOKKOS_INLINE_FUNCTION explicit ViewDimension##R(size_t) {}             \
    ViewDimension##R()                        = default;                    \
    ViewDimension##R(const ViewDimension##R&) = default;                    \
    ViewDimension##R& operator=(const ViewDimension##R&) = default;         \
  };                                                                        \
  template <size_t V, unsigned RD>                                          \
  constexpr size_t ViewDimension##R<V, RD>::ArgN##R;                        \
  template <size_t V, unsigned RD>                                          \
  constexpr size_t ViewDimension##R<V, RD>::N##R;                           \
  template <unsigned RD>                                                    \
  struct ViewDimension##R<0u, RD> {                                         \
    static constexpr size_t ArgN##R = 0;                                    \
    std::conditional_t<(RD < 3), size_t, unsigned> N##R;                    \
    ViewDimension##R()                        = default;                    \
    ViewDimension##R(const ViewDimension##R&) = default;                    \
    ViewDimension##R& operator=(const ViewDimension##R&) = default;         \
    KOKKOS_INLINE_FUNCTION explicit ViewDimension##R(size_t V) : N##R(V) {} \
  };                                                                        \
  template <unsigned RD>                                                    \
  constexpr size_t ViewDimension##R<0u, RD>::ArgN##R;

KOKKOS_IMPL_VIEW_DIMENSION(0)
KOKKOS_IMPL_VIEW_DIMENSION(1)
KOKKOS_IMPL_VIEW_DIMENSION(2)
KOKKOS_IMPL_VIEW_DIMENSION(3)
KOKKOS_IMPL_VIEW_DIMENSION(4)
KOKKOS_IMPL_VIEW_DIMENSION(5)
KOKKOS_IMPL_VIEW_DIMENSION(6)
KOKKOS_IMPL_VIEW_DIMENSION(7)

#undef KOKKOS_IMPL_VIEW_DIMENSION

// MSVC does not do empty base class optimization by default.
// Per standard it is required for standard layout types
template <size_t... Vals>
struct KOKKOS_IMPL_ENFORCE_EMPTY_BASE_OPTIMIZATION ViewDimension
    : public ViewDimension0<variadic_size_t<0u, Vals...>::value,
                            rank_dynamic<Vals...>::value>,
      public ViewDimension1<variadic_size_t<1u, Vals...>::value,
                            rank_dynamic<Vals...>::value>,
      public ViewDimension2<variadic_size_t<2u, Vals...>::value,
                            rank_dynamic<Vals...>::value>,
      public ViewDimension3<variadic_size_t<3u, Vals...>::value,
                            rank_dynamic<Vals...>::value>,
      public ViewDimension4<variadic_size_t<4u, Vals...>::value,
                            rank_dynamic<Vals...>::value>,
      public ViewDimension5<variadic_size_t<5u, Vals...>::value,
                            rank_dynamic<Vals...>::value>,
      public ViewDimension6<variadic_size_t<6u, Vals...>::value,
                            rank_dynamic<Vals...>::value>,
      public ViewDimension7<variadic_size_t<7u, Vals...>::value,
                            rank_dynamic<Vals...>::value> {
  using D0 = ViewDimension0<variadic_size_t<0U, Vals...>::value,
                            rank_dynamic<Vals...>::value>;
  using D1 = ViewDimension1<variadic_size_t<1U, Vals...>::value,
                            rank_dynamic<Vals...>::value>;
  using D2 = ViewDimension2<variadic_size_t<2U, Vals...>::value,
                            rank_dynamic<Vals...>::value>;
  using D3 = ViewDimension3<variadic_size_t<3U, Vals...>::value,
                            rank_dynamic<Vals...>::value>;
  using D4 = ViewDimension4<variadic_size_t<4U, Vals...>::value,
                            rank_dynamic<Vals...>::value>;
  using D5 = ViewDimension5<variadic_size_t<5U, Vals...>::value,
                            rank_dynamic<Vals...>::value>;
  using D6 = ViewDimension6<variadic_size_t<6U, Vals...>::value,
                            rank_dynamic<Vals...>::value>;
  using D7 = ViewDimension7<variadic_size_t<7U, Vals...>::value,
                            rank_dynamic<Vals...>::value>;

  using D0::ArgN0;
  using D1::ArgN1;
  using D2::ArgN2;
  using D3::ArgN3;
  using D4::ArgN4;
  using D5::ArgN5;
  using D6::ArgN6;
  using D7::ArgN7;

  using D0::N0;
  using D1::N1;
  using D2::N2;
  using D3::N3;
  using D4::N4;
  using D5::N5;
  using D6::N6;
  using D7::N7;

  static constexpr unsigned rank         = sizeof...(Vals);
  static constexpr unsigned rank_dynamic = Impl::rank_dynamic<Vals...>::value;

  ViewDimension()                     = default;
  ViewDimension(const ViewDimension&) = default;
  ViewDimension& operator=(const ViewDimension&) = default;

  KOKKOS_INLINE_FUNCTION
  constexpr ViewDimension(size_t n0, size_t n1, size_t n2, size_t n3, size_t n4,
                          size_t n5, size_t n6, size_t n7)
      : D0(n0 == KOKKOS_INVALID_INDEX ? 1 : n0),
        D1(n1 == KOKKOS_INVALID_INDEX ? 1 : n1),
        D2(n2 == KOKKOS_INVALID_INDEX ? 1 : n2),
        D3(n3 == KOKKOS_INVALID_INDEX ? 1 : n3),
        D4(n4 == KOKKOS_INVALID_INDEX ? 1 : n4),
        D5(n5 == KOKKOS_INVALID_INDEX ? 1 : n5),
        D6(n6 == KOKKOS_INVALID_INDEX ? 1 : n6),
        D7(n7 == KOKKOS_INVALID_INDEX ? 1 : n7) {}

  KOKKOS_INLINE_FUNCTION
  constexpr size_t extent(const unsigned r) const noexcept {
    return r == 0
               ? N0
               : (r == 1
                      ? N1
                      : (r == 2
                             ? N2
                             : (r == 3
                                    ? N3
                                    : (r == 4
                                           ? N4
                                           : (r == 5
                                                  ? N5
                                                  : (r == 6
                                                         ? N6
                                                         : (r == 7 ? N7
                                                                   : 0)))))));
  }

  static KOKKOS_INLINE_FUNCTION constexpr size_t static_extent(
      const unsigned r) noexcept {
    return r == 0
               ? ArgN0
               : (r == 1
                      ? ArgN1
                      : (r == 2
                             ? ArgN2
                             : (r == 3
                                    ? ArgN3
                                    : (r == 4
                                           ? ArgN4
                                           : (r == 5
                                                  ? ArgN5
                                                  : (r == 6
                                                         ? ArgN6
                                                         : (r == 7 ? ArgN7
                                                                   : 0)))))));
  }

  template <size_t N>
  struct prepend {
    using type = ViewDimension<N, Vals...>;
  };

  template <size_t N>
  struct append {
    using type = ViewDimension<Vals..., N>;
  };
};

template <class A, class B>
struct ViewDimensionJoin;

template <size_t... A, size_t... B>
struct ViewDimensionJoin<ViewDimension<A...>, ViewDimension<B...>> {
  using type = ViewDimension<A..., B...>;
};

//----------------------------------------------------------------------------

template <class DstDim, class SrcDim>
struct ViewDimensionAssignable;

template <size_t... DstArgs, size_t... SrcArgs>
struct ViewDimensionAssignable<ViewDimension<DstArgs...>,
                               ViewDimension<SrcArgs...>> {
  using dst = ViewDimension<DstArgs...>;
  using src = ViewDimension<SrcArgs...>;

  enum {
    value = unsigned(dst::rank) == unsigned(src::rank) &&
            (
                // Compile time check that potential static dimensions match
                ((1 > dst::rank_dynamic && 1 > src::rank_dynamic)
                     ? (size_t(dst::ArgN0) == size_t(src::ArgN0))
                     : true) &&
                ((2 > dst::rank_dynamic && 2 > src::rank_dynamic)
                     ? (size_t(dst::ArgN1) == size_t(src::ArgN1))
                     : true) &&
                ((3 > dst::rank_dynamic && 3 > src::rank_dynamic)
                     ? (size_t(dst::ArgN2) == size_t(src::ArgN2))
                     : true) &&
                ((4 > dst::rank_dynamic && 4 > src::rank_dynamic)
                     ? (size_t(dst::ArgN3) == size_t(src::ArgN3))
                     : true) &&
                ((5 > dst::rank_dynamic && 5 > src::rank_dynamic)
                     ? (size_t(dst::ArgN4) == size_t(src::ArgN4))
                     : true) &&
                ((6 > dst::rank_dynamic && 6 > src::rank_dynamic)
                     ? (size_t(dst::ArgN5) == size_t(src::ArgN5))
                     : true) &&
                ((7 > dst::rank_dynamic && 7 > src::rank_dynamic)
                     ? (size_t(dst::ArgN6) == size_t(src::ArgN6))
                     : true) &&
                ((8 > dst::rank_dynamic && 8 > src::rank_dynamic)
                     ? (size_t(dst::ArgN7) == size_t(src::ArgN7))
                     : true))
  };
};

/** \brief  Given a value type and dimension generate the View data type */
template <class T, class Dim>
struct ViewDataType;

template <class T>
struct ViewDataType<T, ViewDimension<>> {
  using type = T;
};

template <class T, size_t... Args>
struct ViewDataType<T, ViewDimension<0, Args...>> {
  using type = typename ViewDataType<T*, ViewDimension<Args...>>::type;
};

template <class T, size_t N, size_t... Args>
struct ViewDataType<T, ViewDimension<N, Args...>> {
  using type = typename ViewDataType<T, ViewDimension<Args...>>::type[N];
};

/**\brief  Analysis of View data type.
 *
 *  Data type conforms to one of the following patterns :
 *    {const} value_type [][#][#][#]
 *    {const} value_type ***[#][#][#]
 *  Where the sum of counts of '*' and '[#]' is at most ten.
 *
 *  Provide alias for ViewDimension<...> and value_type.
 */
template <class T>
struct ViewArrayAnalysis {
  using value_type           = T;
  using const_value_type     = std::add_const_t<T>;
  using non_const_value_type = std::remove_const_t<T>;
  using static_dimension     = ViewDimension<>;
  using dynamic_dimension    = ViewDimension<>;
  using dimension            = ViewDimension<>;
};

template <class T, size_t N>
struct ViewArrayAnalysis<T[N]> {
 private:
  using nested = ViewArrayAnalysis<T>;

 public:
  using value_type           = typename nested::value_type;
  using const_value_type     = typename nested::const_value_type;
  using non_const_value_type = typename nested::non_const_value_type;

  using static_dimension =
      typename nested::static_dimension::template prepend<N>::type;

  using dynamic_dimension = typename nested::dynamic_dimension;

  using dimension =
      typename ViewDimensionJoin<dynamic_dimension, static_dimension>::type;
};

template <class T>
struct ViewArrayAnalysis<T[]> {
 private:
  using nested           = ViewArrayAnalysis<T>;
  using nested_dimension = typename nested::dimension;

 public:
  using value_type           = typename nested::value_type;
  using const_value_type     = typename nested::const_value_type;
  using non_const_value_type = typename nested::non_const_value_type;

  using dynamic_dimension =
      typename nested::dynamic_dimension::template prepend<0>::type;

  using static_dimension = typename nested::static_dimension;

  using dimension =
      typename ViewDimensionJoin<dynamic_dimension, static_dimension>::type;
};

template <class T>
struct ViewArrayAnalysis<T*> {
 private:
  using nested = ViewArrayAnalysis<T>;

 public:
  using value_type           = typename nested::value_type;
  using const_value_type     = typename nested::const_value_type;
  using non_const_value_type = typename nested::non_const_value_type;

  using dynamic_dimension =
      typename nested::dynamic_dimension::template prepend<0>::type;

  using static_dimension = typename nested::static_dimension;

  using dimension =
      typename ViewDimensionJoin<dynamic_dimension, static_dimension>::type;
};

template <class DataType, class ArrayLayout, class ValueType>
struct ViewDataAnalysis {
 private:
  using array_analysis = ViewArrayAnalysis<DataType>;

  // ValueType is opportunity for partial specialization.
  // Must match array analysis when this default template is used.
  static_assert(
      std::is_same<ValueType,
                   typename array_analysis::non_const_value_type>::value);

 public:
  using specialize = void;  // No specialization

  using dimension            = typename array_analysis::dimension;
  using value_type           = typename array_analysis::value_type;
  using const_value_type     = typename array_analysis::const_value_type;
  using non_const_value_type = typename array_analysis::non_const_value_type;

  // Generate analogous multidimensional array specification type.
  using type       = typename ViewDataType<value_type, dimension>::type;
  using const_type = typename ViewDataType<const_value_type, dimension>::type;
  using non_const_type =
      typename ViewDataType<non_const_value_type, dimension>::type;

  // Generate "flattened" multidimensional array specification type.
  using scalar_array_type           = type;
  using const_scalar_array_type     = const_type;
  using non_const_scalar_array_type = non_const_type;
};

template <class Dimension, class Layout, class Enable = void>
struct ViewOffset {
  using is_mapping_plugin = std::false_type;
};
}  // namespace Kokkos::Impl

#endif  // KOKKOS_VIEW_DATA_ANALYSIS_HPP
