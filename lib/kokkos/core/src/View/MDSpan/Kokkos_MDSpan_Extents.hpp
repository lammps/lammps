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

#ifndef KOKKOS_EXPERIMENTAL_MDSPAN_EXTENTS_HPP
#define KOKKOS_EXPERIMENTAL_MDSPAN_EXTENTS_HPP

#include "Kokkos_MDSpan_Header.hpp"

namespace Kokkos::Impl {

// Forward declarations from impl/Kokkos_ViewMapping.hpp
// We cannot include directly since ViewMapping is used elsewhere in View.
// After View is fully moved to mdspan we can include it only from here.
template <class DataType>
struct ViewArrayAnalysis;

template <std::size_t... Vals>
struct ViewDimension;

template <class T, class Dim>
struct ViewDataType;
}  // namespace Kokkos::Impl

namespace Kokkos::Experimental::Impl {

// A few things to note --
// - mdspan allows for 0-rank extents similarly to View, so we don't need
// special handling of this case
// - View dynamic dimensions must be appear before static dimensions. This isn't
// a requirement in mdspan but won't cause an issue here
template <std::size_t N>
struct ExtentFromDimension {
  static constexpr std::size_t value = N;
};

// Kokkos uses a dimension of '0' to denote a dynamic dimension.
template <>
struct ExtentFromDimension<std::size_t{0}> {
  static constexpr std::size_t value = dynamic_extent;
};

template <std::size_t N>
struct DimensionFromExtent {
  static constexpr std::size_t value = N;
};

template <>
struct DimensionFromExtent<dynamic_extent> {
  static constexpr std::size_t value = std::size_t{0};
};

template <class IndexType, class Dimension, class Indices>
struct ExtentsFromDimension;

template <class IndexType, class Dimension, std::size_t... Indices>
struct ExtentsFromDimension<IndexType, Dimension,
                            std::index_sequence<Indices...>> {
  using type =
      extents<IndexType,
              ExtentFromDimension<Dimension::static_extent(Indices)>::value...>;
};

template <class Extents, class Indices>
struct DimensionsFromExtent;

template <class Extents, std::size_t... Indices>
struct DimensionsFromExtent<Extents, std::index_sequence<Indices...>> {
  using type = ::Kokkos::Impl::ViewDimension<
      DimensionFromExtent<Extents::static_extent(Indices)>::value...>;
};

template <class IndexType, class DataType>
struct ExtentsFromDataType {
  using array_analysis = ::Kokkos::Impl::ViewArrayAnalysis<DataType>;
  using dimension_type = typename array_analysis::dimension;

  using type = typename ExtentsFromDimension<
      IndexType, dimension_type,
      std::make_index_sequence<dimension_type::rank>>::type;
};

template <class T, class Extents>
struct DataTypeFromExtents {
  using extents_type   = Extents;
  using dimension_type = typename DimensionsFromExtent<
      Extents, std::make_index_sequence<extents_type::rank()>>::type;

  // Will cause a compile error if it is malformed (i.e. dynamic after static)
  using type = typename ::Kokkos::Impl::ViewDataType<T, dimension_type>::type;
};
}  // namespace Kokkos::Experimental::Impl

#endif  // KOKKOS_EXPERIMENTAL_MDSPAN_EXTENTS_HPP
