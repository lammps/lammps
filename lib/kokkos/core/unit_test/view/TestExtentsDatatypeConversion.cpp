// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
import kokkos.core_impl;
#else
#include <Kokkos_Core.hpp>
#endif
#include <type_traits>

#ifdef KOKKOS_ENABLE_IMPL_MDSPAN

namespace {

// Helper to make static tests more succinct
template <typename DataType, typename Extent>
constexpr bool datatype_matches_extent = std::is_same_v<
    typename Kokkos::Impl::ExtentsFromDataType<std::size_t, DataType>::type,
    Extent>;

template <typename DataType, typename BaseType, typename Extents>
constexpr bool extent_matches_datatype =
    std::is_same_v<DataType, typename Kokkos::Impl::DataTypeFromExtents<
                                 BaseType, Extents>::type>;

// Conversion from DataType to extents
// 0-rank view
static_assert(datatype_matches_extent<double, Kokkos::extents<std::size_t>>);

// Only dynamic
static_assert(datatype_matches_extent<
              double***,
              Kokkos::extents<std::size_t, Kokkos::dynamic_extent,
                              Kokkos::dynamic_extent, Kokkos::dynamic_extent>>);
// Only static
static_assert(
    datatype_matches_extent<double[2][3][17],
                            Kokkos::extents<std::size_t, std::size_t{2},
                                            std::size_t{3}, std::size_t{17}>>);

// Both dynamic and static
static_assert(datatype_matches_extent<
              double** [3][2][8],
              Kokkos::extents<std::size_t, Kokkos::dynamic_extent,
                              Kokkos::dynamic_extent, std::size_t { 3 },
                              std::size_t{2}, std::size_t{8}>>);

// Conversion from extents to DataType
// 0-rank extents
static_assert(
    extent_matches_datatype<double, double, Kokkos::extents<std::size_t>>);

// only dynamic
static_assert(extent_matches_datatype<
              double****, double,
              Kokkos::extents<std::size_t, Kokkos::dynamic_extent,
                              Kokkos::dynamic_extent, Kokkos::dynamic_extent,
                              Kokkos::dynamic_extent>>);

// only static
static_assert(extent_matches_datatype<double[7][5][3], double,
                                      Kokkos::extents<std::size_t, 7, 5, 3>>);

// both dynamic and static
static_assert(
    extent_matches_datatype<double*** [20][45], double,
                            Kokkos::extents<std::size_t, Kokkos::dynamic_extent,
                                            Kokkos::dynamic_extent,
                                            Kokkos::dynamic_extent, 20, 45>>);
}  // namespace

#endif  // KOKKOS_ENABLE_IMPL_MDSPAN
