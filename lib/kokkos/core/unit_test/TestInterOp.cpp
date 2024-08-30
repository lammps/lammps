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

#include <Kokkos_Core.hpp>
#include <Kokkos_DynRankView.hpp>
#include <KokkosExp_InterOp.hpp>

// View
static_assert(
    std::is_same<
        Kokkos::Experimental::python_view_type_t<Kokkos::View<double*>>,
        Kokkos::View<double*,
                     typename Kokkos::DefaultExecutionSpace::array_layout,
                     typename Kokkos::DefaultExecutionSpace::memory_space,
                     Kokkos::Experimental::DefaultViewHooks>>::value,
    "Error! Unexpected python_view_type for: View");

// DynRankView
static_assert(
    std::is_same<
        Kokkos::Experimental::python_view_type_t<Kokkos::DynRankView<double>>,
        Kokkos::DynRankView<
            double, typename Kokkos::DefaultExecutionSpace::array_layout,
            typename Kokkos::DefaultExecutionSpace::memory_space>>::value,
    "Error! Unexpected python_view_type for: DynRankView");

// View + Execution Space
static_assert(
    std::is_same<
        Kokkos::Experimental::python_view_type_t<
            Kokkos::View<double*, Kokkos::DefaultExecutionSpace>>,
        Kokkos::View<double*,
                     typename Kokkos::DefaultExecutionSpace::array_layout,
                     typename Kokkos::DefaultExecutionSpace::memory_space,
                     Kokkos::Experimental::DefaultViewHooks>>::value,
    "Error! Unexpected python_view_type for: View + Execution Space");

// DynRankView + Execution Space
static_assert(
    std::is_same<
        Kokkos::Experimental::python_view_type_t<
            Kokkos::DynRankView<double, Kokkos::DefaultExecutionSpace>>,
        Kokkos::DynRankView<
            double, typename Kokkos::DefaultExecutionSpace::array_layout,
            typename Kokkos::DefaultExecutionSpace::memory_space>>::value,
    "Error! Unexpected python_view_type for: DynRankView + Execution Space");

// View + Memory space
static_assert(
    std::is_same<Kokkos::Experimental::python_view_type_t<
                     Kokkos::View<int64_t*, Kokkos::HostSpace>>,
                 Kokkos::View<int64_t*, Kokkos::LayoutRight, Kokkos::HostSpace,
                              Kokkos::Experimental::DefaultViewHooks>>::value,
    "Error! Unexpected python_view_type for: View + Memory space");

// DynRankView + Memory space
static_assert(
    std::is_same<Kokkos::Experimental::python_view_type_t<
                     Kokkos::DynRankView<int16_t, Kokkos::HostSpace>>,
                 Kokkos::DynRankView<int16_t, Kokkos::LayoutRight,
                                     Kokkos::HostSpace>>::value,
    "Error! Unexpected python_view_type for: DynRankView + Memory space");

// View + Layout + Execution space
static_assert(
    std::is_same<
        Kokkos::Experimental::python_view_type_t<Kokkos::View<
            int**, Kokkos::LayoutLeft, Kokkos::DefaultExecutionSpace>>,
        Kokkos::View<int**, Kokkos::LayoutLeft,
                     typename Kokkos::DefaultExecutionSpace::memory_space,
                     Kokkos::Experimental::DefaultViewHooks>>::value,
    "Error! Unexpected python_view_type for: View + Layout + Execution space");

// DynRankView + Layout + Execution space
static_assert(
    std::is_same<Kokkos::Experimental::python_view_type_t<Kokkos::DynRankView<
                     int, Kokkos::LayoutLeft, Kokkos::DefaultExecutionSpace>>,
                 Kokkos::DynRankView<int, Kokkos::LayoutLeft,
                                     typename Kokkos::DefaultExecutionSpace::
                                         memory_space>>::value,
    "Error! Unexpected python_view_type for: DynRankView + Layout + Execution "
    "space");

// View + Layout + Memory Space
static_assert(
    std::is_same<Kokkos::Experimental::python_view_type_t<Kokkos::View<
                     uint32_t**, Kokkos::LayoutLeft, Kokkos::HostSpace>>,
                 Kokkos::View<uint32_t**, Kokkos::LayoutLeft, Kokkos::HostSpace,
                              Kokkos::Experimental::DefaultViewHooks>>::value,
    "Error! Unexpected python_view_type for: View + Layout + Memory Space");

// DynRankView + Layout + Memory Space
static_assert(
    std::is_same<Kokkos::Experimental::python_view_type_t<Kokkos::DynRankView<
                     uint64_t, Kokkos::LayoutLeft, Kokkos::HostSpace>>,
                 Kokkos::DynRankView<uint64_t, Kokkos::LayoutLeft,
                                     Kokkos::HostSpace>>::value,
    "Error! Unexpected python_view_type for: DynRankView + Layout + Memory "
    "Space");

// View + Layout + Execution space + Memory Trait
static_assert(
    std::is_same<
        Kokkos::Experimental::python_view_type_t<Kokkos::View<
            float***, Kokkos::LayoutLeft, Kokkos::DefaultHostExecutionSpace,
            Kokkos::MemoryTraits<Kokkos::RandomAccess>>>,
        Kokkos::View<float***, Kokkos::LayoutLeft,
                     typename Kokkos::DefaultHostExecutionSpace::memory_space,
                     Kokkos::Experimental::DefaultViewHooks,
                     Kokkos::MemoryTraits<Kokkos::RandomAccess>>>::value,
    "Error! Unexpected python_view_type for: View + Layout + Execution space + "
    "Memory Trait");

// DynRankView + Layout + Execution space  + Memory trait
static_assert(
    std::is_same<
        Kokkos::Experimental::python_view_type_t<Kokkos::DynRankView<
            float, Kokkos::LayoutLeft, Kokkos::DefaultHostExecutionSpace,
            Kokkos::MemoryTraits<Kokkos::Atomic>>>,
        Kokkos::DynRankView<
            float, Kokkos::LayoutLeft,
            typename Kokkos::DefaultHostExecutionSpace::memory_space,
            Kokkos::MemoryTraits<Kokkos::Atomic>>>::value,
    "Error! Unexpected python_view_type for: DynRankView + Layout + Execution "
    "space  + Memory trait");
