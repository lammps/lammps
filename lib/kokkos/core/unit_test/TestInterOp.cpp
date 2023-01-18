/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

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
