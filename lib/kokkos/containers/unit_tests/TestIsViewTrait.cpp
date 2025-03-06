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
#include <Kokkos_DualView.hpp>
#include <Kokkos_DynRankView.hpp>
#include <Kokkos_DynamicView.hpp>
#include <Kokkos_OffsetView.hpp>
#include <Kokkos_ScatterView.hpp>

namespace {

using view_t          = Kokkos::View<int*>;
using dual_view_t     = Kokkos::DualView<int*>;
using dyn_rank_view_t = Kokkos::DynRankView<int*>;
using dynamic_view_t  = Kokkos::Experimental::DynamicView<int*>;
using offset_view_t   = Kokkos::Experimental::OffsetView<int*>;
using scatter_view_t  = Kokkos::Experimental::ScatterView<int*>;

static_assert(Kokkos::is_dual_view_v<dual_view_t>);
static_assert(!Kokkos::is_dyn_rank_view_v<dual_view_t>);
static_assert(!Kokkos::is_dynamic_view_v<dual_view_t>);
static_assert(!Kokkos::Experimental::is_offset_view_v<dual_view_t>);
static_assert(!Kokkos::Experimental::is_scatter_view_v<dual_view_t>);
static_assert(!Kokkos::is_view_v<dual_view_t>);

static_assert(!Kokkos::is_dual_view_v<dyn_rank_view_t>);
static_assert(Kokkos::is_dyn_rank_view_v<dyn_rank_view_t>);
static_assert(!Kokkos::is_dynamic_view_v<dyn_rank_view_t>);
static_assert(!Kokkos::Experimental::is_offset_view_v<dyn_rank_view_t>);
static_assert(!Kokkos::Experimental::is_scatter_view_v<dyn_rank_view_t>);
static_assert(!Kokkos::is_view_v<dyn_rank_view_t>);

static_assert(!Kokkos::is_dual_view_v<dynamic_view_t>);
static_assert(!Kokkos::is_dyn_rank_view_v<dynamic_view_t>);
static_assert(Kokkos::is_dynamic_view_v<dynamic_view_t>);
static_assert(!Kokkos::Experimental::is_offset_view_v<dynamic_view_t>);
static_assert(!Kokkos::Experimental::is_scatter_view_v<dynamic_view_t>);
static_assert(!Kokkos::is_view_v<dynamic_view_t>);

static_assert(!Kokkos::is_dual_view_v<offset_view_t>);
static_assert(!Kokkos::is_dyn_rank_view_v<offset_view_t>);
static_assert(!Kokkos::is_dynamic_view_v<offset_view_t>);
static_assert(Kokkos::Experimental::is_offset_view_v<offset_view_t>);
static_assert(!Kokkos::Experimental::is_scatter_view_v<offset_view_t>);
static_assert(!Kokkos::is_view_v<offset_view_t>);

static_assert(!Kokkos::is_dual_view_v<scatter_view_t>);
static_assert(!Kokkos::is_dyn_rank_view_v<scatter_view_t>);
static_assert(!Kokkos::is_dynamic_view_v<scatter_view_t>);
static_assert(!Kokkos::Experimental::is_offset_view_v<scatter_view_t>);
static_assert(Kokkos::Experimental::is_scatter_view_v<scatter_view_t>);
static_assert(!Kokkos::is_view_v<scatter_view_t>);

}  // namespace
