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
#include <type_traits>

namespace Foo {

struct FooVal {
  int value;
};

struct BarVal {
  int value;
};
// Customization point to control mdspan arguments from view arguments
// Default implementation returns void to indicate no customization
template <class LayoutType, class DeviceType, class MemoryTraits>
constexpr auto customize_view_arguments(
    Kokkos::Impl::ViewArguments<BarVal, LayoutType, DeviceType, MemoryTraits>) {
  using mem_space_t = typename DeviceType::memory_space;
  return Kokkos::Impl::ViewCustomArguments<
      int, Kokkos::Impl::SpaceAwareAccessor<
               mem_space_t, Kokkos::default_accessor<BarVal>>>{};
}
template <class LayoutType, class DeviceType, class MemoryTraits>
constexpr auto customize_view_arguments(
    Kokkos::Impl::ViewArguments<const BarVal, LayoutType, DeviceType,
                                MemoryTraits>) {
  using mem_space_t = typename DeviceType::memory_space;
  return Kokkos::Impl::ViewCustomArguments<
      unsigned, Kokkos::Impl::SpaceAwareAccessor<
                    mem_space_t, Kokkos::default_accessor<const BarVal>>>{};
}

}  // namespace Foo

using view_fooval_t  = Kokkos::View<Foo::FooVal*>;
using view_barval_t  = Kokkos::View<Foo::BarVal*>;
using view_cbarval_t = Kokkos::View<const Foo::BarVal*>;

static_assert(!view_fooval_t::traits::impl_is_customized);
static_assert(view_barval_t::traits::impl_is_customized);
static_assert(view_cbarval_t::traits::impl_is_customized);

static_assert(
    std::is_same_v<typename view_fooval_t::extents_type::index_type, size_t>);
static_assert(
    std::is_same_v<typename view_barval_t::extents_type::index_type, int>);
static_assert(std::is_same_v<typename view_cbarval_t::extents_type::index_type,
                             unsigned>);

using mem_space_t = typename Kokkos::DefaultExecutionSpace::memory_space;
static_assert(std::is_same_v<typename view_fooval_t::accessor_type,
                             Kokkos::Impl::CheckedReferenceCountedAccessor<
                                 Foo::FooVal, mem_space_t>>);
static_assert(
    std::is_same_v<typename view_barval_t::accessor_type,
                   Kokkos::Impl::SpaceAwareAccessor<
                       mem_space_t, Kokkos::default_accessor<Foo::BarVal>>>);
static_assert(std::is_same_v<
              typename view_cbarval_t::accessor_type,
              Kokkos::Impl::SpaceAwareAccessor<
                  mem_space_t, Kokkos::default_accessor<const Foo::BarVal>>>);
