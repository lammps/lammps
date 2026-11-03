// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
import kokkos.core_impl;
#else
#include <Kokkos_Core.hpp>
#endif

#include <cstdlib>
#include <type_traits>

namespace {

template <class T, class Extents, class Layout, class MemSpace, class MemTraits>
constexpr bool test_create_mirror_types() {
  using view_t =
      Kokkos::View<T, Extents, Layout,
                   Kokkos::Experimental::Accessor<T, MemSpace, MemTraits>>;
  constexpr bool host_accessible =
      Kokkos::SpaceAccessibility<Kokkos::DefaultHostExecutionSpace,
                                 MemSpace>::accessible;

  constexpr bool const_T = std::is_const_v<T>;

  using host_mirror_space_t =
      std::conditional_t<host_accessible, MemSpace, Kokkos::HostSpace>;

  using create_mirror_t = std::conditional_t<
      host_accessible && !const_T, view_t,
      Kokkos::View<std::remove_const_t<T>, Extents, Layout,
                   Kokkos::Experimental::Accessor<std::remove_const_t<T>,
                                                  host_mirror_space_t,
                                                  Kokkos::MemoryTraits<>>>>;

  using create_mirror_view_t = std::conditional_t<
      host_accessible, view_t,
      Kokkos::View<std::remove_const_t<T>, Extents, Layout,
                   Kokkos::Experimental::Accessor<std::remove_const_t<T>,
                                                  host_mirror_space_t,
                                                  Kokkos::MemoryTraits<>>>>;

  static_assert(std::is_same_v<create_mirror_t, decltype(Kokkos::create_mirror(
                                                    std::declval<view_t>()))>);
  static_assert(
      std::is_same_v<create_mirror_t, decltype(Kokkos::create_mirror_view(
                                          std::declval<view_t>()))>);
  static_assert(std::is_same_v<create_mirror_view_t,
                               decltype(Kokkos::create_mirror_view_and_copy(
                                   std::declval<view_t>()))>);

  using host_create_mirror_t = Kokkos::View<
      std::remove_const_t<T>, Extents, Layout,
      Kokkos::Experimental::Accessor<std::remove_const_t<T>, Kokkos::HostSpace,
                                     Kokkos::MemoryTraits<>>>;
  static_assert(std::is_same_v<host_create_mirror_t,
                               decltype(Kokkos::create_mirror(
                                   Kokkos::HostSpace(), view_t()))>);
  using host_create_mirror_view_t = std::conditional_t<
      std::is_same_v<Kokkos::HostSpace, MemSpace>, view_t,
      Kokkos::View<std::remove_const_t<T>, Extents, Layout,
                   Kokkos::Experimental::Accessor<std::remove_const_t<T>,
                                                  Kokkos::HostSpace,
                                                  Kokkos::MemoryTraits<>>>>;
  static_assert(std::is_same_v<host_create_mirror_view_t,
                               decltype(Kokkos::create_mirror_view(
                                   Kokkos::HostSpace(), view_t()))>);
  static_assert(std::is_same_v<host_create_mirror_view_t,
                               decltype(Kokkos::create_mirror_view_and_copy(
                                   Kokkos::HostSpace(), view_t()))>);

#ifdef KOKKOS_HAS_SHARED_SPACE
  using shared_create_mirror_t = Kokkos::View<
      std::remove_const_t<T>, Extents, Layout,
      Kokkos::Experimental::Accessor<
          std::remove_const_t<T>, Kokkos::SharedSpace, Kokkos::MemoryTraits<>>>;
  using shared_create_mirror_view_t = std::conditional_t<
      std::is_same_v<MemSpace, Kokkos::SharedSpace>, view_t,
      Kokkos::View<std::remove_const_t<T>, Extents, Layout,
                   Kokkos::Experimental::Accessor<std::remove_const_t<T>,
                                                  Kokkos::SharedSpace,
                                                  Kokkos::MemoryTraits<>>>>;
  static_assert(std::is_same_v<shared_create_mirror_t,
                               decltype(Kokkos::create_mirror(
                                   Kokkos::SharedSpace(), view_t()))>);
  static_assert(std::is_same_v<shared_create_mirror_view_t,
                               decltype(Kokkos::create_mirror_view(
                                   Kokkos::SharedSpace(), view_t()))>);
  static_assert(std::is_same_v<shared_create_mirror_view_t,
                               decltype(Kokkos::create_mirror_view_and_copy(
                                   Kokkos::SharedSpace(), view_t()))>);
#endif
  return true;
}

template <class T, class Extents, class Layout, class MemSpace, class MemTraits>
constexpr bool test_typedefs() {
  using view_t =
      Kokkos::View<T, Extents, Layout,
                   Kokkos::Experimental::Accessor<T, MemSpace, MemTraits>>;
  static_assert(std::is_same_v<view_t, typename view_t::type>);

  using const_view_t = Kokkos::View<
      const T, Extents, Layout,
      Kokkos::Experimental::Accessor<const T, MemSpace, MemTraits>>;
  static_assert(std::is_same_v<const_view_t, typename view_t::const_type>);

  using non_const_view_t =
      Kokkos::View<std::remove_const_t<T>, Extents, Layout,
                   Kokkos::Experimental::Accessor<std::remove_const_t<T>,
                                                  MemSpace, MemTraits>>;
  static_assert(
      std::is_same_v<non_const_view_t, typename view_t::non_const_type>);

  return true;
}

template <class T, class Extents, class Layout, class MemSpace, class MemTraits>
constexpr bool test_derived_types() {
  constexpr bool typedef_result =
      test_typedefs<T, Extents, Layout, MemSpace, MemTraits>();

  if constexpr (MemTraits::is_atomic) {
    return typedef_result;
  } else {
    return typedef_result &&
           test_create_mirror_types<T, Extents, Layout, MemSpace, MemTraits>();
  }
}

using LL           = Kokkos::layout_left;
using LR           = Kokkos::layout_right;
using SH           = Kokkos::HostSpace;
using SD           = typename Kokkos::DefaultExecutionSpace::memory_space;
using f            = float;
using cf           = const float;
constexpr size_t d = Kokkos::dynamic_extent;

// clang-format off

static_assert(test_derived_types<f,  Kokkos::dextents<unsigned, 0>,         LL, SD, Kokkos::MemoryTraits<>>());
static_assert(test_derived_types<cf, Kokkos::dextents<unsigned, 2>,         LR, SD, Kokkos::MemoryTraits<Kokkos::Unmanaged>>());
static_assert(test_derived_types<f,  Kokkos::extents<unsigned, 2, 3>,       LL, SD, Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::Atomic>>());
static_assert(test_derived_types<cf, Kokkos::dextents<unsigned, 8>,         LR, SD, Kokkos::MemoryTraits<Kokkos::Atomic | Kokkos::RandomAccess>>());
static_assert(test_derived_types<f,  Kokkos::extents<unsigned, d, d, 2, 3>, LL, SD, Kokkos::MemoryTraits<Kokkos::RandomAccess>>());

// FIXME_OPENACC, FIXME_NVHPC This particular case causes internal compiler errors with NVHPC 23.7
#ifndef KOKKOS_ENABLE_OPENACC
static_assert(test_derived_types<f,  Kokkos::dextents<unsigned, 0>,         LL, SH, Kokkos::MemoryTraits<>>());
#endif
static_assert(test_derived_types<cf, Kokkos::dextents<unsigned, 2>,         LR, SH, Kokkos::MemoryTraits<Kokkos::Unmanaged>>());
static_assert(test_derived_types<f,  Kokkos::extents<unsigned, 2, 3>,       LL, SH, Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::Atomic>>());
static_assert(test_derived_types<cf, Kokkos::dextents<unsigned, 8>,         LR, SH, Kokkos::MemoryTraits<Kokkos::Atomic | Kokkos::RandomAccess>>());
static_assert(test_derived_types<f,  Kokkos::extents<unsigned, d, d, 2, 3>, LL, SH, Kokkos::MemoryTraits<Kokkos::RandomAccess>>());

#ifdef KOKKOS_HAS_SHARED_SPACE
using SS = Kokkos::SharedSpace;

static_assert(test_derived_types<f,  Kokkos::dextents<unsigned, 0>,         LL, SS, Kokkos::MemoryTraits<>>());
static_assert(test_derived_types<cf, Kokkos::dextents<unsigned, 2>,         LR, SS, Kokkos::MemoryTraits<Kokkos::Unmanaged>>());
static_assert(test_derived_types<f,  Kokkos::extents<unsigned, 2, 3>,       LL, SS, Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::Atomic>>());
static_assert(test_derived_types<cf, Kokkos::dextents<unsigned, 8>,         LR, SS, Kokkos::MemoryTraits<Kokkos::Atomic | Kokkos::RandomAccess>>());
static_assert(test_derived_types<f,  Kokkos::extents<unsigned, d, d, 2, 3>, LL, SS, Kokkos::MemoryTraits<Kokkos::RandomAccess>>());
#endif

// clang-format on

}  // namespace
