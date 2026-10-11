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

#define CHECK_SAME(type) \
  static_assert(         \
      std::is_same_v<typename new_view_t::type, typename old_view_t::type>)

template <class T, class Extents, class DataType, class ArrayLayout,
          class Space, class MemTraits>
constexpr bool test_equivalence() {
  using layout_type =
      typename Kokkos::Impl::LayoutFromArrayLayout<ArrayLayout>::type;
  using new_view_t =
      Kokkos::View<T, Extents, layout_type,
                   Kokkos::Experimental::Accessor<
                       T, typename Space::memory_space, MemTraits>>;
  using old_view_t = Kokkos::View<DataType, ArrayLayout, Space, MemTraits>;

  // Primary test: this is equivalent to checking that they derive from the same
  // BasicView implementation, since mdspan_type comes from BasicView and uses
  // directly the template args of BasicView
  CHECK_SAME(mdspan_type);

  // Checking typedefs which come from ViewTraits and not from BasicView
  // Not currently the same
  // CHECK_SAME(traits);
  CHECK_SAME(value_type);
  CHECK_SAME(const_value_type);
  CHECK_SAME(non_const_value_type);
  CHECK_SAME(data_type);
  CHECK_SAME(const_data_type);
  CHECK_SAME(non_const_data_type);
  // Not currently the same
  // CHECK_SAME(view_tracker_type);
  CHECK_SAME(array_layout);
  CHECK_SAME(device_type);
  CHECK_SAME(execution_space);
  CHECK_SAME(memory_space);
  // Not currently the same
  // CHECK_SAME(memory_traits);
  CHECK_SAME(host_mirror_space);

  return true;
}

using LL           = Kokkos::LayoutLeft;
using LR           = Kokkos::LayoutRight;
using LS           = Kokkos::LayoutStride;
using ED           = Kokkos::DefaultExecutionSpace;
using EH           = Kokkos::DefaultHostExecutionSpace;
using SH           = Kokkos::HostSpace;
using SD           = typename Kokkos::DefaultExecutionSpace::memory_space;
using DD           = typename Kokkos::DefaultExecutionSpace::device_type;
using f            = float;
using cf           = const float;
constexpr size_t d = Kokkos::dynamic_extent;

// clang-format off

static_assert(test_equivalence<f,  Kokkos::dextents<size_t, 0>,         f,          LL, ED, Kokkos::MemoryTraits<>>());
static_assert(test_equivalence<cf, Kokkos::dextents<size_t, 2>,         cf**,       LS, EH, Kokkos::MemoryTraits<Kokkos::Unmanaged>>());
static_assert(test_equivalence<f,  Kokkos::extents<size_t, 2, 3>,       f[2][3],    LR, SH, Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::Atomic>>());
static_assert(test_equivalence<cf, Kokkos::dextents<size_t, 8>,         cf********, LS, SD, Kokkos::MemoryTraits<Kokkos::Atomic | Kokkos::RandomAccess>>());
static_assert(test_equivalence<f,  Kokkos::extents<size_t, d, d, 2, 3>, f**[2][3],  LR, DD, Kokkos::MemoryTraits<Kokkos::RandomAccess>>());
// This can't give the same
//static_assert(test_equivalence<f,  Kokkos::extents<size_t, 3, d, 2, 3>, f**[2][3],     LR, HS, Kokkos::MemoryTraits<Kokkos::RandomAccess>>());

// clang-format on

}  // namespace
