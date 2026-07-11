// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
import kokkos.dual_view;
#else
#include <Kokkos_Core.hpp>
#include <Kokkos_DualView.hpp>
#endif

#include <type_traits>

namespace {

template <class DataType, class Arg1Type = void, class Arg2Type = void,
          class Arg3Type = void>
void not_supported_anymore(
    Kokkos::DualView<DataType, Arg1Type, Arg2Type, Arg2Type> x) {
  static_assert(Kokkos::is_dual_view_v<decltype(x)>);
}

template <class DataType, class... Properties>
void prefer_instead(Kokkos::DualView<DataType, Properties...> x) {
  static_assert(Kokkos::is_dual_view_v<decltype(x)>);
}

using KDV = Kokkos::DualView<int*>;

#ifdef KOKKOS_ENABLE_DEPRECATED_CODE_4
static_assert(
    std::is_void_v<decltype(not_supported_anymore(std::declval<KDV>()))>);
#endif

static_assert(std::is_void_v<decltype(prefer_instead(std::declval<KDV>()))>);

}  // namespace
