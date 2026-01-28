// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Core.hpp>
#include <type_traits>

// Checking requirement of explict type conversion to View

namespace {

template <class From, class To>
inline constexpr bool conversion_is_explicit =
    std::is_constructible_v<To, From> && !std::is_convertible_v<From, To>;

using T = int;
static_assert(conversion_is_explicit<T*, Kokkos::View<T>>);
static_assert(conversion_is_explicit<T*, Kokkos::View<T*>>);
static_assert(conversion_is_explicit<T*, Kokkos::View<T[4]>>);
static_assert(conversion_is_explicit<std::nullptr_t, Kokkos::View<T[4]>>);

}  // namespace
