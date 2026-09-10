// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#ifdef KOKKOS_ENABLE_IMPL_VIEW_LEGACY
import kokkos.core_impl;
#endif
#else
#include <Kokkos_Core.hpp>
#endif

#include <desul/atomics.hpp>

#include <type_traits>

namespace {

using test_atomic_view =
    Kokkos::View<double *, Kokkos::MemoryTraits<Kokkos::Atomic>>;
static_assert(
    std::is_same_v<
        decltype(std::declval<test_atomic_view>()(std::declval<int>())),
#ifdef KOKKOS_ENABLE_IMPL_VIEW_LEGACY
        Kokkos::Impl::AtomicDataElement<
            Kokkos::ViewTraits<double *, Kokkos::MemoryTraits<Kokkos::Atomic>>>
#else
        desul::AtomicRef<double, desul::MemoryOrderRelaxed,
                         desul::MemoryScopeDevice>
#endif
        >);

}  // namespace
