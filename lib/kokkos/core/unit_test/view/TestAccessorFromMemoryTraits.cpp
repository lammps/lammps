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

// Checking requirement of explict type conversion to View

namespace {

template <class T, class Space>
constexpr bool test_equivalence() {
  // Checking that the aliases lead to expected accessor
  static_assert(std::is_same_v<
                Kokkos::Experimental::Accessor<T>,
                Kokkos::Impl::CheckedReferenceCountedAccessor<
                    T, typename Kokkos::DefaultExecutionSpace::memory_space>>);
  static_assert(std::is_same_v<Kokkos::Experimental::Accessor<T, Space>,
                               Kokkos::Impl::CheckedReferenceCountedAccessor<
                                   T, typename Space::memory_space>>);
  static_assert(
      std::is_same_v<
          Kokkos::Experimental::Accessor<T, Space, Kokkos::MemoryTraits<>>,
          Kokkos::Impl::CheckedReferenceCountedAccessor<
              T, typename Space::memory_space>>);
  static_assert(std::is_same_v<
                Kokkos::Experimental::Accessor<
                    T, Space,
                    Kokkos::MemoryTraits<Kokkos::Unmanaged | Kokkos::Atomic>>,
                Kokkos::Impl::CheckedRelaxedAtomicAccessor<
                    T, typename Space::memory_space>>);
  static_assert(
      std::is_same_v<Kokkos::Experimental::Accessor<
                         T, Space, Kokkos::MemoryTraits<Kokkos::Atomic>>,
                     Kokkos::Impl::CheckedReferenceCountedRelaxedAtomicAccessor<
                         T, typename Space::memory_space>>);
  static_assert(std::is_same_v<
                Kokkos::Experimental::Accessor<
                    T, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>,
                Kokkos::Impl::SpaceAwareAccessor<typename Space::memory_space,
                                                 Kokkos::default_accessor<T>>>);
  static_assert(
      std::is_same_v<
          Kokkos::Experimental::Accessor<
              T, Space,
              Kokkos::MemoryTraits<Kokkos::Atomic | Kokkos::RandomAccess>>,
          Kokkos::Impl::CheckedReferenceCountedRelaxedAtomicAccessor<
              T, typename Space::memory_space>>);

  // Checking MemoryTraits behavior of our accessors
  // Default ones should stay default
  static_assert(
      std::is_same_v<
          Kokkos::Impl::MemoryTraitsFromAccessor<
              Kokkos::Experimental::Accessor<T, Space, Kokkos::MemoryTraits<>>>,
          Kokkos::MemoryTraits<>>);
  // Unmanaged and Atomic are propagated now
  static_assert(
      std::is_same_v<
          Kokkos::Impl::MemoryTraitsFromAccessor<Kokkos::Experimental::Accessor<
              T, Space, Kokkos::MemoryTraits<Kokkos::Unmanaged>>>,
          Kokkos::MemoryTraits<Kokkos::Unmanaged>>);
  static_assert(
      std::is_same_v<
          Kokkos::Impl::MemoryTraitsFromAccessor<Kokkos::Experimental::Accessor<
              T, Space, Kokkos::MemoryTraits<Kokkos::Atomic>>>,
          Kokkos::MemoryTraits<Kokkos::Atomic>>);
  // RandomAccess is dropped, since no accessor currently implements this
  static_assert(
      std::is_same_v<
          Kokkos::Impl::MemoryTraitsFromAccessor<Kokkos::Experimental::Accessor<
              T, Space,
              Kokkos::MemoryTraits<Kokkos::Atomic | Kokkos::RandomAccess>>>,
          Kokkos::MemoryTraits<Kokkos::Atomic>>);
  return true;
}

static_assert(test_equivalence<float, Kokkos::DefaultExecutionSpace>());
static_assert(
    test_equivalence<const float, Kokkos::DefaultHostExecutionSpace>());
static_assert(test_equivalence<Kokkos::complex<float>, Kokkos::HostSpace>());
}  // namespace
