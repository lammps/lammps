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

#include <gtest/gtest.h>
#include <concepts>

template <class T>
struct funky_data_handle {
  T* val;

  KOKKOS_FUNCTION
  operator T*() { return val; }
  KOKKOS_FUNCTION
  operator const T*() const { return val; }
};

template <class ElementType>
struct FunkyAcc {
  using element_type     = ElementType;
  using reference        = std::conditional_t<std::is_const_v<element_type>,
                                       element_type, element_type&>;
  using data_handle_type = funky_data_handle<element_type>;
  using offset_policy    = Kokkos::default_accessor<element_type>;
  KOKKOS_FUNCTION
  reference access(data_handle_type p, size_t i) const { return p.val[i]; }
  KOKKOS_FUNCTION
  element_type* offset(data_handle_type p, size_t i) const { return p.val + i; }
};

template <class T, class ExecutionSpace,
          class MemorySpace = typename ExecutionSpace::memory_space>
void test_space_aware_accessor() {
  using memory_space_t = MemorySpace;
  using value_type     = std::remove_const_t<T>;
  Kokkos::View<value_type*, ExecutionSpace> v("V", 100);

  Kokkos::parallel_for(
      Kokkos::RangePolicy<ExecutionSpace>(0, v.extent(0)),
      KOKKOS_LAMBDA(int i) { v(i) = i; });

  int errors;
  using acc_t = Kokkos::Impl::SpaceAwareAccessor<memory_space_t, FunkyAcc<T>>;
  acc_t acc{};
  typename acc_t::data_handle_type ptr{v.data()};

  Kokkos::parallel_reduce(
      Kokkos::RangePolicy<ExecutionSpace>(0, v.extent(0)),
      KOKKOS_LAMBDA(int i, int& error) {
        if (acc.access(ptr, i) != ptr[i]) error++;
        if (acc.offset(ptr, i) != ptr + i) error++;
        static_assert(std::is_same_v<typename acc_t::element_type, T>);
        if constexpr (std::is_const_v<T>) {
          static_assert(std::is_same_v<typename acc_t::reference, T>);
        } else {
          static_assert(std::is_same_v<typename acc_t::reference, T&>);
        }
        static_assert(std::is_same_v<typename acc_t::data_handle_type,
                                     funky_data_handle<T>>);
        static_assert(
            std::is_same_v<typename acc_t::offset_policy,
                           Kokkos::Impl::SpaceAwareAccessor<
                               memory_space_t, Kokkos::default_accessor<T>>>);
        if constexpr (std::is_const_v<T>) {
          static_assert(std::is_same_v<decltype(acc.access(ptr, i)),
                                       std::remove_const_t<T>>);
        } else {
          static_assert(std::is_same_v<decltype(acc.access(ptr, i)), T&>);
        }
        static_assert(std::is_same_v<decltype(acc.offset(ptr, i)), T*>);
        static_assert(std::is_same_v<decltype(acc.nested_accessor()),
                                     const FunkyAcc<T>&>);
        static_assert(std::is_nothrow_move_constructible_v<acc_t>);
        static_assert(std::is_nothrow_move_assignable_v<acc_t>);
        static_assert(std::is_nothrow_swappable_v<acc_t>);
        static_assert(
            std::is_same_v<typename acc_t::memory_space, memory_space_t>);
        static_assert(
            std::is_same_v<typename acc_t::nested_accessor_type, FunkyAcc<T>>);
        static_assert(std::copyable<acc_t>);
// Windows and nvcc >= 12.9 (separately) don't treat no-unique-address correctly
#if !defined(_WIN32) &&                                               \
    !(defined(KOKKOS_ENABLE_CUDA) && defined(KOKKOS_COMPILER_NVCC) && \
      (KOKKOS_COMPILER_NVCC >= 1290))
        static_assert(std::is_empty_v<acc_t>);
#endif
      },
      errors);
  ASSERT_EQ(errors, 0);
}

void test_space_aware_accessor_conversion() {
  using ExecutionSpace = TEST_EXECSPACE;
  using memory_space_t = typename ExecutionSpace::memory_space;
  using T              = float;
  using acc_t          = Kokkos::Impl::SpaceAwareAccessor<memory_space_t,
                                                 Kokkos::default_accessor<T>>;
  using const_acc_t =
      Kokkos::Impl::SpaceAwareAccessor<memory_space_t,
                                       Kokkos::default_accessor<const T>>;
  using int_acc_t =
      Kokkos::Impl::SpaceAwareAccessor<memory_space_t,
                                       Kokkos::default_accessor<int>>;
  using host_acc_t =
      Kokkos::Impl::SpaceAwareAccessor<Kokkos::HostSpace,
                                       Kokkos::default_accessor<T>>;
  using anon_acc_t =
      Kokkos::Impl::SpaceAwareAccessor<Kokkos::AnonymousSpace,
                                       Kokkos::default_accessor<T>>;

  Kokkos::parallel_for(
      Kokkos::RangePolicy<ExecutionSpace>(0, 1), KOKKOS_LAMBDA(int) {
        static_assert(std::is_constructible_v<const_acc_t, acc_t>);
        static_assert(std::is_convertible_v<acc_t, const_acc_t>);
        static_assert(!std::is_constructible_v<acc_t, const_acc_t>);
        static_assert(!std::is_constructible_v<acc_t, int_acc_t>);
        static_assert(
            std::is_constructible_v<acc_t, host_acc_t> ==
            Kokkos::Impl::MemorySpaceAccess<memory_space_t,
                                            Kokkos::HostSpace>::assignable);
        static_assert(
            std::is_constructible_v<host_acc_t, acc_t> ==
            Kokkos::Impl::MemorySpaceAccess<Kokkos::HostSpace,
                                            memory_space_t>::assignable);
        static_assert(std::is_constructible_v<anon_acc_t, acc_t>);
        static_assert(std::is_constructible_v<acc_t, anon_acc_t>);
        static_assert(std::is_convertible_v<anon_acc_t, acc_t>);
        static_assert(std::is_convertible_v<acc_t, anon_acc_t>);
      });
}

TEST(TEST_CATEGORY, mdspan_space_aware_accessor) {
  using ExecutionSpace = TEST_EXECSPACE;
  test_space_aware_accessor<int, ExecutionSpace>();
  test_space_aware_accessor<double, ExecutionSpace>();
  test_space_aware_accessor<const int, ExecutionSpace>();
  test_space_aware_accessor<const double, ExecutionSpace>();
  test_space_aware_accessor<double, ExecutionSpace, Kokkos::AnonymousSpace>();
  test_space_aware_accessor<const int, ExecutionSpace,
                            Kokkos::AnonymousSpace>();
  test_space_aware_accessor_conversion();
}
