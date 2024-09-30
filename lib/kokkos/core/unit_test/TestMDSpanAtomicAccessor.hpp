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

#include <gtest/gtest.h>
#ifndef KOKKOS_ENABLE_CXX17
#include <concepts>
#endif

template <class T, class ExecutionSpace>
void test_atomic_accessor() {
  using value_type = std::remove_const_t<T>;
  Kokkos::View<value_type*, ExecutionSpace> v("V", 100);

  Kokkos::parallel_for(
      Kokkos::RangePolicy<ExecutionSpace>(0, v.extent(0)),
      KOKKOS_LAMBDA(int i) { v(i) = i; });

  int errors;
  using acc_t = Kokkos::Impl::AtomicAccessorRelaxed<T>;
  acc_t acc{};
  typename acc_t::data_handle_type ptr = v.data();

  Kokkos::parallel_reduce(
      Kokkos::RangePolicy<ExecutionSpace>(0, v.extent(0)),
      KOKKOS_LAMBDA(int i, int& error) {
        if (acc.access(ptr, i) != ptr[i]) error++;
        if (acc.offset(ptr, i) != ptr + i) error++;
        static_assert(std::is_same_v<typename acc_t::element_type, T>);
        static_assert(
            std::is_same_v<typename acc_t::reference,
                           desul::AtomicRef<T, desul::MemoryOrderRelaxed,
                                            desul::MemoryScopeDevice>>);
        static_assert(std::is_same_v<typename acc_t::data_handle_type, T*>);
        static_assert(std::is_same_v<typename acc_t::offset_policy, acc_t>);
        static_assert(std::is_same_v<decltype(acc.access(ptr, i)),
                                     typename acc_t::reference>);
        static_assert(std::is_same_v<decltype(acc.offset(ptr, i)), T*>);
        static_assert(std::is_nothrow_move_constructible_v<acc_t>);
        static_assert(std::is_nothrow_move_assignable_v<acc_t>);
        static_assert(std::is_nothrow_swappable_v<acc_t>);
        static_assert(std::is_trivially_copyable_v<acc_t>);
        static_assert(std::is_trivially_default_constructible_v<acc_t>);
        static_assert(std::is_trivially_constructible_v<acc_t>);
        static_assert(std::is_trivially_move_constructible_v<acc_t>);
        static_assert(std::is_trivially_assignable_v<acc_t, acc_t>);
        static_assert(std::is_trivially_move_assignable_v<acc_t>);
#ifndef KOKKOS_ENABLE_CXX17
        static_assert(std::copyable<acc_t>);
        static_assert(std::is_empty_v<acc_t>);
#endif
      },
      errors);
  ASSERT_EQ(errors, 0);
}

void test_atomic_accessor_conversion() {
  using ExecutionSpace = TEST_EXECSPACE;
  using T              = float;
  using acc_t          = Kokkos::Impl::AtomicAccessorRelaxed<T>;
  using const_acc_t    = Kokkos::Impl::AtomicAccessorRelaxed<const T>;
  using int_acc_t      = Kokkos::Impl::AtomicAccessorRelaxed<int>;
  using defacc_t       = Kokkos::default_accessor<T>;
  using const_defacc_t = Kokkos::default_accessor<const T>;
  using int_defacc_t   = Kokkos::default_accessor<int>;

  Kokkos::parallel_for(
      Kokkos::RangePolicy<ExecutionSpace>(0, 1), KOKKOS_LAMBDA(int) {
        static_assert(std::is_constructible_v<const_acc_t, acc_t>);
        static_assert(std::is_convertible_v<acc_t, const_acc_t>);
        static_assert(!std::is_constructible_v<acc_t, const_acc_t>);
        static_assert(!std::is_constructible_v<acc_t, int_acc_t>);
        static_assert(std::is_constructible_v<defacc_t, acc_t>);
        static_assert(std::is_constructible_v<acc_t, defacc_t>);
        static_assert(!std::is_constructible_v<int_defacc_t, acc_t>);
        static_assert(!std::is_constructible_v<int_acc_t, defacc_t>);
        static_assert(std::is_constructible_v<const_defacc_t, const_acc_t>);
        static_assert(std::is_constructible_v<const_acc_t, const_defacc_t>);
        static_assert(std::is_constructible_v<const_defacc_t, acc_t>);
        static_assert(std::is_constructible_v<const_acc_t, defacc_t>);
        static_assert(!std::is_constructible_v<defacc_t, const_acc_t>);
        static_assert(!std::is_constructible_v<acc_t, const_defacc_t>);
        static_assert(std::is_convertible_v<acc_t, const_acc_t>);
        static_assert(std::is_convertible_v<defacc_t, acc_t>);
        static_assert(std::is_convertible_v<defacc_t, const_acc_t>);
        static_assert(std::is_convertible_v<const_defacc_t, const_acc_t>);
        static_assert(!std::is_convertible_v<acc_t, defacc_t>);
        static_assert(!std::is_convertible_v<acc_t, const_defacc_t>);
        static_assert(!std::is_convertible_v<const_acc_t, const_defacc_t>);
      });
}

TEST(TEST_CATEGORY, mdspan_atomic_accessor) {
  using ExecutionSpace = TEST_EXECSPACE;
  test_atomic_accessor<int, ExecutionSpace>();
  test_atomic_accessor<double, ExecutionSpace>();
}
