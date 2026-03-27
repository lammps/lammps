// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
import kokkos.core_impl;
#else
#include <Kokkos_Core.hpp>
#endif

#include <gtest/gtest.h>

namespace Foo {
// This accessor doesn't do anything useful.
// I just need one to test that it gets properly constructed from accessor_arg.
template <class ElementType, class MemorySpace>
struct TestAccessor {
  using element_type = ElementType;
  using reference    = element_type&;
  using data_handle_type =
      Kokkos::Impl::ReferenceCountedDataHandle<ElementType, MemorySpace>;
  using offset_policy = TestAccessor;

  // View expects this from accessors right now
  using memory_space = MemorySpace;

  KOKKOS_DEFAULTED_FUNCTION
  constexpr TestAccessor() = default;

  template <class OtherElementType,
            std::enable_if_t<std::is_constructible_v<
                                 Kokkos::default_accessor<ElementType>,
                                 Kokkos::default_accessor<OtherElementType>>,
                             int> = 0>
  KOKKOS_FUNCTION constexpr TestAccessor(
      const TestAccessor<OtherElementType, MemorySpace>& other) noexcept
      : value(other.value) {}

  KOKKOS_FUNCTION
  TestAccessor(const size_t val) : value(val) {}

  KOKKOS_FUNCTION
  constexpr reference access(
#ifndef KOKKOS_ENABLE_OPENACC
      const data_handle_type& p,
#else
      // FIXME OPENACC: illegal address when passing by reference
      data_handle_type p,
#endif
      size_t i) const noexcept {
    return *(p + (i % value));
  }

  KOKKOS_FUNCTION
  constexpr typename offset_policy::data_handle_type offset(
#ifndef KOKKOS_ENABLE_OPENACC
      const data_handle_type& p,
#else
      // FIXME OPENACC: illegal address when passing by reference
      data_handle_type p,
#endif
      size_t i) const noexcept {
    return p + i;
  }

  size_t value{};
};

struct Bar {
  int val;
};

// Use the customization point to inject the custom accessor
template <class LayoutType, class DeviceType, class MemoryTraits>
constexpr auto customize_view_arguments(
    Kokkos::Impl::ViewArguments<Bar, LayoutType, DeviceType, MemoryTraits>) {
  return Kokkos::Impl::ViewCustomArguments<
      size_t, TestAccessor<Bar, typename DeviceType::memory_space>>{};
}
}  // namespace Foo

void test_accessor_arg() {
  using view_t = Kokkos::View<Foo::Bar*, TEST_EXECSPACE>;

  // Make sure I got the accessor I execpt
  static_assert(
      std::is_same_v<
          view_t::accessor_type,
          Foo::TestAccessor<Foo::Bar, typename TEST_EXECSPACE::memory_space>>);

  // accessor will be default constructed
  view_t a("A", 10);
  // Cheap way of checking I get the expected accessor.
  ASSERT_EQ(a.accessor().value, 0ul);

  view_t b(a.data(), 10);
  ASSERT_EQ(b.accessor().value, 0ul);

  // accessor will be constructed from AccessorArg_t
  view_t c(Kokkos::view_alloc("C", Kokkos::Impl::AccessorArg_t{5ul}), 10);
  ASSERT_EQ(c.accessor().value, 5ul);

  // Test copy ctor to make sure the customize_view_arguments thing doesn't
  // interfere
  view_t c_copy = c;
  ASSERT_EQ(c_copy.accessor().value, 5ul);

  view_t d(Kokkos::view_wrap(c.data(), Kokkos::Impl::AccessorArg_t{7ul}), 10);
  ASSERT_EQ(d.accessor().value, 7ul);

  // Test unmanaged ctors on GPU too (if GPU is enabled)
  int num_error = 0;
  Kokkos::parallel_reduce(
      "test_accessor_arg", 1,
      KOKKOS_LAMBDA(int, int& errors) {
        view_t e(a.data(), 10);
        if (e.accessor().value != 0lu) errors++;
        view_t f(Kokkos::view_wrap(e.data(), Kokkos::Impl::AccessorArg_t{7ul}),
                 10);
        if (f.accessor().value != 7lu) errors++;
        // Test copy ctor to make sure the customize_view_arguments thing
        // doesn't interfere
        view_t f_copy = f;
        if (f_copy.accessor().value != 7lu) errors++;
      },
      num_error);
  ASSERT_EQ(num_error, 0);
}

TEST(TEST_CATEGORY, view_customization_accessor_arg) { test_accessor_arg(); }
