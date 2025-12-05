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

#include <gtest/gtest.h>

namespace Foo {

// ElementType which is actually just a tag type
struct Bar {
  float vals[5];
  KOKKOS_FUNCTION
  float& operator[](size_t idx) { return vals[idx % 5]; }
  KOKKOS_FUNCTION
  const float& operator[](size_t idx) const { return vals[idx % 5]; }
};

// Reference type for Bar
template <class T>
struct BarRef {
  T* ptr{nullptr};
  size_t size{0ul};

  KOKKOS_FUNCTION
  T& operator[](size_t idx) const { return ptr[idx]; }

  KOKKOS_FUNCTION
  operator Bar() const {
    Bar val;
    for (size_t i = 0; i < size; i++) val[i] = ptr[i];
    return val;
  }

  KOKKOS_FUNCTION
  auto operator=(const Bar& val) {
    for (size_t i = 0; i < size; i++) ptr[i] = val[i];
  }
};

// A TestAccessor mimicking some of what Sacado does
// Specifically:
//   * returns a proxy reference
//   * the underlying storage is of some basic scalar type
//   * the size of the allocation does not actually come from
//   sizeof(ElementType)
template <class ElementType, class MemorySpace>
struct TestAccessor {
  static_assert(std::is_same_v<std::remove_cv_t<ElementType>, Bar>);
  using value_type =
      std::conditional_t<std::is_const_v<ElementType>, const double, double>;

  using element_type = ElementType;
  using reference    = BarRef<value_type>;
  using data_handle_type =
      Kokkos::Impl::ReferenceCountedDataHandle<value_type, MemorySpace>;
  using offset_policy = TestAccessor;

  // View expects this from accessors right now
  using memory_space = MemorySpace;

  KOKKOS_DEFAULTED_FUNCTION
  constexpr TestAccessor() = default;

  template <
      class OtherElementType, class OtherMemorySpace,
      std::enable_if_t<std::is_constructible_v<element_type, OtherElementType>,
                       int> = 0>
  KOKKOS_FUNCTION constexpr TestAccessor(
      const TestAccessor<OtherElementType, OtherMemorySpace>& other) noexcept
      : size(other.size) {}

  KOKKOS_FUNCTION
  TestAccessor(const size_t val) : size(val) {}

  KOKKOS_FUNCTION
  constexpr reference access(
#ifndef KOKKOS_ENABLE_OPENACC
      const data_handle_type& p,
#else
      // FIXME OPENACC: illegal address when passing by reference
      data_handle_type p,
#endif
      size_t i) const noexcept {
    return reference{(p.get() + i * size), size};
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
    return p + i * size;
  }

  size_t size{0lu};
};

// Use the customization point to inject the custom accessor
template <class LayoutType, class DeviceType, class MemoryTraits>
constexpr auto customize_view_arguments(
    Kokkos::Impl::ViewArguments<Bar, LayoutType, DeviceType, MemoryTraits>) {
  return Kokkos::Impl::ViewCustomArguments<
      size_t, TestAccessor<Bar, typename DeviceType::memory_space>>{};
}

template <class LayoutType, class DeviceType, class MemoryTraits>
constexpr auto customize_view_arguments(
    Kokkos::Impl::ViewArguments<const Bar, LayoutType, DeviceType,
                                MemoryTraits>) {
  return Kokkos::Impl::ViewCustomArguments<
      size_t, TestAccessor<const Bar, typename DeviceType::memory_space>>{};
}

// Customization point to compute allocation sizes
template <class MappingType, class ElementType, class MemorySpace>
size_t allocation_size_from_mapping_and_accessor(
    const MappingType& map, const TestAccessor<ElementType, MemorySpace>& acc) {
  return map.required_span_size() * acc.size;
}
}  // namespace Foo

void test_allocation_type() {
  size_t ext0 = 10;
  size_t ext1 = 11;
  size_t size = 5;

  using view_t = Kokkos::View<Foo::Bar**, TEST_EXECSPACE>;

  // Make sure I got the accessor I expect
  static_assert(
      std::is_same_v<
          view_t::accessor_type,
          Foo::TestAccessor<Foo::Bar, typename TEST_EXECSPACE::memory_space>>);
  static_assert(std::is_same_v<view_t::pointer_type, double*>);

  using c_view_t = Kokkos::View<const Foo::Bar**, TEST_EXECSPACE>;
  static_assert(
      std::is_same_v<c_view_t::accessor_type,
                     Foo::TestAccessor<const Foo::Bar,
                                       typename TEST_EXECSPACE::memory_space>>);
  static_assert(std::is_same_v<c_view_t::pointer_type, const double*>);

  // accessor will be constructed from AccessorArg_t
  view_t a(Kokkos::view_alloc("A", Kokkos::Impl::AccessorArg_t{size}), ext0,
           ext1);
  ASSERT_EQ(a.accessor().size, size);
  static_assert(std::is_same_v<decltype(a.data()), double*>);

  // Test copy ctor to make sure the customize_view_arguments thing doesn't
  // interfere
  view_t a_copy = a;
  ASSERT_EQ(a_copy.accessor().size, size);

  c_view_t const_a = a;
  ASSERT_EQ(const_a.accessor().size, size);
  ASSERT_EQ(const_a.data(), a.data());

  view_t b(Kokkos::view_wrap(a.data(), Kokkos::Impl::AccessorArg_t{size}), ext0,
           ext1);
  ASSERT_EQ(b.accessor().size, size);

  // Get a compatible mapping for address calculation in the kernel
  using mapping_t = typename Kokkos::View<
      int**, typename TEST_EXECSPACE::memory_space>::mdspan_type::mapping_type;
  mapping_t map(Kokkos::dextents<size_t, 2>{ext0, ext1});

  // Test unmanaged ctors on GPU too (if GPU is enabled)
  int num_error = 0;
  Kokkos::parallel_reduce(
      "test_accessor_arg",
      Kokkos::MDRangePolicy<Kokkos::Rank<2>, TEST_EXECSPACE>({0, 0},
                                                             {ext0, ext1}),
      KOKKOS_LAMBDA(int i, int j, int& errors) {
        view_t c(Kokkos::view_wrap(a.data(), Kokkos::Impl::AccessorArg_t{size}),
                 ext0, ext1);
        if (c.accessor().size != size) errors++;
        // Test copy ctor to make sure the customize_view_arguments thing
        // doesn't interfere
        view_t c_copy = c;
        if (c_copy.accessor().size != size) errors++;
        for (size_t k = 0; k < size; k++) {
          if (&a(i, j)[k] != a.data() + map(i, j) * size + k) errors++;
          if (&c_copy(i, j)[k] != a.data() + map(i, j) * size + k) errors++;
        }
      },
      num_error);
  ASSERT_EQ(num_error, 0);
}

TEST(TEST_CATEGORY, view_customization_allocation_type) {
  test_allocation_type();
}

void test_create_mirror() {
  size_t ext0 = 10;
  size_t ext1 = 11;
  size_t size = 5;

  using view_t = Kokkos::View<Foo::Bar**, TEST_EXECSPACE>;

  view_t a(Kokkos::view_alloc("A", Kokkos::Impl::AccessorArg_t{size}), ext0,
           ext1);

  auto host_a1 = Kokkos::create_mirror_view(a);
  ASSERT_EQ(host_a1.accessor().size, size);
  static_assert(std::is_same_v<decltype(host_a1.data()), double*>);

  auto host_a2 = Kokkos::create_mirror(a);
  ASSERT_EQ(host_a2.accessor().size, size);
  static_assert(std::is_same_v<decltype(host_a2.data()), double*>);
}

TEST(TEST_CATEGORY, view_customization_mirror) { test_create_mirror(); }

void test_deep_copy() {
  size_t ext0 = 10;
  size_t ext1 = 11;
  size_t size = 5;

  using view_t = Kokkos::View<Foo::Bar**, TEST_EXECSPACE>;

  view_t a(Kokkos::view_alloc("A", Kokkos::Impl::AccessorArg_t{size}), ext0,
           ext1);

  auto host_a = Kokkos::create_mirror_view(a);

  for (size_t i = 0; i < ext0; i++)
    for (size_t j = 0; j < ext1; j++)
      for (size_t k = 0; k < size; k++) host_a(i, j)[k] = i * 100 + j + 0.1 * k;

  // Contiguous deep_copy, potentially host to device
  Kokkos::deep_copy(a, host_a);

  Kokkos::MDRangePolicy<Kokkos::Rank<2>, TEST_EXECSPACE> policy2d({0, 0},
                                                                  {ext0, ext1});

  int num_errors = 0;
  Kokkos::parallel_reduce(
      "view_customization_deep_copy_a", policy2d,
      KOKKOS_LAMBDA(int i, int j, int& error) {
        for (size_t k = 0; k < size; k++) {
          if (Kokkos::abs(a(i, j)[k] - (i * 100 + j + 0.1 * k)) >
              Kokkos::Experimental::epsilon_v<float>)
            error++;
          a(i, j)[k] = i * 200 + j + 0.1 * k;
        }
      },
      num_errors);
  ASSERT_EQ(num_errors, 0);

  // Contiguous deep_copy, potentially device to host
  Kokkos::deep_copy(host_a, a);
  for (size_t i = 0; i < ext0; i++)
    for (size_t j = 0; j < ext1; j++)
      for (size_t k = 0; k < size; k++)
        ASSERT_FLOAT_EQ(host_a(i, j)[k], i * 200 + j + 0.1 * k);
  auto b = Kokkos::create_mirror(TEST_EXECSPACE::memory_space(), a);

  num_errors = 0;
  Kokkos::parallel_reduce(
      "view_customization_check_pre_deep_copy_b", policy2d,
      KOKKOS_LAMBDA(int i, int j, int& error) {
        for (size_t k = 0; k < size; k++) {
          // Note b is value initalized for the scalar value type for the
          // allocation
          if (b(i, j)[k] != 0) {
            error++;
          }
        }
      },
      num_errors);
  ASSERT_EQ(num_errors, 0);

  // Contiguous deep_copy, same execspace
  Kokkos::deep_copy(b, a);

  num_errors = 0;
  Kokkos::parallel_reduce(
      "view_customization_check_pre_deep_copy_b", policy2d,
      KOKKOS_LAMBDA(int i, int j, int& error) {
        for (size_t k = 0; k < size; k++) {
          if (Kokkos::abs(b(i, j)[k] - (i * 200 + j + 0.1 * k)) >
              Kokkos::Experimental::epsilon_v<float>)
            error++;
        }
      },
      num_errors);
  ASSERT_EQ(num_errors, 0);
}

TEST(TEST_CATEGORY, view_customization_deep_copy) { test_deep_copy(); }
