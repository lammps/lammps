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

// ElementType which is actually just a tag type
struct BarStrided {};

// Reference type for Bar
struct BarRefStrided {
  double* ptr{nullptr};
  size_t size{0ul};
  size_t stride{0ul};
  KOKKOS_FUNCTION
  double& operator[](size_t idx) const { return ptr[idx * stride]; }
};

// A TestAccessor mimicking some of what Sacado does
// Specifically:
//   * returns a proxy reference
//   * the underlying storage is of some basic scalar type
//   * the size of the allocation does not actually come from
//   sizeof(ElementType)
template <class ElementType, class MemorySpace>
struct TestAccessorStrided {
  static_assert(std::is_same_v<std::remove_cv_t<ElementType>, BarStrided>);
  using element_type     = ElementType;
  using reference        = BarRefStrided;
  using data_handle_type = Kokkos::Impl::ReferenceCountedDataHandle<
      std::conditional_t<std::is_const_v<ElementType>, const double, double>,
      MemorySpace>;
  using offset_policy = TestAccessorStrided;

  // View expects this from accessors right now
  using memory_space = MemorySpace;

  KOKKOS_DEFAULTED_FUNCTION
  constexpr TestAccessorStrided() = default;

  template <class OtherElementType,
            std::enable_if_t<std::is_constructible_v<
                                 Kokkos::default_accessor<element_type>,
                                 Kokkos::default_accessor<OtherElementType>>,
                             int> = 0>
  KOKKOS_FUNCTION constexpr TestAccessorStrided(
      const TestAccessorStrided<OtherElementType, MemorySpace>& other) noexcept
      : size(other.size), stride(other.stride) {}

  KOKKOS_FUNCTION
  TestAccessorStrided(const size_t val, const size_t s_val)
      : size(val), stride(s_val) {}

  KOKKOS_FUNCTION
  constexpr reference access(
#ifndef KOKKOS_ENABLE_OPENACC
      const data_handle_type& p,
#else
      // FIXME OPENACC: illegal address when passing by reference
      data_handle_type p,
#endif
      size_t i) const noexcept {
    return BarRefStrided{(p.get() + i), size, stride};
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

  size_t size{0lu};
  size_t stride{0lu};
};

// Use the customization point to inject the custom accessor
template <class LayoutType, class DeviceType, class MemoryTraits>
constexpr auto customize_view_arguments(
    Kokkos::Impl::ViewArguments<BarStrided, LayoutType, DeviceType,
                                MemoryTraits>) {
  return Kokkos::Impl::ViewCustomArguments<
      size_t,
      TestAccessorStrided<BarStrided, typename DeviceType::memory_space>>{};
}

template <class LayoutType, class DeviceType, class MemoryTraits>
constexpr auto customize_view_arguments(
    Kokkos::Impl::ViewArguments<const BarStrided, LayoutType, DeviceType,
                                MemoryTraits>) {
  return Kokkos::Impl::ViewCustomArguments<
      size_t, TestAccessorStrided<const BarStrided,
                                  typename DeviceType::memory_space>>{};
}

// Customization point to compute allocation sizes
template <class MappingType, class ElementType, class MemorySpace>
KOKKOS_INLINE_FUNCTION size_t allocation_size_from_mapping_and_accessor(
    const MappingType& map,
    const TestAccessorStrided<ElementType, MemorySpace>& acc) {
  return map.required_span_size() * acc.size;
}

// Customization point to create accessor from AccessorArg and Mapping
template <class ElementType, class MemorySpace, class MappingType>
KOKKOS_INLINE_FUNCTION constexpr auto accessor_from_mapping_and_accessor_arg(
    Kokkos::Impl::AccessorTypeTag<
        TestAccessorStrided<ElementType, MemorySpace>>,
    const MappingType& map, const Kokkos::Impl::AccessorArg_t& acc_arg) {
  return TestAccessorStrided<ElementType, MemorySpace>{
      acc_arg.value, map.required_span_size()};
}
}  // namespace Foo

void test_accessor_from_mapping() {
  size_t ext0   = 10;
  size_t ext1   = 11;
  size_t size   = 5;
  size_t stride = ext0 * ext1;

  using view_t = Kokkos::View<Foo::BarStrided**, TEST_EXECSPACE>;

  // Make sure I got the accessor I expect
  static_assert(
      std::is_same_v<
          view_t::accessor_type,
          Foo::TestAccessorStrided<Foo::BarStrided,
                                   typename TEST_EXECSPACE::memory_space>>);
  static_assert(std::is_same_v<view_t::pointer_type, double*>);

  using c_view_t = Kokkos::View<const Foo::BarStrided**, TEST_EXECSPACE>;
  static_assert(
      std::is_same_v<
          c_view_t::accessor_type,
          Foo::TestAccessorStrided<const Foo::BarStrided,
                                   typename TEST_EXECSPACE::memory_space>>);
  static_assert(std::is_same_v<c_view_t::pointer_type, const double*>);

  // accessor will be constructed from AccessorArg_t
  view_t a(Kokkos::view_alloc("A", Kokkos::Impl::AccessorArg_t{size}), ext0,
           ext1);
  ASSERT_EQ(a.accessor().size, size);
  ASSERT_EQ(a.accessor().stride, stride);
  static_assert(std::is_same_v<decltype(a.data()), double*>);

  // Test copy ctor to make sure the customize_view_arguments thing doesn't
  // interfere
  view_t a_copy = a;
  ASSERT_EQ(a_copy.accessor().size, size);
  ASSERT_EQ(a_copy.accessor().stride, stride);

  c_view_t const_a = a;
  ASSERT_EQ(const_a.accessor().size, size);
  ASSERT_EQ(const_a.accessor().stride, stride);
  ASSERT_EQ(const_a.data(), a.data());

  view_t b(Kokkos::view_wrap(a.data(), Kokkos::Impl::AccessorArg_t{size}), ext0,
           ext1);
  ASSERT_EQ(b.accessor().size, size);
  ASSERT_EQ(b.accessor().stride, stride);

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
        if (c.accessor().stride != stride) errors++;
        // Test copy ctor to make sure the customize_view_arguments thing
        // doesn't interfere
        view_t c_copy = c;
        if (c_copy.accessor().size != size) errors++;
        if (c_copy.accessor().stride != stride) errors++;
        for (size_t k = 0; k < size; k++) {
          if (&a(i, j)[k] != a.data() + map(i, j) + k * stride) errors++;
          if (&c_copy(i, j)[k] != a.data() + map(i, j) + k * stride) errors++;
        }
      },
      num_error);
  ASSERT_EQ(num_error, 0);
}

TEST(TEST_CATEGORY, view_customization_accessor_from_mapping) {
  test_accessor_from_mapping();
}

// This tests the ability to pass in the accessor arg
// as an additional integral argument to constructor and shmem_size

template <class ViewType, size_t... Idx>
bool test_device_side_ctor(ViewType a, std::index_sequence<Idx...>) {
  int error = 0;
  Kokkos::parallel_reduce(
      Kokkos::RangePolicy<TEST_EXECSPACE>(0, 1),
      KOKKOS_LAMBDA(int, int& errors) {
        ViewType a2(a.data(), a.extent(Idx)..., a.accessor().size);
        if (a2.accessor().size != a.accessor().size) errors++;
        if (a2.accessor().stride != a.accessor().stride) errors++;
      },
      error);
  return error == 0;
}

TEST(TEST_CATEGORY, view_customization_extra_int_arg) {
  // Rank 0
  {
    using view_t = Kokkos::View<Foo::BarStrided, TEST_EXECSPACE>;
    view_t a("A", 5);
    ASSERT_EQ(a.accessor().size, size_t(5));
    ASSERT_EQ(a.accessor().stride, size_t(1));
    view_t b(a.data(), 5);
    ASSERT_EQ(b.accessor().size, size_t(5));
    ASSERT_EQ(b.accessor().stride, size_t(1));
    size_t shmem               = view_t::shmem_size(5);
    size_t expected_shmem_size = 5lu * sizeof(double) + sizeof(double);
    ASSERT_EQ(shmem, expected_shmem_size);
    ASSERT_TRUE(test_device_side_ctor(a, std::make_index_sequence<0>()));
  }
  // Rank 3
  {
    using view_t = Kokkos::View<Foo::BarStrided***, TEST_EXECSPACE>;
    view_t a("A", 3, 7, 11, 5);
    ASSERT_EQ(a.accessor().size, size_t(5));
    ASSERT_EQ(a.accessor().stride, size_t(3 * 7 * 11));
    view_t b(a.data(), 3, 7, 11, 5);
    ASSERT_EQ(b.accessor().size, size_t(5));
    ASSERT_EQ(b.accessor().stride, size_t(3 * 7 * 11));
    size_t shmem = view_t::shmem_size(3, 7, 11, 5);
    size_t expected_shmem_size =
        3lu * 7lu * 11lu * 5lu * sizeof(double) + sizeof(double);
    ASSERT_EQ(shmem, expected_shmem_size);
    ASSERT_TRUE(test_device_side_ctor(a, std::make_index_sequence<3>()));
  }
  // Rank 6
  {
    using view_t = Kokkos::View<Foo::BarStrided******, TEST_EXECSPACE>;
    view_t a("A", 2, 3, 2, 7, 2, 11, 5);
    ASSERT_EQ(a.accessor().size, size_t(5));
    ASSERT_EQ(a.accessor().stride, size_t(8 * 3 * 7 * 11));
    view_t b(a.data(), 2, 3, 2, 7, 2, 11, 5);
    ASSERT_EQ(b.accessor().size, size_t(5));
    ASSERT_EQ(b.accessor().stride, size_t(8 * 3 * 7 * 11));
    size_t shmem = view_t::shmem_size(2, 3, 2, 7, 2, 11, 5);
    size_t expected_shmem_size =
        2lu * 3lu * 2lu * 7lu * 2lu * 11lu * 5lu * sizeof(double) +
        sizeof(double);
    ASSERT_EQ(shmem, expected_shmem_size);
    ASSERT_TRUE(test_device_side_ctor(a, std::make_index_sequence<6>()));
  }
  // Rank 3
  {
    using view_t = Kokkos::View<Foo::BarStrided** [11], TEST_EXECSPACE>;
    view_t a("A", 3, 7, 11, 5);
    ASSERT_EQ(a.accessor().size, size_t(5));
    ASSERT_EQ(a.accessor().stride, size_t(3 * 7 * 11));
    view_t b(a.data(), 3, 7, 11, 5);
    ASSERT_EQ(b.accessor().size, size_t(5));
    ASSERT_EQ(b.accessor().stride, size_t(3 * 7 * 11));
    size_t shmem = view_t::shmem_size(3, 7, 5);
    size_t expected_shmem_size =
        3lu * 7lu * 11lu * 5lu * sizeof(double) + sizeof(double);
    ASSERT_EQ(shmem, expected_shmem_size);
    ASSERT_TRUE(test_device_side_ctor(a, std::make_index_sequence<3>()));
  }
  // Rank 6
  {
    using view_t = Kokkos::View<Foo::BarStrided***** [11], TEST_EXECSPACE>;
    view_t a("A", 2, 3, 2, 7, 2, 11, 5);
    ASSERT_EQ(a.accessor().size, size_t(5));
    ASSERT_EQ(a.accessor().stride, size_t(8 * 3 * 7 * 11));
    view_t b(a.data(), 2, 3, 2, 7, 2, 11, 5);
    ASSERT_EQ(b.accessor().size, size_t(5));
    ASSERT_EQ(b.accessor().stride, size_t(8 * 3 * 7 * 11));
    size_t shmem = view_t::shmem_size(2, 3, 2, 7, 2, 5);
    size_t expected_shmem_size =
        2lu * 3lu * 2lu * 7lu * 2lu * 11lu * 5lu * sizeof(double) +
        sizeof(double);
    ASSERT_EQ(shmem, expected_shmem_size);
    ASSERT_TRUE(test_device_side_ctor(a, std::make_index_sequence<6>()));
  }
}

void test_scratch_memory_allocation() {
  using view_t = Kokkos::View<Foo::BarStrided*, TEST_EXECSPACE>;
  size_t size  = 5;

  using policy_t = Kokkos::TeamPolicy<TEST_EXECSPACE>;
  policy_t p(1, Kokkos::AUTO());

  using team_t = typename policy_t::member_type;

  size_t shmem = view_t::shmem_size(2, size);
  p.set_scratch_size(0, Kokkos::PerTeam(3 * shmem), Kokkos::PerThread(shmem));
  int reported_errors;
  Kokkos::parallel_reduce(
      "TestScratch", p,
      KOKKOS_LAMBDA(const team_t& team, int& errors) {
        view_t tshmem1(team.team_shmem(), 2, size);
        view_t tshmem2(team.team_scratch(0), 2, size);
        if (tshmem1.accessor().size != size) errors++;
        if (tshmem2.accessor().size != size) errors++;
        if (tshmem1.data() + 2 * size != tshmem2.data()) errors++;

        view_t my_data(team.thread_scratch(0), 2, size);
        if (my_data.accessor().size != size) errors++;
        view_t tshmem3(team.team_scratch(0), 2, size);
        if (my_data.data() + 2 * size > tshmem3.data()) errors++;
      },
      reported_errors);
  ASSERT_EQ(reported_errors, 0);
}

TEST(TEST_CATEGORY, view_customization_scratch_memory) {
#ifdef KOKKOS_ENABLE_HPX  // FIXME_HPX
  if (std::is_same_v<Kokkos::Experimental::HPX, TEST_EXECSPACE>)
    GTEST_SKIP() << "HPX backend fails this test intermittently, since its the "
                    "only backend failing disabling the test for now";
#endif
  test_scratch_memory_allocation();
}
