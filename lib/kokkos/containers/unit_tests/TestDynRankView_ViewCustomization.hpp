// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif
#include <Kokkos_DynRankView.hpp>

// Duplicate from
// core/unit_test/view/TestViewCustomizationAccessorFromMapping.hpp
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

  template <class OtherElementType, class OtherMemorySpace,
            std::enable_if_t<std::is_constructible_v<
                                 Kokkos::default_accessor<element_type>,
                                 Kokkos::default_accessor<OtherElementType>>,
                             int> = 0>
  KOKKOS_FUNCTION constexpr TestAccessorStrided(
      const TestAccessorStrided<OtherElementType, OtherMemorySpace>&
          other) noexcept
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

// This tests the ability to pass in the accessor arg
// as an additional integral argument to constructor and shmem_size
TEST(TEST_CATEGORY, view_customization_extra_int_arg) {
  using view_t = Kokkos::DynRankView<Foo::BarStrided, TEST_EXECSPACE>;
  // Rank 0
  {
    view_t a("A", 5);
    ASSERT_EQ(a.rank(), 0lu);
    ASSERT_EQ(a.accessor().size, size_t(5));
    ASSERT_EQ(a.accessor().stride, size_t(1));

    view_t b(a.data(), 5);
    ASSERT_EQ(b.rank(), 0lu);
    ASSERT_EQ(b.accessor().size, 5lu);
    ASSERT_EQ(b.accessor().stride, 1lu);
    size_t shmem               = view_t::shmem_size(5);
    size_t expected_shmem_size = 5lu * sizeof(double) + sizeof(double);
    ASSERT_EQ(shmem, expected_shmem_size);
  }
  // Rank 3
  {
    view_t a("A", 3, 7, 11, 5);
    ASSERT_EQ(a.rank(), 3lu);
    ASSERT_EQ(a.extent(0), 3lu);
    ASSERT_EQ(a.extent(1), 7lu);
    ASSERT_EQ(a.extent(2), 11lu);
    ASSERT_EQ(a.accessor().size, 5lu);
    ASSERT_EQ(a.accessor().stride, size_t(3 * 7 * 11));

    view_t b(a.data(), 3, 7, 11, 5);
    ASSERT_EQ(b.rank(), 3lu);
    ASSERT_EQ(b.extent(0), 3lu);
    ASSERT_EQ(b.extent(1), 7lu);
    ASSERT_EQ(b.extent(2), 11lu);
    ASSERT_EQ(b.accessor().size, 5lu);
    ASSERT_EQ(b.accessor().stride, size_t(3 * 7 * 11));
    size_t shmem = view_t::shmem_size(3, 7, 11, 5);
    size_t expected_shmem_size =
        3lu * 7lu * 11lu * 5lu * sizeof(double) + sizeof(double);
    ASSERT_EQ(shmem, expected_shmem_size);
  }
  // Rank 6
  {
    view_t a("A", 2, 3, 2, 7, 2, 11, 5);
    ASSERT_EQ(a.rank(), 6lu);
    ASSERT_EQ(a.extent(0), 2lu);
    ASSERT_EQ(a.extent(1), 3lu);
    ASSERT_EQ(a.extent(2), 2lu);
    ASSERT_EQ(a.extent(3), 7lu);
    ASSERT_EQ(a.extent(4), 2lu);
    ASSERT_EQ(a.extent(5), 11lu);
    ASSERT_EQ(a.accessor().size, 5lu);
    ASSERT_EQ(a.accessor().stride, size_t(8 * 3 * 7 * 11));
    view_t b(a.data(), 2, 3, 2, 7, 2, 11, 5);
    ASSERT_EQ(b.rank(), 6lu);
    ASSERT_EQ(b.extent(0), 2lu);
    ASSERT_EQ(b.extent(1), 3lu);
    ASSERT_EQ(b.extent(2), 2lu);
    ASSERT_EQ(b.extent(3), 7lu);
    ASSERT_EQ(b.extent(4), 2lu);
    ASSERT_EQ(b.extent(5), 11lu);
    ASSERT_EQ(b.accessor().size, 5lu);
    ASSERT_EQ(b.accessor().stride, size_t(8 * 3 * 7 * 11));
    size_t shmem = view_t::shmem_size(2, 3, 2, 7, 2, 11, 5);
    size_t expected_shmem_size =
        2lu * 3lu * 2lu * 7lu * 2lu * 11lu * 5lu * sizeof(double) +
        sizeof(double);
    ASSERT_EQ(shmem, expected_shmem_size);
  }
  // Rank 7
  {
    view_t a("A", 2, 3, 2, 7, 2, 11, 2, 5);
    ASSERT_EQ(a.rank(), 7lu);
    ASSERT_EQ(a.extent(0), 2lu);
    ASSERT_EQ(a.extent(1), 3lu);
    ASSERT_EQ(a.extent(2), 2lu);
    ASSERT_EQ(a.extent(3), 7lu);
    ASSERT_EQ(a.extent(4), 2lu);
    ASSERT_EQ(a.extent(5), 11lu);
    ASSERT_EQ(a.extent(2), 2lu);
    ASSERT_EQ(a.accessor().size, 5lu);
    ASSERT_EQ(a.accessor().stride, size_t(16 * 3 * 7 * 11));
    view_t b(a.data(), 2, 3, 2, 7, 2, 11, 2, 5);
    ASSERT_EQ(b.rank(), 7lu);
    ASSERT_EQ(b.extent(0), 2lu);
    ASSERT_EQ(b.extent(1), 3lu);
    ASSERT_EQ(b.extent(2), 2lu);
    ASSERT_EQ(b.extent(3), 7lu);
    ASSERT_EQ(b.extent(4), 2lu);
    ASSERT_EQ(b.extent(5), 11lu);
    ASSERT_EQ(b.extent(6), 2lu);
    ASSERT_EQ(b.accessor().size, 5lu);
    ASSERT_EQ(b.accessor().stride, size_t(16 * 3 * 7 * 11));
    size_t shmem = view_t::shmem_size(2, 3, 2, 7, 2, 11, 2, 5);
    size_t expected_shmem_size =
        2lu * 3lu * 2lu * 7lu * 2lu * 11lu * 2lu * 5lu * sizeof(double) +
        sizeof(double);
    ASSERT_EQ(shmem, expected_shmem_size);
  }
  // With accessor arg and no label
  {
    // This should not interpret the last argument (11) as an accessor arg since
    // we are providing AccessorArg_t explicitly
    // Note that with and without
    // labels are two separate cases
    view_t a(Kokkos::view_alloc(Kokkos::Impl::AccessorArg_t{5ul}), 3, 7, 11);
    ASSERT_EQ(a.rank(), 3lu);
    ASSERT_EQ(a.extent(0), 3lu);
    ASSERT_EQ(a.extent(1), 7lu);
    ASSERT_EQ(a.extent(2), 11lu);
    ASSERT_EQ(a.accessor().size, 5lu);
    ASSERT_EQ(a.accessor().stride, size_t(3 * 7 * 11));
  }
  // With accessor arg and label
  {
    // This should not interpret the last argument (11) as an accessor arg since
    // we are providing AccessorArg_t explicitly
    // Note that with and without
    // labels are two separate cases
    view_t a(Kokkos::view_alloc("A", Kokkos::Impl::AccessorArg_t{5ul}), 3, 7,
             11);
    ASSERT_EQ(a.rank(), 3lu);
    ASSERT_EQ(a.extent(0), 3lu);
    ASSERT_EQ(a.extent(1), 7lu);
    ASSERT_EQ(a.extent(2), 11lu);
    ASSERT_EQ(a.accessor().size, 5lu);
    ASSERT_EQ(a.accessor().stride, size_t(3 * 7 * 11));
  }
  // Create mirror
  {
    view_t a("A", 3, 7, 11, 5);

    auto b = Kokkos::create_mirror(a);
    ASSERT_EQ(b.rank(), 3lu);
    ASSERT_EQ(b.extent(0), 3lu);
    ASSERT_EQ(b.extent(1), 7lu);
    ASSERT_EQ(b.extent(2), 11lu);
    ASSERT_EQ(b.accessor().size, 5lu);
    ASSERT_EQ(b.accessor().stride, size_t(3 * 7 * 11));
  }
}
