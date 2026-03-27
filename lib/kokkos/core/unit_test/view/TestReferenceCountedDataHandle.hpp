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

namespace {
using element_t = float;
using mem_t     = typename TEST_EXECSPACE::memory_space;
using data_handle_t =
    Kokkos::Impl::ReferenceCountedDataHandle<element_t, mem_t>;
using const_data_handle_t =
    Kokkos::Impl::ReferenceCountedDataHandle<const element_t, mem_t>;
using data_handle_anonym_t =
    Kokkos::Impl::ReferenceCountedDataHandle<element_t, Kokkos::AnonymousSpace>;
using const_data_handle_anonym_t =
    Kokkos::Impl::ReferenceCountedDataHandle<const element_t,
                                             Kokkos::AnonymousSpace>;

}  // namespace

TEST(TEST_CATEGORY, RefCountedDataHandle_Typedefs) {
  static_assert(std::is_same_v<data_handle_t::value_type, element_t>);
  static_assert(std::is_same_v<data_handle_t::pointer, element_t*>);
  static_assert(std::is_same_v<data_handle_t::reference, element_t&>);
  static_assert(std::is_same_v<data_handle_t::memory_space, mem_t>);
}

template <class DataHandleType, class ConstDataHandleType>
void test_ref_counted_data_handle() {
  auto shared_alloc =
      Kokkos::Impl::make_shared_allocation_record<element_t, mem_t,
                                                  TEST_EXECSPACE>(
          100, "Test", mem_t(), std::optional<TEST_EXECSPACE>(std::nullopt),
          std::bool_constant<true>(),    // init
          std::bool_constant<false>());  // sequential_host_init

  element_t* ptr         = static_cast<element_t*>(shared_alloc->data());
  const element_t* c_ptr = ptr;
  DataHandleType dh(shared_alloc);
  ASSERT_EQ(dh.use_count(), 1);
  ASSERT_EQ(dh.get_label(), std::string("Test"));
  ASSERT_EQ(dh.get(), ptr);
  ASSERT_EQ(dh.has_record(), true);
  {
    element_t* ptr_tmp(dh);
    ASSERT_EQ(ptr_tmp, ptr);
    static_assert(!std::is_convertible_v<data_handle_t, element_t*>);
  }
  {
    ConstDataHandleType c_dh(dh);
    ASSERT_EQ(dh.use_count(), 2);
    ASSERT_EQ(c_dh.use_count(), 2);
  }
  ASSERT_EQ(dh.use_count(), 1);

  DataHandleType um_dh(ptr);
  ASSERT_EQ(um_dh.get(), ptr);
  ASSERT_EQ(um_dh.has_record(), false);

  DataHandleType dh_offset(dh, ptr + 5);
  ASSERT_EQ(dh_offset.use_count(), 2);
  ASSERT_EQ(dh_offset.get(), ptr + 5);
  ASSERT_EQ(dh_offset.get_label(), std::string("Test"));
  ASSERT_EQ(dh_offset.has_record(), true);
  {
    element_t* ptr_tmp(dh_offset);
    ASSERT_EQ(ptr_tmp, ptr + 5);
  }
  Kokkos::View<int, TEST_EXECSPACE> errors("Errors");

  // clang-format screws the following pieces up for some reason
  // Tested with 16 and with 18 to the same effect
  // clang-format off
  Kokkos::parallel_for(Kokkos::RangePolicy<TEST_EXECSPACE>(0, 1), KOKKOS_LAMBDA(int) {

    // default ctor and non-const to const
    {
      DataHandleType dh2(dh);
      if(dh2.get() != ptr) errors() += 1;
      ConstDataHandleType c_dh2(dh);
      ConstDataHandleType c_dh3(c_dh2);
      static_assert(!std::is_constructible_v<data_handle_t, const_data_handle_t>);
    }

    {
      // ctor from pointer
      DataHandleType dh2(ptr);
      if (dh2.get() != ptr) errors() += 2;
      ConstDataHandleType c_dh1(ptr);
      if (c_dh1.get() != ptr) errors() += 4;
      ConstDataHandleType c_dh2(c_ptr);
      if (c_dh2.get() != ptr) errors() += 8;
      static_assert(!std::is_constructible_v<data_handle_t, decltype(c_ptr)>);
    }

    // ctor for subviews
    {
      DataHandleType dh2(dh, ptr + 5);
      if (dh2.get() != ptr + 5) errors() += 16;
    }
  });

  int h_errors = 0;
  Kokkos::deep_copy(h_errors, errors);
  ASSERT_FALSE(h_errors & 1);
  ASSERT_FALSE(h_errors & 2);
  ASSERT_FALSE(h_errors & 4);
  ASSERT_FALSE(h_errors & 8);
  ASSERT_FALSE(h_errors & 16);
}
// clang-format on

TEST(TEST_CATEGORY, RefCountedDataHandle) {
  test_ref_counted_data_handle<data_handle_t, const_data_handle_t>();
}

TEST(TEST_CATEGORY, RefCountedDataHandleAnonym) {
  test_ref_counted_data_handle<data_handle_anonym_t,
                               const_data_handle_anonym_t>();
}

template <class T>
KOKKOS_FUNCTION void unused_variable_sink(T) {}

void test_ref_counted_data_handle_conversion() {
  auto shared_alloc1 =
      Kokkos::Impl::make_shared_allocation_record<element_t, mem_t,
                                                  TEST_EXECSPACE>(
          100, "Test1", mem_t(), std::optional<TEST_EXECSPACE>(std::nullopt),
          std::bool_constant<true>(),    // init
          std::bool_constant<false>());  // sequential_host_init

  element_t* ptr1         = static_cast<element_t*>(shared_alloc1->data());
  const element_t* c_ptr1 = ptr1;
  unused_variable_sink(c_ptr1);

  data_handle_t dh(shared_alloc1);
  ASSERT_EQ(dh.use_count(), 1);
  ASSERT_EQ(dh.get_label(), std::string("Test1"));
  ASSERT_EQ(dh.get(), ptr1);
  ASSERT_EQ(dh.has_record(), true);

  auto shared_alloc2 =
      Kokkos::Impl::make_shared_allocation_record<element_t, mem_t,
                                                  TEST_EXECSPACE>(
          100, "Test2", mem_t(), std::optional<TEST_EXECSPACE>(std::nullopt),
          std::bool_constant<true>(),    // init
          std::bool_constant<false>());  // sequential_host_init

  element_t* ptr2         = static_cast<element_t*>(shared_alloc2->data());
  const element_t* c_ptr2 = ptr2;
  unused_variable_sink(c_ptr2);

  data_handle_anonym_t dha(shared_alloc2);
  ASSERT_EQ(dha.use_count(), 1);
  ASSERT_EQ(dha.get_label(), std::string("Test2"));
  ASSERT_EQ(dha.get(), ptr2);
  ASSERT_EQ(dha.has_record(), true);

  {
    data_handle_anonym_t dha2(dh);
    ASSERT_EQ(dha2.use_count(), 2);
    ASSERT_EQ(dha2.get_label(), std::string("Test1"));
    ASSERT_EQ(dha2.get(), ptr1);
    ASSERT_EQ(dha2.has_record(), true);

    data_handle_t dh2(dha);
    ASSERT_EQ(dh2.use_count(), 2);
    ASSERT_EQ(dh2.get_label(), std::string("Test2"));
    ASSERT_EQ(dh2.get(), ptr2);
    ASSERT_EQ(dh2.has_record(), true);

    dha2 = dh2;
    ASSERT_EQ(dha2.use_count(), 3);
    ASSERT_EQ(dha2.get_label(), std::string("Test2"));
    ASSERT_EQ(dha2.get(), ptr2);
    ASSERT_EQ(dha2.has_record(), true);
  }

  ASSERT_EQ(dh.use_count(), 1);
  ASSERT_EQ(dh.get_label(), std::string("Test1"));
  ASSERT_EQ(dh.get(), ptr1);
  ASSERT_EQ(dh.has_record(), true);

  ASSERT_EQ(dha.use_count(), 1);
  ASSERT_EQ(dha.get_label(), std::string("Test2"));
  ASSERT_EQ(dha.get(), ptr2);
  ASSERT_EQ(dha.has_record(), true);
}

TEST(TEST_CATEGORY, RefCountedDataHandleConversion) {
  test_ref_counted_data_handle_conversion();
}
