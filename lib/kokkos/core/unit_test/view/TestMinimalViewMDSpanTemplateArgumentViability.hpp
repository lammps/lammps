// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <gtest/gtest.h>

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
import kokkos.core_impl;
#else
#include <Kokkos_Core.hpp>
#endif
#include <cstddef>

using view_t_64bit_idx =
    Kokkos::View<int, Kokkos::dextents<int64_t, 4>, Kokkos::layout_left,
                 Kokkos::Experimental::Accessor<int, Kokkos::HostSpace>>;
using view_t_32bit_idx =
    Kokkos::View<int, Kokkos::dextents<int32_t, 4>, Kokkos::layout_left,
                 Kokkos::Experimental::Accessor<int, Kokkos::HostSpace>>;

// Just a minimal check that the object storage actually changes
// Doing exact numbers would require a bunch of compiler ifdefs
// due to missing full support of no-unique-address in some cases
static_assert(sizeof(view_t_64bit_idx) > sizeof(view_t_32bit_idx));

template <class ElementType, class Exts,
          class LayoutType = typename Kokkos::View<int>::layout_type>
struct ViewTestHarness {
  // using the default accessor type
  using accessor_type =
      typename Kokkos::View<ElementType, TEST_EXECSPACE>::accessor_type;
  using new_view_t = Kokkos::View<ElementType, Exts, LayoutType, accessor_type>;

  // check that the index_type (and thus the template argument style) propagates
  // this is just a sanity check, the real check is in
  // TestViewMDSpanTemplateArgumentTypeDefs.cpp
  static_assert(std::is_same_v<typename new_view_t::type::index_type,
                               typename new_view_t::index_type>);
  static_assert(std::is_same_v<typename new_view_t::const_type::index_type,
                               typename new_view_t::index_type>);
  static_assert(std::is_same_v<typename new_view_t::non_const_type::index_type,
                               typename new_view_t::index_type>);
  static_assert(
      std::is_same_v<typename new_view_t::host_mirror_type::index_type,
                     typename new_view_t::index_type>);

  using old_view_t = Kokkos::View<
      typename new_view_t::data_type, typename new_view_t::array_layout,
      typename new_view_t::device_type, typename new_view_t::memory_traits>;
  using rank_indicies =
      decltype(std::make_index_sequence<new_view_t::rank()>());

  // function to get zeros in a fold expression
  // , operator and ternary result in warnings
  template <class T>
  static long unsigned make_zero(const T&) {
    return 0lu;
  }

  template <class ViewT, class... Extents>
  static void init_view(ViewT a, size_t extra_val, Extents... extents) {
    using exec_t = typename ViewT::execution_space;
    auto p = Kokkos::MDRangePolicy<exec_t, Kokkos::Rank<new_view_t::rank()>>(
        {(make_zero(extents))...}, {static_cast<long unsigned>(extents)...});
    Kokkos::parallel_for(
        p,
        KOKKOS_LAMBDA(Extents... idx) { a(idx...) = (idx + ... + extra_val); });
  }

  template <class ViewT, class... Extents>
  static size_t check_view(ViewT a, size_t extra_val, Extents... extents) {
    using exec_t = typename ViewT::execution_space;
    auto p = Kokkos::MDRangePolicy<exec_t, Kokkos::Rank<new_view_t::rank()>>(
        {(make_zero(extents))...}, {static_cast<long unsigned>(extents)...});
    size_t errors = 0;
    Kokkos::parallel_reduce(
        p,
        KOKKOS_LAMBDA(Extents... idx, size_t & lerr) {
          if (a(idx...) != static_cast<typename ViewT::element_type>(
                               (idx + ... + extra_val)))
            lerr++;
        },
        errors);
    return errors;
  }

  template <class... Extents>
  static void access(Extents... extents) {
    new_view_t a("A", extents...);

    old_view_t b(a);

    ASSERT_EQ(a.data(), b.data());
    ASSERT_EQ(a.extents(), b.extents());

    init_view(a, 1, extents...);
    size_t errors = check_view(b, 1, extents...);
    ASSERT_EQ(errors, 0lu);
  }

  template <class... Extents>
  static void deep_copy(Extents... extents) {
    new_view_t a("A", extents...);

    // New View to new View
    {
      new_view_t b("B", extents...);

      init_view(a, 1, extents...);
      Kokkos::deep_copy(b, a);
      size_t errors = check_view(b, 1, extents...);
      ASSERT_EQ(errors, 0lu);
    }

    // New View to old View
    {
      old_view_t b("B", extents...);

      init_view(a, 1, extents...);
      Kokkos::deep_copy(b, a);
      size_t errors_to = check_view(b, 1, extents...);
      ASSERT_EQ(errors_to, 0lu);

      init_view(b, 2, extents...);
      Kokkos::deep_copy(a, b);
      size_t errors_from = check_view(a, 2, extents...);
      ASSERT_EQ(errors_from, 0lu);
    }

    // New View to new View::host_mirror_type
    {
      typename new_view_t::host_mirror_type b("B", extents...);

      init_view(a, 1, extents...);
      Kokkos::deep_copy(b, a);
      size_t errors_to = check_view(b, 1, extents...);
      ASSERT_EQ(errors_to, 0lu);

      init_view(b, 2, extents...);
      Kokkos::deep_copy(a, b);
      size_t errors_from = check_view(a, 2, extents...);
      ASSERT_EQ(errors_from, 0lu);
    }

    // New View to old View::host_mirror_type
    {
      typename old_view_t::host_mirror_type b("B", extents...);

      init_view(a, 1, extents...);
      Kokkos::deep_copy(b, a);
      size_t errors_to = check_view(b, 1, extents...);
      ASSERT_EQ(errors_to, 0lu);

      init_view(b, 2, extents...);
      Kokkos::deep_copy(a, b);
      size_t errors_from = check_view(a, 2, extents...);
      ASSERT_EQ(errors_from, 0lu);
    }
  }

  // This checks that create_mirror works
  // correct return types are checked in
  // TestViewMDSpanTemplateArgumentTypeDefs.cpp
  template <class... Extents>
  static void create_mirror(Extents... extents) {
    new_view_t a("A", extents...);

    // Note: only for APU builds - i.e. making default GPU memory space
    // host accessible the two bools below are not the same in this test.
    constexpr bool is_host_space =
        std::is_same_v<typename new_view_t::memory_space, Kokkos::HostSpace>;
    constexpr bool host_accessible = Kokkos::SpaceAccessibility<
        Kokkos::DefaultHostExecutionSpace,
        typename new_view_t::memory_space>::accessible;

    {
      auto h_a = Kokkos::create_mirror(a);
      ASSERT_EQ(a.extents(), h_a.extents());
      ASSERT_EQ(a.data() == h_a.data(), false);
    }
    {
      typename new_view_t::const_type ac = a;
      auto h_ac                          = Kokkos::create_mirror(ac);
      ASSERT_EQ(a.extents(), h_ac.extents());
      ASSERT_EQ(a.data() == h_ac.data(), false);
    }
    {
      auto h_a = Kokkos::create_mirror(Kokkos::HostSpace(), a);
      ASSERT_EQ(a.extents(), h_a.extents());
      ASSERT_EQ(a.data() == h_a.data(), false);
    }
    {
      typename new_view_t::const_type ac = a;
      auto h_ac = Kokkos::create_mirror(Kokkos::HostSpace(), ac);
      ASSERT_EQ(a.extents(), h_ac.extents());
      ASSERT_EQ(a.data() == h_ac.data(), false);
    }
    {
      auto h_a = Kokkos::create_mirror_view(a);
      ASSERT_EQ(a.extents(), h_a.extents());
      ASSERT_EQ(a.data() == h_a.data(), host_accessible);
    }
    // There is an inconsistency here where for constant element type
    // create_mirror_view without space arg always returns a new allocation
    // while if you give a space arg it doesn't.
    {
      typename new_view_t::const_type ac = a;
      auto h_ac                          = Kokkos::create_mirror_view(ac);
      ASSERT_EQ(a.extents(), h_ac.extents());
      ASSERT_EQ(a.data() == h_ac.data(), false);
    }
    {
      auto h_a = Kokkos::create_mirror_view(Kokkos::HostSpace(), a);
      ASSERT_EQ(a.extents(), h_a.extents());
      ASSERT_EQ(a.data() == h_a.data(), is_host_space);
    }
    {
      typename new_view_t::const_type ac = a;
      auto h_ac = Kokkos::create_mirror_view(Kokkos::HostSpace(), ac);
      ASSERT_EQ(a.extents(), h_ac.extents());
      ASSERT_EQ(a.data() == h_ac.data(), is_host_space);
    }
    {
      init_view(a, 1, extents...);
      auto h_a = Kokkos::create_mirror_view_and_copy(a);
      ASSERT_EQ(a.extents(), h_a.extents());
      size_t errors = check_view(h_a, 1, extents...);
      ASSERT_EQ(errors, 0lu);
      ASSERT_EQ(a.data() == h_a.data(), host_accessible);
    }
    {
      init_view(a, 1, extents...);
      typename new_view_t::const_type ac = a;
      auto h_ac = Kokkos::create_mirror_view_and_copy(ac);
      ASSERT_EQ(a.extents(), h_ac.extents());
      size_t errors = check_view(h_ac, 1, extents...);
      ASSERT_EQ(errors, 0lu);
      ASSERT_EQ(a.data() == h_ac.data(), host_accessible);
    }
    {
      init_view(a, 1, extents...);
      auto h_a = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), a);
      ASSERT_EQ(a.extents(), h_a.extents());
      size_t errors = check_view(h_a, 1, extents...);
      ASSERT_EQ(errors, 0lu);
      ASSERT_EQ(a.data() == h_a.data(), is_host_space);
    }
    {
      init_view(a, 1, extents...);
      typename new_view_t::const_type ac = a;
      auto h_ac = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), ac);
      ASSERT_EQ(a.extents(), h_ac.extents());
      size_t errors = check_view(h_ac, 1, extents...);
      ASSERT_EQ(errors, 0lu);
      ASSERT_EQ(a.data() == h_ac.data(), is_host_space);
    }
    // Sanity check that the classical View args show
    // the same inconsistency in the return value for create_mirror_view
    // with respect to constant element types.
    {
      Kokkos::View<int, Kokkos::HostSpace> b("A");
      Kokkos::View<const int, Kokkos::HostSpace> bc(b);
      {
        auto h_b  = create_mirror_view(b);
        auto h_bc = create_mirror_view(bc);
        ASSERT_EQ(b.data(), h_b.data());
        ASSERT_NE(b.data(), h_bc.data());
      }
      {
        auto h_b  = create_mirror_view(Kokkos::HostSpace(), b);
        auto h_bc = create_mirror_view(Kokkos::HostSpace(), bc);
        ASSERT_EQ(b.data(), h_b.data());
        ASSERT_EQ(b.data(), h_bc.data());
      }
      {
        auto h_b  = create_mirror_view_and_copy(b);
        auto h_bc = create_mirror_view_and_copy(bc);
        ASSERT_EQ(b.data(), h_b.data());
        ASSERT_EQ(b.data(), h_bc.data());
      }
      {
        auto h_b  = create_mirror_view_and_copy(Kokkos::HostSpace(), b);
        auto h_bc = create_mirror_view_and_copy(Kokkos::HostSpace(), bc);
        ASSERT_EQ(b.data(), h_b.data());
        ASSERT_EQ(b.data(), h_bc.data());
      }
    }
  }
};

TEST(TEST_CATEGORY, view_minimal_mdspan_args_access) {
  ViewTestHarness<int, Kokkos::extents<unsigned, 3, 7>>::access(3, 7);
  ViewTestHarness<float, Kokkos::extents<unsigned, Kokkos::dynamic_extent,
                                         Kokkos::dynamic_extent, 7>>::access(3,
                                                                             5,
                                                                             7);
  ViewTestHarness<int, Kokkos::dextents<size_t, 6>,
                  Kokkos::Experimental::layout_right_padded<
                      Kokkos::dynamic_extent>>::access(3, 5lu, 7, 9, 2, 1);
}

TEST(TEST_CATEGORY, view_minimal_mdspan_args_deep_copy) {
  ViewTestHarness<int, Kokkos::extents<unsigned, 3, 7>>::deep_copy(3, 7);
  ViewTestHarness<float,
                  Kokkos::extents<unsigned, Kokkos::dynamic_extent,
                                  Kokkos::dynamic_extent, 7>>::deep_copy(3, 5,
                                                                         7);
  ViewTestHarness<int, Kokkos::dextents<size_t, 6>,
                  Kokkos::Experimental::layout_right_padded<
                      Kokkos::dynamic_extent>>::deep_copy(3, 5lu, 7, 9, 2, 1);
}

TEST(TEST_CATEGORY, view_minimal_mdspan_args_create_mirror) {
  ViewTestHarness<int, Kokkos::extents<unsigned, 3, 7>>::create_mirror(3, 7);
  ViewTestHarness<float,
                  Kokkos::extents<unsigned, Kokkos::dynamic_extent,
                                  Kokkos::dynamic_extent, 7>>::create_mirror(3,
                                                                             5,
                                                                             7);
  ViewTestHarness<int, Kokkos::dextents<size_t, 6>,
                  Kokkos::Experimental::layout_right_padded<
                      Kokkos::dynamic_extent>>::create_mirror(3, 5lu, 7, 9, 2,
                                                              1);
}
