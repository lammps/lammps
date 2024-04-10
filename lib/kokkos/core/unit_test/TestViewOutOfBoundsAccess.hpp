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
#include <sstream>

#include <gtest/gtest.h>

namespace {

TEST(TEST_CATEGORY, append_formatted_multidimensional_index) {
  using Kokkos::Impl::append_formatted_multidimensional_index;
  {
    char buffer[64] = "my prefix ";
    append_formatted_multidimensional_index(buffer, 1);
    EXPECT_STREQ(buffer, "my prefix [1]");
  }
  {
    char buffer[64] = "I was here";
    append_formatted_multidimensional_index(buffer, 1, 2, 3);
    EXPECT_STREQ(buffer, "I was here[1,2,3]");
  }
  {
    char buffer[64] = "with mixed integer types ";
    append_formatted_multidimensional_index(buffer, 1u, -2);
    EXPECT_STREQ(buffer, "with mixed integer types [1,-2]");
  }
}

#ifdef KOKKOS_ENABLE_DEBUG_BOUNDS_CHECK

template <class View, class ExecutionSpace>
struct TestViewOutOfBoundAccess {
  View v;
  static constexpr auto rank = View::rank;

  template <std::size_t... Is>
  KOKKOS_FUNCTION decltype(auto) bad_access(std::index_sequence<Is...>) const {
    return v((Is * 1 + Is == 0 ? v.extent(Is) + 3 : 0)...);
  }

  KOKKOS_FUNCTION void operator()(int) const {
    ++bad_access(std::make_index_sequence<rank>{});
  }

  template <std::size_t... Is>
  std::string get_details(std::index_sequence<Is...>) {
    std::stringstream ss;
    ss << "with indices \\[";
    ((ss << (Is == 0 ? v.extent(Is) + 3 : 0)
         << (Is == View::rank() - 1 ? "\\]" : ",")),
     ...);
    ss << " but extents \\[";
    ((ss << v.extent(Is) << (Is == View::rank() - 1 ? "\\]" : ",")), ...);
    return ss.str();
  }

  auto get_details() {
    return get_details(std::make_index_sequence<View::rank()>());
  }

  TestViewOutOfBoundAccess(View w, ExecutionSpace const& s, std::string matcher)
      : v(std::move(w)) {
    constexpr bool view_accessible_from_execution_space =
        Kokkos::SpaceAccessibility<
            /*AccessSpace=*/ExecutionSpace,
            /*MemorySpace=*/typename View::memory_space>::accessible;
    EXPECT_TRUE(view_accessible_from_execution_space);

    matcher += ".*" + get_details();

    EXPECT_DEATH(
        {
          Kokkos::parallel_for(Kokkos::RangePolicy<ExecutionSpace>(s, 0, 1),
                               *this);
          Kokkos::fence();
        },
        matcher);
  }
};

template <class View, class LblOrPtr, std::size_t... Is>
auto make_view_impl(LblOrPtr x, std::index_sequence<Is...>) {
  return View(x, (Is + 1)...);
}

template <class View, class LblOrPtr>
auto make_view(LblOrPtr x) {
  return make_view_impl<View>(std::move(x),
                              std::make_index_sequence<View::rank>());
}

template <class ExecutionSpace>
void test_view_out_of_bounds_access() {
  ExecutionSpace const exec_space{};
  // clang-format off
  using V1 = Kokkos::View<int*,        ExecutionSpace>;
  using V2 = Kokkos::View<int**,       ExecutionSpace>;
  using V3 = Kokkos::View<int***,      ExecutionSpace>;
  using V4 = Kokkos::View<int****,     ExecutionSpace>;
  using V5 = Kokkos::View<int*****,    ExecutionSpace>;
  using V6 = Kokkos::View<int******,   ExecutionSpace>;
  using V7 = Kokkos::View<int*******,  ExecutionSpace>;
  using V8 = Kokkos::View<int********, ExecutionSpace>;
  std::string const prefix = "Kokkos::View ERROR: out of bounds access";
  std::string const lbl = "my_label";
  TestViewOutOfBoundAccess(make_view<V1>(lbl), exec_space, prefix + ".*" + lbl);
  TestViewOutOfBoundAccess(make_view<V2>(lbl), exec_space, prefix + ".*" + lbl);
  TestViewOutOfBoundAccess(make_view<V3>(lbl), exec_space, prefix + ".*" + lbl);
  TestViewOutOfBoundAccess(make_view<V4>(lbl), exec_space, prefix + ".*" + lbl);
  TestViewOutOfBoundAccess(make_view<V5>(lbl), exec_space, prefix + ".*" + lbl);
  TestViewOutOfBoundAccess(make_view<V6>(lbl), exec_space, prefix + ".*" + lbl);
  TestViewOutOfBoundAccess(make_view<V7>(lbl), exec_space, prefix + ".*" + lbl);
  TestViewOutOfBoundAccess(make_view<V8>(lbl), exec_space, prefix + ".*" + lbl);
  int* const ptr = nullptr;
  TestViewOutOfBoundAccess(make_view<V1>(ptr), exec_space, prefix + ".*UNMANAGED");
  TestViewOutOfBoundAccess(make_view<V2>(ptr), exec_space, prefix + ".*UNMANAGED");
  TestViewOutOfBoundAccess(make_view<V3>(ptr), exec_space, prefix + ".*UNMANAGED");
  TestViewOutOfBoundAccess(make_view<V4>(ptr), exec_space, prefix + ".*UNMANAGED");
  TestViewOutOfBoundAccess(make_view<V5>(ptr), exec_space, prefix + ".*UNMANAGED");
  TestViewOutOfBoundAccess(make_view<V6>(ptr), exec_space, prefix + ".*UNMANAGED");
  TestViewOutOfBoundAccess(make_view<V7>(ptr), exec_space, prefix + ".*UNMANAGED");
  TestViewOutOfBoundAccess(make_view<V8>(ptr), exec_space, prefix + ".*UNMANAGED");
  // clang-format on
}

TEST(TEST_CATEGORY_DEATH, view_out_of_bounds_access) {
  ::testing::FLAGS_gtest_death_test_style = "threadsafe";

  using ExecutionSpace = TEST_EXECSPACE;

  if (false && Kokkos::SpaceAccessibility<
                   /*AccessSpace=*/ExecutionSpace,
                   /*MemorySpace=*/Kokkos::HostSpace>::accessible) {
    GTEST_SKIP() << "skipping since no memory access violation would occur";
  }

#if defined(KOKKOS_ENABLE_SYCL) && defined(NDEBUG)  // FIXME_SYCL
  if (std::is_same_v<ExecutionSpace, Kokkos::Experimental::SYCL>) {
    GTEST_SKIP() << "skipping SYCL device-side abort does not work when NDEBUG "
                    "is defined";
  }
#endif
#if defined(KOKKOS_ENABLE_OPENMPTARGET)  // FIXME_OPENMPTARGET
  if (std::is_same_v<ExecutionSpace, Kokkos::Experimental::OpenMPTarget>) {
    GTEST_SKIP() << "skipping because OpenMPTarget backend is currently not "
                    "able to abort from the device";
  }
#endif
#if defined(KOKKOS_ENABLE_OPENACC)  // FIXME_OPENACC
  if (std::is_same<ExecutionSpace, Kokkos::Experimental::OpenACC>::value) {
    GTEST_SKIP() << "skipping because OpenACC backend is currently not "
                    "able to abort from the device";
  }
#endif

  test_view_out_of_bounds_access<ExecutionSpace>();
}

#endif

}  // namespace
