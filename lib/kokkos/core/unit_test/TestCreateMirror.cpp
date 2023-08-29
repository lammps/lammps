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

template <typename TestView, typename MemorySpace>
void check_memory_space(TestView, MemorySpace) {
  static_assert(std::is_same_v<typename TestView::memory_space, MemorySpace>);
}

template <class View>
auto host_mirror_test_space(View) {
  return std::conditional_t<
      Kokkos::SpaceAccessibility<Kokkos::HostSpace,
                                 typename View::memory_space>::accessible,
      typename View::memory_space, Kokkos::HostSpace>{};
}

template <typename View>
void test_create_mirror_properties(const View& view) {
  using namespace Kokkos;
  using DeviceMemorySpace = typename DefaultExecutionSpace::memory_space;

  // clang-format off
  
  // create_mirror
  check_memory_space(create_mirror(WithoutInitializing,                          view), host_mirror_test_space(view));
  check_memory_space(create_mirror(                                              view), host_mirror_test_space(view));
  check_memory_space(create_mirror(WithoutInitializing, DefaultExecutionSpace{}, view), DeviceMemorySpace{});
  check_memory_space(create_mirror(                     DefaultExecutionSpace{}, view), DeviceMemorySpace{});

  // create_mirror_view
  check_memory_space(create_mirror_view(WithoutInitializing,                          view), host_mirror_test_space(view));
  check_memory_space(create_mirror_view(                                              view), host_mirror_test_space(view));
  check_memory_space(create_mirror_view(WithoutInitializing, DefaultExecutionSpace{}, view), DeviceMemorySpace{});
  check_memory_space(create_mirror_view(                     DefaultExecutionSpace{}, view), DeviceMemorySpace{});

  // create_mirror view_alloc
  check_memory_space(create_mirror(view_alloc(WithoutInitializing),                      view), host_mirror_test_space(view));
  check_memory_space(create_mirror(view_alloc(),                                         view), host_mirror_test_space(view));
  check_memory_space(create_mirror(view_alloc(WithoutInitializing, DeviceMemorySpace{}), view), DeviceMemorySpace{});
  check_memory_space(create_mirror(view_alloc(                     DeviceMemorySpace{}), view), DeviceMemorySpace{});

  // create_mirror_view view_alloc
  check_memory_space(create_mirror_view(view_alloc(WithoutInitializing),                      view), host_mirror_test_space(view));
  check_memory_space(create_mirror_view(view_alloc(),                                         view), host_mirror_test_space(view));
  check_memory_space(create_mirror_view(view_alloc(WithoutInitializing, DeviceMemorySpace{}), view), DeviceMemorySpace{});
  check_memory_space(create_mirror_view(view_alloc(                     DeviceMemorySpace{}), view), DeviceMemorySpace{});

  // create_mirror view_alloc + execution space
  check_memory_space(create_mirror(view_alloc(DefaultHostExecutionSpace{}, WithoutInitializing),                      view), host_mirror_test_space(view));
  check_memory_space(create_mirror(view_alloc(DefaultHostExecutionSpace{}),                                           view), host_mirror_test_space(view));
  check_memory_space(create_mirror(view_alloc(DefaultExecutionSpace{},     WithoutInitializing, DeviceMemorySpace{}), view), DeviceMemorySpace{});
  check_memory_space(create_mirror(view_alloc(DefaultExecutionSpace{},                          DeviceMemorySpace{}), view), DeviceMemorySpace{});

  // create_mirror_view view_alloc + execution space
  check_memory_space(create_mirror_view(view_alloc(DefaultHostExecutionSpace{}, WithoutInitializing),                      view), host_mirror_test_space(view));
  check_memory_space(create_mirror_view(view_alloc(DefaultHostExecutionSpace{}),                                           view), host_mirror_test_space(view));
  check_memory_space(create_mirror_view(view_alloc(DefaultExecutionSpace{},     WithoutInitializing, DeviceMemorySpace{}), view), DeviceMemorySpace{});
  check_memory_space(create_mirror_view(view_alloc(DefaultExecutionSpace{},                          DeviceMemorySpace{}), view), DeviceMemorySpace{});

  // create_mirror_view_and_copy
  check_memory_space(create_mirror_view_and_copy(HostSpace{},         view), HostSpace{});
  check_memory_space(create_mirror_view_and_copy(DeviceMemorySpace{}, view), DeviceMemorySpace{});

  // create_mirror_view_and_copy view_alloc
  check_memory_space(create_mirror_view_and_copy(view_alloc(HostSpace{}),         view), HostSpace{});
  check_memory_space(create_mirror_view_and_copy(view_alloc(DeviceMemorySpace{}), view), DeviceMemorySpace{});

  // create_mirror_view_and_copy view_alloc + execution space
  check_memory_space(create_mirror_view_and_copy(view_alloc(HostSpace{},         DefaultHostExecutionSpace{}),   view), HostSpace{});
  check_memory_space(create_mirror_view_and_copy(view_alloc(DeviceMemorySpace{}, DefaultExecutionSpace{}),       view), DeviceMemorySpace{});

  // clang-format on
}

void test() {
  Kokkos::View<int*, Kokkos::DefaultExecutionSpace> device_view("device view",
                                                                10);
  Kokkos::View<int*, Kokkos::HostSpace> host_view("host view", 10);

  test_create_mirror_properties(device_view);
  test_create_mirror_properties(host_view);
}
