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
#include <TestSYCL_Category.hpp>

namespace Test {

TEST(sycl, space_access) {
  static_assert(Kokkos::Impl::MemorySpaceAccess<Kokkos::HostSpace,
                                                Kokkos::HostSpace>::assignable);

  static_assert(
      Kokkos::Impl::MemorySpaceAccess<Kokkos::HostSpace,
                                      Kokkos::SYCLHostUSMSpace>::assignable);

  static_assert(
      !Kokkos::Impl::MemorySpaceAccess<Kokkos::HostSpace,
                                       Kokkos::SYCLDeviceUSMSpace>::assignable);

  static_assert(
      !Kokkos::Impl::MemorySpaceAccess<Kokkos::HostSpace,
                                       Kokkos::SYCLDeviceUSMSpace>::accessible);

  static_assert(
      !Kokkos::Impl::MemorySpaceAccess<Kokkos::HostSpace,
                                       Kokkos::SYCLSharedUSMSpace>::assignable);

  static_assert(
      Kokkos::Impl::MemorySpaceAccess<Kokkos::HostSpace,
                                      Kokkos::SYCLSharedUSMSpace>::accessible);

  //--------------------------------------

  static_assert(
      Kokkos::Impl::MemorySpaceAccess<Kokkos::SYCLDeviceUSMSpace,
                                      Kokkos::SYCLDeviceUSMSpace>::assignable);

  static_assert(
      Kokkos::Impl::MemorySpaceAccess<Kokkos::SYCLDeviceUSMSpace,
                                      Kokkos::SYCLSharedUSMSpace>::assignable);

  static_assert(
      !Kokkos::Impl::MemorySpaceAccess<Kokkos::SYCLDeviceUSMSpace,
                                       Kokkos::SYCLHostUSMSpace>::assignable);

  static_assert(
      Kokkos::Impl::MemorySpaceAccess<Kokkos::SYCLDeviceUSMSpace,
                                      Kokkos::SYCLHostUSMSpace>::accessible);

  static_assert(
      !Kokkos::Impl::MemorySpaceAccess<Kokkos::SYCLDeviceUSMSpace,
                                       Kokkos::HostSpace>::assignable);

  static_assert(
      !Kokkos::Impl::MemorySpaceAccess<Kokkos::SYCLDeviceUSMSpace,
                                       Kokkos::HostSpace>::accessible);

  //--------------------------------------

  static_assert(
      Kokkos::Impl::MemorySpaceAccess<Kokkos::SYCLSharedUSMSpace,
                                      Kokkos::SYCLSharedUSMSpace>::assignable);

  static_assert(
      !Kokkos::Impl::MemorySpaceAccess<Kokkos::SYCLSharedUSMSpace,
                                       Kokkos::SYCLDeviceUSMSpace>::assignable);

  static_assert(
      Kokkos::Impl::MemorySpaceAccess<Kokkos::SYCLSharedUSMSpace,
                                      Kokkos::SYCLDeviceUSMSpace>::accessible);

  static_assert(
      !Kokkos::Impl::MemorySpaceAccess<Kokkos::SYCLSharedUSMSpace,
                                       Kokkos::HostSpace>::assignable);

  static_assert(
      !Kokkos::Impl::MemorySpaceAccess<Kokkos::SYCLSharedUSMSpace,
                                       Kokkos::HostSpace>::accessible);

  static_assert(
      !Kokkos::Impl::MemorySpaceAccess<Kokkos::SYCLSharedUSMSpace,
                                       Kokkos::SYCLHostUSMSpace>::assignable);

  static_assert(
      Kokkos::Impl::MemorySpaceAccess<Kokkos::SYCLSharedUSMSpace,
                                      Kokkos::SYCLHostUSMSpace>::accessible);

  //--------------------------------------

  static_assert(
      Kokkos::Impl::MemorySpaceAccess<Kokkos::SYCLHostUSMSpace,
                                      Kokkos::SYCLHostUSMSpace>::assignable);

  static_assert(
      !Kokkos::Impl::MemorySpaceAccess<Kokkos::SYCLHostUSMSpace,
                                       Kokkos::HostSpace>::assignable);

  static_assert(Kokkos::Impl::MemorySpaceAccess<Kokkos::SYCLHostUSMSpace,
                                                Kokkos::HostSpace>::accessible);

  static_assert(
      !Kokkos::Impl::MemorySpaceAccess<Kokkos::SYCLHostUSMSpace,
                                       Kokkos::SYCLDeviceUSMSpace>::assignable);

  static_assert(
      !Kokkos::Impl::MemorySpaceAccess<Kokkos::SYCLHostUSMSpace,
                                       Kokkos::SYCLDeviceUSMSpace>::accessible);

  static_assert(
      !Kokkos::Impl::MemorySpaceAccess<Kokkos::SYCLHostUSMSpace,
                                       Kokkos::SYCLSharedUSMSpace>::assignable);

  static_assert(
      Kokkos::Impl::MemorySpaceAccess<Kokkos::SYCLHostUSMSpace,
                                      Kokkos::SYCLSharedUSMSpace>::accessible);

  //--------------------------------------

  static_assert(
      !Kokkos::SpaceAccessibility<Kokkos::SYCL, Kokkos::HostSpace>::accessible);

  static_assert(
      Kokkos::SpaceAccessibility<Kokkos::SYCL,
                                 Kokkos::SYCLDeviceUSMSpace>::accessible);

  static_assert(
      Kokkos::SpaceAccessibility<Kokkos::SYCL,
                                 Kokkos::SYCLSharedUSMSpace>::accessible);

  static_assert(
      Kokkos::SpaceAccessibility<Kokkos::SYCL,
                                 Kokkos::SYCLHostUSMSpace>::accessible);

  static_assert(
      !Kokkos::SpaceAccessibility<Kokkos::HostSpace,
                                  Kokkos::SYCLDeviceUSMSpace>::accessible);

  static_assert(
      Kokkos::SpaceAccessibility<Kokkos::HostSpace,
                                 Kokkos::SYCLSharedUSMSpace>::accessible);

  static_assert(
      Kokkos::SpaceAccessibility<Kokkos::HostSpace,
                                 Kokkos::SYCLHostUSMSpace>::accessible);

  static_assert(std::is_same_v<
                Kokkos::Impl::HostMirror<Kokkos::SYCLDeviceUSMSpace>::Space,
                Kokkos::HostSpace>);

  static_assert(std::is_same_v<
                Kokkos::Impl::HostMirror<Kokkos::SYCLSharedUSMSpace>::Space,
                Kokkos::Device<Kokkos::HostSpace::execution_space,
                               Kokkos::SYCLSharedUSMSpace>>);

  static_assert(Kokkos::Impl::MemorySpaceAccess<Kokkos::SYCLHostUSMSpace,
                                                Kokkos::HostSpace>::accessible);

  static_assert(
      Kokkos::Impl::MemorySpaceAccess<Kokkos::HostSpace,
                                      Kokkos::SYCLHostUSMSpace>::accessible);

  static_assert(
      std::is_same_v<Kokkos::Impl::HostMirror<Kokkos::SYCLHostUSMSpace>::Space,
                     Kokkos::SYCLHostUSMSpace>);

  static_assert(
      std::is_same_v<Kokkos::Device<Kokkos::HostSpace::execution_space,
                                    Kokkos::SYCLSharedUSMSpace>,
                     Kokkos::Device<Kokkos::HostSpace::execution_space,
                                    Kokkos::SYCLSharedUSMSpace>>);

  static_assert(
      Kokkos::SpaceAccessibility<Kokkos::Impl::HostMirror<Kokkos::SYCL>::Space,
                                 Kokkos::HostSpace>::accessible);

  static_assert(Kokkos::SpaceAccessibility<
                Kokkos::Impl::HostMirror<Kokkos::SYCLDeviceUSMSpace>::Space,
                Kokkos::HostSpace>::accessible);

  static_assert(Kokkos::SpaceAccessibility<
                Kokkos::Impl::HostMirror<Kokkos::SYCLSharedUSMSpace>::Space,
                Kokkos::HostSpace>::accessible);

  static_assert(Kokkos::SpaceAccessibility<
                Kokkos::Impl::HostMirror<Kokkos::SYCLHostUSMSpace>::Space,
                Kokkos::HostSpace>::accessible);
}

TEST(sycl, uvm) {
  int *uvm_ptr =
      static_cast<int *>(Kokkos::kokkos_malloc<Kokkos::SYCLSharedUSMSpace>(
          "uvm_ptr", sizeof(int)));

  *uvm_ptr = 42;

  Kokkos::SYCL().fence();
  Kokkos::parallel_for(
      Kokkos::RangePolicy<Kokkos::SYCL>(0, 1), KOKKOS_LAMBDA(int) {
        if (*uvm_ptr == 42) {
          *uvm_ptr = 2 * 42;
        }
      });
  Kokkos::SYCL().fence();

  EXPECT_EQ(*uvm_ptr, int(2 * 42));

  Kokkos::kokkos_free<Kokkos::SYCLSharedUSMSpace>(uvm_ptr);
}

template <class MemSpace, class ExecSpace>
struct TestViewSYCLAccessible {
  enum { N = 1000 };

  using V = Kokkos::View<double *, MemSpace>;

  V m_base;

  struct TagInit {};
  struct TagTest {};

  KOKKOS_INLINE_FUNCTION
  void operator()(const TagInit &, const int i) const { m_base[i] = i + 1; }

  KOKKOS_INLINE_FUNCTION
  void operator()(const TagTest &, const int i, long &error_count) const {
    if (m_base[i] != i + 1) ++error_count;
  }

  TestViewSYCLAccessible() : m_base("base", N) {}

  static void run() {
    TestViewSYCLAccessible self;
    Kokkos::parallel_for(
        Kokkos::RangePolicy<typename MemSpace::execution_space, TagInit>(0, N),
        self);
    typename MemSpace::execution_space().fence();

    // Next access is a different execution space, must complete prior kernel.
    long error_count = -1;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<ExecSpace, TagTest>(0, N), self,
                            error_count);
    EXPECT_EQ(error_count, 0);
  }
};

TEST(sycl, impl_view_accessible) {
  TestViewSYCLAccessible<Kokkos::SYCLDeviceUSMSpace, Kokkos::SYCL>::run();

  TestViewSYCLAccessible<Kokkos::SYCLSharedUSMSpace, Kokkos::SYCL>::run();
  TestViewSYCLAccessible<Kokkos::SYCLSharedUSMSpace,
                         Kokkos::HostSpace::execution_space>::run();

  TestViewSYCLAccessible<Kokkos::SYCLHostUSMSpace, Kokkos::SYCL>::run();
  TestViewSYCLAccessible<Kokkos::SYCLHostUSMSpace,
                         Kokkos::HostSpace::execution_space>::run();
}

}  // namespace Test
