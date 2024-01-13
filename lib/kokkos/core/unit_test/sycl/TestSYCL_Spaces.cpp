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
                                                Kokkos::HostSpace>::assignable,
                "");

  static_assert(Kokkos::Impl::MemorySpaceAccess<
                    Kokkos::HostSpace,
                    Kokkos::Experimental::SYCLHostUSMSpace>::assignable,
                "");

  static_assert(!Kokkos::Impl::MemorySpaceAccess<
                    Kokkos::HostSpace,
                    Kokkos::Experimental::SYCLDeviceUSMSpace>::assignable,
                "");

  static_assert(!Kokkos::Impl::MemorySpaceAccess<
                    Kokkos::HostSpace,
                    Kokkos::Experimental::SYCLDeviceUSMSpace>::accessible,
                "");

  static_assert(!Kokkos::Impl::MemorySpaceAccess<
                    Kokkos::HostSpace,
                    Kokkos::Experimental::SYCLSharedUSMSpace>::assignable,
                "");

  static_assert(Kokkos::Impl::MemorySpaceAccess<
                    Kokkos::HostSpace,
                    Kokkos::Experimental::SYCLSharedUSMSpace>::accessible,
                "");

  //--------------------------------------

  static_assert(Kokkos::Impl::MemorySpaceAccess<
                    Kokkos::Experimental::SYCLDeviceUSMSpace,
                    Kokkos::Experimental::SYCLDeviceUSMSpace>::assignable,
                "");

  static_assert(Kokkos::Impl::MemorySpaceAccess<
                    Kokkos::Experimental::SYCLDeviceUSMSpace,
                    Kokkos::Experimental::SYCLSharedUSMSpace>::assignable,
                "");

  static_assert(!Kokkos::Impl::MemorySpaceAccess<
                    Kokkos::Experimental::SYCLDeviceUSMSpace,
                    Kokkos::Experimental::SYCLHostUSMSpace>::assignable,
                "");

  static_assert(Kokkos::Impl::MemorySpaceAccess<
                    Kokkos::Experimental::SYCLDeviceUSMSpace,
                    Kokkos::Experimental::SYCLHostUSMSpace>::accessible,
                "");

  static_assert(
      !Kokkos::Impl::MemorySpaceAccess<Kokkos::Experimental::SYCLDeviceUSMSpace,
                                       Kokkos::HostSpace>::assignable,
      "");

  static_assert(
      !Kokkos::Impl::MemorySpaceAccess<Kokkos::Experimental::SYCLDeviceUSMSpace,
                                       Kokkos::HostSpace>::accessible,
      "");

  //--------------------------------------

  static_assert(Kokkos::Impl::MemorySpaceAccess<
                    Kokkos::Experimental::SYCLSharedUSMSpace,
                    Kokkos::Experimental::SYCLSharedUSMSpace>::assignable,
                "");

  static_assert(!Kokkos::Impl::MemorySpaceAccess<
                    Kokkos::Experimental::SYCLSharedUSMSpace,
                    Kokkos::Experimental::SYCLDeviceUSMSpace>::assignable,
                "");

  static_assert(Kokkos::Impl::MemorySpaceAccess<
                    Kokkos::Experimental::SYCLSharedUSMSpace,
                    Kokkos::Experimental::SYCLDeviceUSMSpace>::accessible,
                "");

  static_assert(
      !Kokkos::Impl::MemorySpaceAccess<Kokkos::Experimental::SYCLSharedUSMSpace,
                                       Kokkos::HostSpace>::assignable,
      "");

  static_assert(
      !Kokkos::Impl::MemorySpaceAccess<Kokkos::Experimental::SYCLSharedUSMSpace,
                                       Kokkos::HostSpace>::accessible,
      "");

  static_assert(!Kokkos::Impl::MemorySpaceAccess<
                    Kokkos::Experimental::SYCLSharedUSMSpace,
                    Kokkos::Experimental::SYCLHostUSMSpace>::assignable,
                "");

  static_assert(Kokkos::Impl::MemorySpaceAccess<
                    Kokkos::Experimental::SYCLSharedUSMSpace,
                    Kokkos::Experimental::SYCLHostUSMSpace>::accessible,
                "");

  //--------------------------------------

  static_assert(Kokkos::Impl::MemorySpaceAccess<
                    Kokkos::Experimental::SYCLHostUSMSpace,
                    Kokkos::Experimental::SYCLHostUSMSpace>::assignable,
                "");

  static_assert(
      !Kokkos::Impl::MemorySpaceAccess<Kokkos::Experimental::SYCLHostUSMSpace,
                                       Kokkos::HostSpace>::assignable,
      "");

  static_assert(
      Kokkos::Impl::MemorySpaceAccess<Kokkos::Experimental::SYCLHostUSMSpace,
                                      Kokkos::HostSpace>::accessible,
      "");

  static_assert(!Kokkos::Impl::MemorySpaceAccess<
                    Kokkos::Experimental::SYCLHostUSMSpace,
                    Kokkos::Experimental::SYCLDeviceUSMSpace>::assignable,
                "");

  static_assert(!Kokkos::Impl::MemorySpaceAccess<
                    Kokkos::Experimental::SYCLHostUSMSpace,
                    Kokkos::Experimental::SYCLDeviceUSMSpace>::accessible,
                "");

  static_assert(!Kokkos::Impl::MemorySpaceAccess<
                    Kokkos::Experimental::SYCLHostUSMSpace,
                    Kokkos::Experimental::SYCLSharedUSMSpace>::assignable,
                "");

  static_assert(Kokkos::Impl::MemorySpaceAccess<
                    Kokkos::Experimental::SYCLHostUSMSpace,
                    Kokkos::Experimental::SYCLSharedUSMSpace>::accessible,
                "");

  //--------------------------------------

  static_assert(!Kokkos::SpaceAccessibility<Kokkos::Experimental::SYCL,
                                            Kokkos::HostSpace>::accessible,
                "");

  static_assert(Kokkos::SpaceAccessibility<
                    Kokkos::Experimental::SYCL,
                    Kokkos::Experimental::SYCLDeviceUSMSpace>::accessible,
                "");

  static_assert(Kokkos::SpaceAccessibility<
                    Kokkos::Experimental::SYCL,
                    Kokkos::Experimental::SYCLSharedUSMSpace>::accessible,
                "");

  static_assert(Kokkos::SpaceAccessibility<
                    Kokkos::Experimental::SYCL,
                    Kokkos::Experimental::SYCLHostUSMSpace>::accessible,
                "");

  static_assert(!Kokkos::SpaceAccessibility<
                    Kokkos::HostSpace,
                    Kokkos::Experimental::SYCLDeviceUSMSpace>::accessible,
                "");

  static_assert(Kokkos::SpaceAccessibility<
                    Kokkos::HostSpace,
                    Kokkos::Experimental::SYCLSharedUSMSpace>::accessible,
                "");

  static_assert(Kokkos::SpaceAccessibility<
                    Kokkos::HostSpace,
                    Kokkos::Experimental::SYCLHostUSMSpace>::accessible,
                "");

  static_assert(
      std::is_same<Kokkos::Impl::HostMirror<
                       Kokkos::Experimental::SYCLDeviceUSMSpace>::Space,
                   Kokkos::HostSpace>::value,
      "");

  static_assert(
      std::is_same<
          Kokkos::Impl::HostMirror<
              Kokkos::Experimental::SYCLSharedUSMSpace>::Space,
          Kokkos::Device<Kokkos::HostSpace::execution_space,
                         Kokkos::Experimental::SYCLSharedUSMSpace>>::value,
      "");

  static_assert(
      Kokkos::Impl::MemorySpaceAccess<Kokkos::Experimental::SYCLHostUSMSpace,
                                      Kokkos::HostSpace>::accessible,
      "");

  static_assert(Kokkos::Impl::MemorySpaceAccess<
                    Kokkos::HostSpace,
                    Kokkos::Experimental::SYCLHostUSMSpace>::accessible,
                "");

  static_assert(std::is_same<Kokkos::Impl::HostMirror<
                                 Kokkos::Experimental::SYCLHostUSMSpace>::Space,
                             Kokkos::Experimental::SYCLHostUSMSpace>::value,
                "");

  static_assert(
      std::is_same<
          Kokkos::Device<Kokkos::HostSpace::execution_space,
                         Kokkos::Experimental::SYCLSharedUSMSpace>,
          Kokkos::Device<Kokkos::HostSpace::execution_space,
                         Kokkos::Experimental::SYCLSharedUSMSpace>>::value,
      "");

  static_assert(Kokkos::SpaceAccessibility<
                    Kokkos::Impl::HostMirror<Kokkos::Experimental::SYCL>::Space,
                    Kokkos::HostSpace>::accessible,
                "");

  static_assert(Kokkos::SpaceAccessibility<
                    Kokkos::Impl::HostMirror<
                        Kokkos::Experimental::SYCLDeviceUSMSpace>::Space,
                    Kokkos::HostSpace>::accessible,
                "");

  static_assert(Kokkos::SpaceAccessibility<
                    Kokkos::Impl::HostMirror<
                        Kokkos::Experimental::SYCLSharedUSMSpace>::Space,
                    Kokkos::HostSpace>::accessible,
                "");

  static_assert(Kokkos::SpaceAccessibility<
                    Kokkos::Impl::HostMirror<
                        Kokkos::Experimental::SYCLHostUSMSpace>::Space,
                    Kokkos::HostSpace>::accessible,
                "");
}

TEST(sycl, uvm) {
  int *uvm_ptr = static_cast<int *>(
      Kokkos::kokkos_malloc<Kokkos::Experimental::SYCLSharedUSMSpace>(
          "uvm_ptr", sizeof(int)));

  *uvm_ptr = 42;

  Kokkos::Experimental::SYCL().fence();
  Kokkos::parallel_for(
      Kokkos::RangePolicy<Kokkos::Experimental::SYCL>(0, 1),
      KOKKOS_LAMBDA(int) {
        if (*uvm_ptr == 42) {
          *uvm_ptr = 2 * 42;
        }
      });
  Kokkos::Experimental::SYCL().fence();

  EXPECT_EQ(*uvm_ptr, int(2 * 42));

  Kokkos::kokkos_free<Kokkos::Experimental::SYCLSharedUSMSpace>(uvm_ptr);
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
  TestViewSYCLAccessible<Kokkos::Experimental::SYCLDeviceUSMSpace,
                         Kokkos::Experimental::SYCL>::run();

  TestViewSYCLAccessible<Kokkos::Experimental::SYCLSharedUSMSpace,
                         Kokkos::Experimental::SYCL>::run();
  TestViewSYCLAccessible<Kokkos::Experimental::SYCLSharedUSMSpace,
                         Kokkos::HostSpace::execution_space>::run();

  TestViewSYCLAccessible<Kokkos::Experimental::SYCLHostUSMSpace,
                         Kokkos::Experimental::SYCL>::run();
  TestViewSYCLAccessible<Kokkos::Experimental::SYCLHostUSMSpace,
                         Kokkos::HostSpace::execution_space>::run();
}

}  // namespace Test
