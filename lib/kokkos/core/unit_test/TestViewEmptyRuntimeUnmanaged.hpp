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

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

namespace {

template <class T>
struct TestEmptyViewRuntimeUnmanaged {
  template <class ExecutionSpace>
  TestEmptyViewRuntimeUnmanaged(ExecutionSpace const& exec) {
    Kokkos::parallel_for(Kokkos::RangePolicy<ExecutionSpace>(exec, 0, 1),
                         *this);
  }

  KOKKOS_FUNCTION void operator()(int) const {
    T d{};
#if defined(_MSC_VER)
#pragma warning(push)
#pragma warning(disable : 4312)
#endif
    auto* p = reinterpret_cast<T*>(0xABADBABE);
#if defined(_MSC_VER)
#pragma warning(push)
#endif

    (void)Kokkos::View<T*>(p, 0);
    (void)Kokkos::View<T*>(&d, 0);
    (void)Kokkos::View<T*>(nullptr, 0);
    (void)Kokkos::View<T*>(NULL, 0);  // NOLINT(modernize-use-nullptr)
    (void)Kokkos::View<T*>(0, 0);     // NOLINT(modernize-use-nullptr)

    (void)Kokkos::View<T**>(p, 0, 0);
    (void)Kokkos::View<T**>(&d, 0, 0);
    (void)Kokkos::View<T**>(nullptr, 0, 0);
    (void)Kokkos::View<T**>(NULL, 0, 0);  // NOLINT(modernize-use-nullptr)
    (void)Kokkos::View<T**>(0, 0, 0);     // NOLINT(modernize-use-nullptr)
  }
};

TEST(TEST_CATEGORY, view_empty_runtime_unmanaged) {
  TEST_EXECSPACE exec;
  (void)TestEmptyViewRuntimeUnmanaged<float>(exec);
  (void)TestEmptyViewRuntimeUnmanaged<const double>(exec);
  (void)TestEmptyViewRuntimeUnmanaged<int>(exec);
  (void)TestEmptyViewRuntimeUnmanaged<char>(exec);
  (void)TestEmptyViewRuntimeUnmanaged<const char>(exec);
}

}  // namespace
