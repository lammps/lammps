// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
// SPDX-FileCopyrightText: Copyright Contributors to the Kokkos project

#include <gtest/gtest.h>

#include <Kokkos_Macros.hpp>
#ifdef KOKKOS_ENABLE_EXPERIMENTAL_CXX20_MODULES
import kokkos.core;
#else
#include <Kokkos_Core.hpp>
#endif

namespace {

template <class T>
struct TestEmptyViewRuntimeUnmanaged {
  template <class ExecutionSpace>
  TestEmptyViewRuntimeUnmanaged(ExecutionSpace const& exec) {
    Kokkos::parallel_for(Kokkos::RangePolicy(exec, 0, 1), *this);
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
